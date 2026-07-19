"""
PyTorch port of the DAM-IC DEC repo's scripts/DECLayer.py (the Student's-t
soft-assignment clustering layer) and the dec_vae_cluster/target_distribution
functions from scripts/clustering_functions.py.

Note on a fix vs. the original repo: clustering_functions.dec_vae_cluster()
is called as `dec_vae_cluster(encoder, ...)` where `encoder` is actually the
whole `xvae` class instance (see cluster_mlp_vae), then does `encoder.output`
/ `Model(inputs=encoder.input, ...)` - but `xvae` instances don't expose a
Keras `.input`/`.output` API, and even the intended `encoder.output` (the
3-tuple [z_mean, z_log_sigma, z] from xvae.encoder) isn't the 2D tensor
DECLayer requires. This looks like a leftover from refactoring. Here, DEC
operates explicitly on the deterministic latent code z_mean (what
xvae.Z_encoder was for in the original), which is the standard and sensible
choice for K-means initialisation and cluster-distance computation.
"""
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.cluster import KMeans

from models.metrics import calculate_metrics


class DECHead(nn.Module):
    """Clustering layer converting latent samples into soft labels via the
    Student's t-distribution, as in t-SNE. Direct port of DECLayer.py."""

    def __init__(self, n_clusters, latent_dim, alpha=1.0):
        super().__init__()
        self.n_clusters = n_clusters
        self.alpha = alpha
        # matches Keras add_weight(initializer='glorot_uniform')
        self.clusters = nn.Parameter(torch.empty(n_clusters, latent_dim))
        nn.init.xavier_uniform_(self.clusters)

    def forward(self, z):
        """z: (n_samples, latent_dim) -> q: (n_samples, n_clusters) soft cluster
        assignment probabilities, using a Student's t-distribution kernel
        (same kernel t-SNE uses) as the similarity measure between each
        sample's latent code and each cluster centroid."""
        dist_sq = torch.sum((z.unsqueeze(1) - self.clusters) ** 2, dim=2)
        q = 1.0 / (1.0 + dist_sq / self.alpha)
        q = q ** ((self.alpha + 1.0) / 2.0)
        q = q / q.sum(dim=1, keepdim=True)
        return q


def target_distribution(q):
    """Auxiliary target distribution P, computed from soft labels Q.

    DEC's self-training trick: square each probability (sharpening toward
    more confident assignments) and re-normalize, with a per-cluster-size
    correction (q.sum(axis=0)) so large clusters don't dominate. The model is
    then trained to make Q match this sharper P, iteratively pulling points
    closer to their assigned cluster's centroid.
    """
    weight = q ** 2 / q.sum(axis=0)
    return (weight.T / weight.sum(axis=1)).T


class XDEC:
    """Runs the X-DEC algorithm: fine-tunes a pretrained xvae encoder jointly
    with a DEC clustering head, on the mixed numerical (s1) + binary (s2)
    input. Port of clustering_functions.dec_vae_cluster / cluster_mlp_vae."""

    def __init__(self, xvae_model, n_clusters, alpha=1.0, seed=5192):
        self.xvae = xvae_model
        self.n_clusters = n_clusters
        self.seed = seed
        self.device = xvae_model.device
        self.head = DECHead(n_clusters, xvae_model.args.ls, alpha).to(self.device)

    def _encode(self, x1, x2):
        self.xvae.net.eval()
        with torch.no_grad():
            z_mean, _ = self.xvae.net.encode(x1, x2)
        return z_mean

    def fit(self, s1, s2, y=None, maxiter=2000, batch_size=64, update_interval=50,
            tol=0.01, lr=0.01, momentum=0.9, n_init=20, verbose=True):
        torch.manual_seed(self.seed)
        np.random.seed(self.seed)

        x1 = torch.tensor(np.asarray(s1), dtype=torch.float32, device=self.device)
        x2 = torch.tensor(np.asarray(s2), dtype=torch.float32, device=self.device)
        n = x1.shape[0]

        # K-means initialisation of the cluster centers in latent space: run
        # the already-pretrained xvae encoder once, K-means-cluster the
        # resulting latent codes, and use those cluster centers as the
        # DECHead's starting `clusters` parameter.
        z_init = self._encode(x1, x2).cpu().numpy()
        kmeans = KMeans(n_clusters=self.n_clusters, n_init=n_init, random_state=self.seed)
        y_pred = kmeans.fit_predict(z_init)
        y_pred_last = np.copy(y_pred)
        with torch.no_grad():
            self.head.clusters.copy_(
                torch.tensor(kmeans.cluster_centers_, dtype=torch.float32, device=self.device))

        # Both the xvae encoder AND the clustering head are trained jointly
        # from here on - the encoder itself gets fine-tuned to produce a more
        # clusterable latent space as training progresses, not just fit a
        # fixed pretrained embedding.
        params = list(self.xvae.net.parameters()) + list(self.head.parameters())
        optimizer = torch.optim.SGD(params, lr=lr, momentum=momentum)

        index_array = np.arange(n)
        index = 0
        p = None
        loss_val = 0.0

        for ite in range(int(maxiter)):
            # Every `update_interval` iterations: recompute soft assignments Q,
            # derive a sharper target P from them, and check whether cluster
            # assignments have stabilized enough to stop early.
            if ite % update_interval == 0:
                self.xvae.net.eval()
                with torch.no_grad():
                    z_mean, _ = self.xvae.net.encode(x1, x2)
                    q = self.head(z_mean).cpu().numpy()
                p = target_distribution(q)

                y_pred = q.argmax(1)
                if y is not None:
                    acc, ari, nmi = calculate_metrics(y, y_pred)
                    if verbose:
                        print('Iter {}: acc = {:.5f}, nmi = {:.5f}, ari = {:.5f} ; loss={:.5f}'.format(
                            ite, acc, nmi, ari, loss_val))

                delta_label = np.sum(y_pred != y_pred_last).astype(np.float32) / y_pred.shape[0]
                y_pred_last = np.copy(y_pred)
                if ite > 0 and delta_label < tol:
                    print('delta_label {} < tol {}'.format(delta_label, tol))
                    print('performed {} iterations.'.format(ite))
                    print('Reached tolerance threshold. Stopping training.')
                    break

            # Mini-batch gradient step: encode this batch (both branches), compute
            # its soft assignments Q, and push Q toward the (fixed until the next
            # update_interval refresh) target P via KL divergence.
            idx = index_array[index * batch_size: min((index + 1) * batch_size, n)]
            if len(idx) < 2:  # BatchNorm needs more than one sample
                index = index + 1 if (index + 1) * batch_size <= n else 0
                continue

            self.xvae.net.train()
            optimizer.zero_grad()
            z_mean, _ = self.xvae.net.encode(x1[idx], x2[idx])
            q_batch = self.head(z_mean)
            p_batch = torch.tensor(p[idx], dtype=torch.float32, device=self.device)
            # KL(p || q), matching Keras model.compile(loss='kld')
            loss = F.kl_div((q_batch + 1e-10).log(), p_batch, reduction='batchmean')
            loss.backward()
            optimizer.step()
            loss_val = loss.item()

            index = index + 1 if (index + 1) * batch_size <= n else 0

        self.xvae.net.eval()
        with torch.no_grad():
            z_mean, _ = self.xvae.net.encode(x1, x2)
            q = self.head(z_mean).cpu().numpy()
        y_pred = q.argmax(1)
        if y is not None:
            acc, ari, nmi = calculate_metrics(y, y_pred)
            print('Final: acc = {:.5f}, nmi = {:.5f}, ari = {:.5f} ; loss={:.5f}'.format(
                acc, nmi, ari, loss_val))

        centroids = self.head.clusters.detach().cpu().numpy()
        return y_pred, q, centroids, z_mean.cpu().numpy()
