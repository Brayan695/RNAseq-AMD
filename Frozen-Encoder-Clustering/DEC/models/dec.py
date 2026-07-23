"""
PyTorch port of the DAM-IC DEC repo's scripts/DECLayer.py (the Student's-t
soft-assignment clustering layer) and the dec_cluster/target_distribution
functions from scripts/clustering_functions.py - the plain, single-input DEC
path (not dec_vae_cluster, which is the X-shaped/mixed-type variant and is
not used here).

Unlike dec_vae_cluster (which passes a whole `xvae` class instance where a
Keras Model with .input/.output is expected - a mismatch), dec_cluster's
`encoder` really is `Model(input_arr, encoded)`, i.e. a real single-input,
single-output model whose `.output` is exactly the 2D latent tensor DECLayer
expects. That's the well-formed code path this module ports.
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


class DEC:
    """FROZEN-ENCODER VARIANT (Frozen-Encoder-Clustering experiment): unlike the
    original DEC (which fine-tunes the MLP autoencoder's encoder jointly with
    the clustering head), this version clusters on top of an already-trained,
    already-loaded encoder WITHOUT changing its weights at all - only the
    DECHead's cluster centers are trained. This answers a different question
    than the original: does the pretrained autoencoder embedding, as-is,
    already cluster well, without DEC's joint encoder fine-tuning helping it
    along? Port of clustering_functions.dec_cluster / cluster_mlp_autoencoder."""

    def __init__(self, autoencoder, n_clusters, alpha=1.0, seed=5192):
        self.autoencoder = autoencoder
        self.n_clusters = n_clusters
        self.seed = seed
        self.device = autoencoder.device
        self.head = DECHead(n_clusters, autoencoder.args.neurons_e, alpha).to(self.device)

    def _encode(self, x):
        self.autoencoder.net.eval()
        with torch.no_grad():
            _, z = self.autoencoder.net.encode(x)
        return z

    def fit(self, X, y=None, maxiter=8000, batch_size=256, update_interval=140,
            tol=0.01, lr=0.01, momentum=0.9, n_init=20, verbose=True):
        torch.manual_seed(self.seed)
        np.random.seed(self.seed)

        x = torch.tensor(np.asarray(X), dtype=torch.float32, device=self.device)
        n = x.shape[0]

        # Freeze the encoder completely: no gradient updates, and it stays in
        # eval() mode for the rest of fit() (never flips to .train()) so even
        # BatchNorm's running statistics don't drift from these batches either.
        # This is the one change from the original DEC: everything below still
        # runs the real DEC self-training algorithm (K-means init, Q/P
        # soft-assignment refinement) - it just never touches the encoder.
        for p in self.autoencoder.net.parameters():
            p.requires_grad_(False)
        self.autoencoder.net.eval()

        # K-means initialisation of the cluster centers in latent space: run
        # the already-pretrained autoencoder's encoder once, K-means-cluster
        # the resulting latent codes, and use those cluster centers as the
        # DECHead's starting `clusters` parameter.
        z_init = self._encode(x).cpu().numpy()
        kmeans = KMeans(n_clusters=self.n_clusters, n_init=n_init, random_state=self.seed)
        y_pred = kmeans.fit_predict(z_init)
        y_pred_last = np.copy(y_pred)
        with torch.no_grad():
            self.head.clusters.copy_(
                torch.tensor(kmeans.cluster_centers_, dtype=torch.float32, device=self.device))

        # Only the clustering head is trained from here on - the encoder's
        # parameters were frozen above and are deliberately left out of
        # `params`, so optimizer.step() never touches them.
        params = list(self.head.parameters())
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
                self.autoencoder.net.eval()
                with torch.no_grad():
                    _, z = self.autoencoder.net.encode(x)
                    q = self.head(z).cpu().numpy()
                p = target_distribution(q)

                y_pred = q.argmax(1)
                if y is not None:
                    acc, ari, nmi = calculate_metrics(y, y_pred)
                    if verbose:
                        print('Iter {}: acc = {:.5f}, nmi = {:.5f}, ari = {:.5f} ; loss={:.5f}'.format(
                            ite, acc, nmi, ari, loss_val))

                delta_label = np.mean(y_pred != y_pred_last).astype(np.float32)
                y_pred_last = np.copy(y_pred)
                if ite > 0 and delta_label < tol:
                    print('delta_label {} < tol {}'.format(delta_label, tol))
                    print('performed {} iterations.'.format(ite))
                    print('Reached tolerance threshold. Stopping training.')
                    break

            # Mini-batch gradient step: encode this batch, compute its soft
            # assignments Q, and push Q toward the (fixed until the next
            # update_interval refresh) target P via KL divergence.
            idx = index_array[index * batch_size: min((index + 1) * batch_size, n)]
            if len(idx) < 1:
                index = index + 1 if (index + 1) * batch_size <= n else 0
                continue

            optimizer.zero_grad()
            _, z_batch = self.autoencoder.net.encode(x[idx])
            q_batch = self.head(z_batch)
            p_batch = torch.tensor(p[idx], dtype=torch.float32, device=self.device)
            # KL(p || q), matching Keras model.compile(loss='kld')
            loss = F.kl_div((q_batch + 1e-10).log(), p_batch, reduction='batchmean')
            loss.backward()
            optimizer.step()
            loss_val = loss.item()

            index = index + 1 if (index + 1) * batch_size <= n else 0  # advance to next mini-batch, wrapping around

        self.autoencoder.net.eval()
        with torch.no_grad():
            _, z = self.autoencoder.net.encode(x)
            q = self.head(z).cpu().numpy()
        y_pred = q.argmax(1)
        if y is not None:
            acc, ari, nmi = calculate_metrics(y, y_pred)
            print('Final: acc = {:.5f}, nmi = {:.5f}, ari = {:.5f} ; loss={:.5f}'.format(
                acc, nmi, ari, loss_val))

        centroids = self.head.clusters.detach().cpu().numpy()
        return y_pred, q, centroids, z.cpu().numpy()
