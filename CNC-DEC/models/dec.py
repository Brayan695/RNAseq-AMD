"""
PyTorch port of the DAM-IC DEC repo's scripts/DECLayer.py (the Student's-t
soft-assignment clustering layer) and the dec_vae_cluster/target_distribution
functions from scripts/clustering_functions.py, adapted for CNC-DEC's single
concatenated-input CNCVAE encoder.

Note: `CNCVAE._CNCVAENet.encode(x)` returns `(z_mean, z_log_sigma)` - the
FIRST element is the deterministic latent code DEC operates on (unlike plain
DEC's `_AENet.encode()`, which returns `(h1, z)` with the code as the SECOND
element - don't mix the two up when reusing this pattern elsewhere).
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
        dist_sq = torch.sum((z.unsqueeze(1) - self.clusters) ** 2, dim=2)
        q = 1.0 / (1.0 + dist_sq / self.alpha)
        q = q ** ((self.alpha + 1.0) / 2.0)
        q = q / q.sum(dim=1, keepdim=True)
        return q


def target_distribution(q):
    """Auxiliary target distribution P, computed from soft labels Q."""
    weight = q ** 2 / q.sum(axis=0)
    return (weight.T / weight.sum(axis=1)).T


class DEC:
    """Runs the DEC algorithm: fine-tunes a pretrained CNCVAE encoder jointly
    with a DEC clustering head. Port of clustering_functions.dec_vae_cluster,
    but for the single concatenated-input CNCVAE rather than the dual-branch
    xvae (see X-DEC/models/dec.py:XDEC for that variant)."""

    def __init__(self, autoencoder, n_clusters, alpha=1.0, seed=5192):
        self.autoencoder = autoencoder
        self.n_clusters = n_clusters
        self.seed = seed
        self.device = autoencoder.device
        self.head = DECHead(n_clusters, autoencoder.args.ls, alpha).to(self.device)

    def _encode(self, x):
        self.autoencoder.net.eval()
        with torch.no_grad():
            z_mean, _ = self.autoencoder.net.encode(x)
        return z_mean

    def fit(self, s1, s2, y=None, maxiter=2000, batch_size=64, update_interval=50,
            tol=0.01, lr=0.01, momentum=0.9, n_init=20, verbose=True):
        torch.manual_seed(self.seed)
        np.random.seed(self.seed)

        x = torch.tensor(
            np.concatenate((np.asarray(s1), np.asarray(s2)), axis=-1),
            dtype=torch.float32, device=self.device)
        n = x.shape[0]

        # K-means initialisation of the cluster centers in latent space
        z_init = self._encode(x).cpu().numpy()
        kmeans = KMeans(n_clusters=self.n_clusters, n_init=n_init, random_state=self.seed)
        y_pred = kmeans.fit_predict(z_init)
        y_pred_last = np.copy(y_pred)
        with torch.no_grad():
            self.head.clusters.copy_(
                torch.tensor(kmeans.cluster_centers_, dtype=torch.float32, device=self.device))

        params = list(self.autoencoder.net.parameters()) + list(self.head.parameters())
        optimizer = torch.optim.SGD(params, lr=lr, momentum=momentum)

        index_array = np.arange(n)
        index = 0
        p = None
        loss_val = 0.0

        for ite in range(int(maxiter)):
            if ite % update_interval == 0:
                self.autoencoder.net.eval()
                with torch.no_grad():
                    z_mean, _ = self.autoencoder.net.encode(x)
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

            idx = index_array[index * batch_size: min((index + 1) * batch_size, n)]
            if len(idx) < 2:  # BatchNorm needs more than one sample
                index = index + 1 if (index + 1) * batch_size <= n else 0
                continue

            self.autoencoder.net.train()
            optimizer.zero_grad()
            z_mean, _ = self.autoencoder.net.encode(x[idx])
            q_batch = self.head(z_mean)
            p_batch = torch.tensor(p[idx], dtype=torch.float32, device=self.device)
            # KL(p || q), matching Keras model.compile(loss='kld')
            loss = F.kl_div((q_batch + 1e-10).log(), p_batch, reduction='batchmean')
            loss.backward()
            optimizer.step()
            loss_val = loss.item()

            index = index + 1 if (index + 1) * batch_size <= n else 0

        self.autoencoder.net.eval()
        with torch.no_grad():
            z_mean, _ = self.autoencoder.net.encode(x)
            q = self.head(z_mean).cpu().numpy()
        y_pred = q.argmax(1)
        if y is not None:
            acc, ari, nmi = calculate_metrics(y, y_pred)
            print('Final: acc = {:.5f}, nmi = {:.5f}, ari = {:.5f} ; loss={:.5f}'.format(
                acc, nmi, ari, loss_val))

        centroids = self.head.clusters.detach().cpu().numpy()
        return y_pred, q, centroids, z_mean.cpu().numpy()
