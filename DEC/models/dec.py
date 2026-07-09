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
    """Runs the DEC algorithm: fine-tunes a pretrained MLP autoencoder's
    encoder jointly with a DEC clustering head. Port of
    clustering_functions.dec_cluster / cluster_mlp_autoencoder."""

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

            idx = index_array[index * batch_size: min((index + 1) * batch_size, n)]
            if len(idx) < 1:
                index = index + 1 if (index + 1) * batch_size <= n else 0
                continue

            self.autoencoder.net.train()
            optimizer.zero_grad()
            _, z_batch = self.autoencoder.net.encode(x[idx])
            q_batch = self.head(z_batch)
            p_batch = torch.tensor(p[idx], dtype=torch.float32, device=self.device)
            # KL(p || q), matching Keras model.compile(loss='kld')
            loss = F.kl_div((q_batch + 1e-10).log(), p_batch, reduction='batchmean')
            loss.backward()
            optimizer.step()
            loss_val = loss.item()

            index = index + 1 if (index + 1) * batch_size <= n else 0

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
