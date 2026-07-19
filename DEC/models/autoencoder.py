"""
PyTorch port of the plain MLP autoencoder built inside the DAM-IC DEC repo's
clustering_functions.cluster_mlp_autoencoder() (scripts/clustering_functions.py).
This is the un-mixed-type "vanilla" DEC autoencoder - a single Dense input,
not the X-shaped dual-branch xvae used for the mixed numeric/binary variant.

Architecture (matches the original Keras layer-for-layer):
    input -> Dense(neurons_h, relu, L1 activity regularizer) -> Dense(neurons_e, relu)  [encoder]
    -> Dense(neurons_h, relu) -> Dense(n_features, sigmoid)                             [decoder]
Loss: MSE reconstruction + L1 penalty on the first hidden layer's activations,
optimized with Adam. Ported to PyTorch for the same reason CNC-VAE was:
TensorFlow/Keras has DLL-loading problems on this Windows machine.
"""
import os

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F


class _AENet(nn.Module):
    """A plain (non-variational) 4-layer autoencoder: input -> hidden -> latent
    -> hidden -> output, all ReLU except a sigmoid output (reconstruction
    target must be scaled to [0, 1] for that to make sense - see run_dec.py's
    normalizeRNA + already-binary clinical dummies)."""

    def __init__(self, n_features, neurons_h=64, neurons_e=8):
        super().__init__()
        self.hidden1 = nn.Linear(n_features, neurons_h)
        self.embed = nn.Linear(neurons_h, neurons_e)     # this is the latent/bottleneck layer
        self.dec_hidden = nn.Linear(neurons_e, neurons_h)
        self.out = nn.Linear(neurons_h, n_features)

    def encode(self, x):
        """Returns (h1, z): h1 is the first hidden layer's activations (used
        only for the L1 sparsity penalty in _loss below), z is the actual
        latent embedding DEC clusters on."""
        h1 = F.relu(self.hidden1(x))
        z = F.relu(self.embed(h1))
        return h1, z

    def decode(self, z):
        h = F.relu(self.dec_hidden(z))
        return torch.sigmoid(self.out(h))

    def forward(self, x):
        h1, z = self.encode(x)
        x_hat = self.decode(z)
        return x_hat, h1, z


class MLPAutoencoder:
    def __init__(self, n_features, neurons_h=64, neurons_e=8, epochs=500, bs=64,
                 l1_coef=10e-5, save_model=False):
        self.args = _Args(n_features=n_features, neurons_h=neurons_h, neurons_e=neurons_e,
                          epochs=epochs, bs=bs, l1_coef=l1_coef, save_model=save_model)
        self.net = None
        self.optimizer = None
        self.device = torch.device('cpu')

    def build_model(self, seed=5192):
        torch.manual_seed(seed)
        np.random.seed(seed)
        a = self.args
        self.net = _AENet(a.n_features, a.neurons_h, a.neurons_e).to(self.device)
        self.optimizer = torch.optim.Adam(self.net.parameters(), lr=0.001, betas=(0.9, 0.999))

    def _loss(self, x_hat, x, h1):
        """Reconstruction MSE + an L1 penalty on the first hidden layer's
        activations (encourages a sparse intermediate representation - this
        is the "L1 activity regularizer" mentioned in the module docstring)."""
        recon = F.mse_loss(x_hat, x, reduction='none').mean(dim=1)
        reg = self.args.l1_coef * h1.abs().sum(dim=1)
        return (recon + reg).mean()

    def train(self, X, seed=5192, verbose=True):
        """Pretrain the autoencoder with plain mini-batch Adam (no held-out
        validation set here, unlike CNC-VAE - this is unsupervised
        reconstruction pretraining before DEC.fit() takes over)."""
        torch.manual_seed(seed)
        np.random.seed(seed)
        a = self.args

        x = torch.tensor(np.asarray(X), dtype=torch.float32, device=self.device)
        n = x.shape[0]
        for epoch in range(a.epochs):
            self.net.train()
            perm = torch.randperm(n)
            epoch_loss, n_batches = 0.0, 0
            for i in range(0, n, a.bs):
                idx = perm[i:i + a.bs]
                b = x[idx]
                self.optimizer.zero_grad()
                x_hat, h1, z = self.net(b)
                loss = self._loss(x_hat, b, h1)
                loss.backward()
                self.optimizer.step()
                epoch_loss += loss.item()
                n_batches += 1
            if verbose and (epoch + 1) % max(1, a.epochs // 20) == 0:
                print('Epoch {}/{} - loss: {:.4f}'.format(epoch + 1, a.epochs, epoch_loss / max(n_batches, 1)))

        if a.save_model:
            self.save_encoder('dec_autoencoder/model_weights.pt')

    def predict(self, X):
        x = torch.tensor(np.asarray(X), dtype=torch.float32, device=self.device)
        self.net.eval()
        with torch.no_grad():
            _, z = self.net.encode(x)
        return z.cpu().numpy()

    def save_encoder(self, path, force_path=True):
        """Saves the weights to `path` and the architecture hyperparameters
        needed to rebuild _AENet to a SEPARATE sidecar CSV (`path +
        '.architecture.csv'`) - unlike CNC-VAE's save_encoder, which bundles
        everything into one .pt file. load_autoencoder() below expects both
        files to be present together."""
        directory = os.path.dirname(path) or '.'
        if force_path and not os.path.exists(directory):
            os.makedirs(directory)

        architecture = pd.DataFrame(index=['n_features', 'neurons_h', 'neurons_e', 'epochs', 'bs',
                                           'l1_coef'], columns=['value'])
        for i in architecture.index:
            architecture.loc[i, 'value'] = getattr(self.args, i)
        architecture.to_csv(path + '.architecture.csv')

        torch.save(self.net.state_dict(), path)


class _Args:
    """Trivial namespace object - turns constructor kwargs into attributes
    (e.g. _Args(epochs=500).epochs == 500) so MLPAutoencoder can store its
    hyperparameters as self.args.<name>."""
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)


def load_autoencoder(path):
    """Load a saved MLPAutoencoder (weights + architecture sidecar CSV saved
    via MLPAutoencoder.save_encoder). Rebuilds a full MLPAutoencoder wrapper,
    not just the raw _AENet."""
    architecture = pd.read_csv(path + '.architecture.csv', index_col=0)
    model = MLPAutoencoder(
        n_features=int(architecture.loc['n_features', 'value']),
        neurons_h=int(architecture.loc['neurons_h', 'value']),
        neurons_e=int(architecture.loc['neurons_e', 'value']),
        l1_coef=float(architecture.loc['l1_coef', 'value']),
    )
    model.build_model()
    model.net.load_state_dict(torch.load(path, map_location=model.device))
    model.net.eval()
    return model
