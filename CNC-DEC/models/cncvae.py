"""
Port of ../CNC-VAE/models/cncvae.py (a single concatenated-input variational
encoder, KL/MMD-regularized, MSE reconstruction) - the "CNC" half of
CNC-DEC. Architecture is unchanged from CNC-VAE; only the class API is
adapted to build its own `self.args` from constructor kwargs (matching the
xvae/MLPAutoencoder style used in ../X-DEC and ../DEC) instead of requiring
an externally-built argparse.Namespace, so it's self-contained.

Unlike the plain DEC autoencoder (`AENet.encode()` returns `(h1, z)` - use
the second element) this net's `encode()` returns `(z_mean, z_log_sigma)` -
use the FIRST element for the deterministic latent code. models/dec.py in
this folder unpacks it accordingly.
"""
import os

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from models.common import kl_regu, mmd, sampling

_ACTIVATIONS = {'elu': nn.ELU, 'relu': nn.ReLU, 'tanh': nn.Tanh, 'sigmoid': nn.Sigmoid}


class _CNCVAENet(nn.Module):
    """Encoder: Linear -> BatchNorm -> activation -> {z_mean, z_log_sigma} -> reparameterize.
    Decoder: Linear -> BatchNorm -> activation -> Dropout -> Linear (no output
    activation - MSE reconstruction over the whole concatenated input, matching
    CNC-VAE's Clin+mRNA integration branch)."""

    def __init__(self, input_size, ds, ls, dropout, act='elu'):
        super().__init__()
        act_cls = _ACTIVATIONS.get(act, nn.ELU)

        self.enc_body = nn.Sequential(
            nn.Linear(input_size, ds),
            nn.BatchNorm1d(ds),
            act_cls(),
        )
        self.fc_mu = nn.Linear(ds, ls)
        self.fc_logvar = nn.Linear(ds, ls)
        nn.init.zeros_(self.fc_logvar.weight)
        nn.init.zeros_(self.fc_logvar.bias)

        self.decoder = nn.Sequential(
            nn.Linear(ls, ds),
            nn.BatchNorm1d(ds),
            act_cls(),
            nn.Dropout(dropout),
            nn.Linear(ds, input_size),
        )

    def encode(self, x):
        h = self.enc_body(x)
        return self.fc_mu(h), self.fc_logvar(h)

    def forward(self, x):
        z_mean, z_log_sigma = self.encode(x)
        z = sampling(z_mean, z_log_sigma)
        x_hat = self.decoder(z)
        return x_hat, z_mean, z_log_sigma, z


class CNCVAE:
    """Single concatenated-input VAE. Inputs s1/s2 are just concatenated
    (order doesn't matter, no type-specific branches/losses) - matching the
    original CNC-VAE, which is why it's the "simple" counterpart to X-DEC's
    dual-branch xvae rather than a mixed-data-type-aware model."""

    def __init__(self, input_size, ds=64, ls=16, act='elu', dropout=0.2,
                 distance='kl', beta=1.0, epochs=150, bs=32, save_model=False):
        self.args = _Args(
            input_size=input_size, ds=ds, ls=ls, act=act, dropout=dropout,
            distance=distance, beta=beta, epochs=epochs, bs=bs, save_model=save_model,
        )
        self.net = None
        self.optimizer = None
        self.device = torch.device('cpu')

    def build_model(self, seed=5192):
        torch.manual_seed(seed)
        np.random.seed(seed)
        a = self.args
        self.net = _CNCVAENet(a.input_size, a.ds, a.ls, a.dropout, a.act).to(self.device)
        self.optimizer = torch.optim.Adam(self.net.parameters(), lr=0.001, betas=(0.9, 0.999))

    def _loss(self, x, x_hat, z_mean, z_log_sigma, z):
        a = self.args
        recon = F.mse_loss(x_hat, x, reduction='mean')
        if a.distance == 'mmd':
            prior = torch.randn_like(z)
            distance = mmd(prior, z)
        elif a.distance == 'kl':
            distance = kl_regu(z_mean, z_log_sigma).mean()
        else:
            raise ValueError('{} not recognised as distance.'.format(a.distance))
        return recon + a.beta * distance

    def train(self, s1_train, s2_train, s1_test=None, s2_test=None, seed=5192, verbose=True):
        torch.manual_seed(seed)
        np.random.seed(seed)
        a = self.args

        train_x = torch.tensor(
            np.concatenate((np.asarray(s1_train), np.asarray(s2_train)), axis=-1),
            dtype=torch.float32, device=self.device)
        test_x = None
        if s1_test is not None and s2_test is not None:
            test_x = torch.tensor(
                np.concatenate((np.asarray(s1_test), np.asarray(s2_test)), axis=-1),
                dtype=torch.float32, device=self.device)

        n = train_x.shape[0]
        for epoch in range(a.epochs):
            self.net.train()
            perm = torch.randperm(n)
            epoch_loss, n_batches = 0.0, 0
            for i in range(0, n, a.bs):
                idx = perm[i:i + a.bs]
                if idx.numel() < 2:  # BatchNorm needs more than one sample
                    continue
                batch = train_x[idx]
                self.optimizer.zero_grad()
                x_hat, z_mean, z_log_sigma, z = self.net(batch)
                loss = self._loss(batch, x_hat, z_mean, z_log_sigma, z)
                loss.backward()
                self.optimizer.step()
                epoch_loss += loss.item()
                n_batches += 1

            if verbose and (epoch + 1) % max(1, a.epochs // 20) == 0:
                if test_x is not None:
                    self.net.eval()
                    with torch.no_grad():
                        x_hat, z_mean, z_log_sigma, z = self.net(test_x)
                        val_loss = self._loss(test_x, x_hat, z_mean, z_log_sigma, z).item()
                    print('Epoch {}/{} - loss: {:.4f} - val_loss: {:.4f}'.format(
                        epoch + 1, a.epochs, epoch_loss / max(n_batches, 1), val_loss))
                else:
                    print('Epoch {}/{} - loss: {:.4f}'.format(
                        epoch + 1, a.epochs, epoch_loss / max(n_batches, 1)))

        if a.save_model:
            self.save_encoder('cncvae_model/model_weights.pt')

    def predict(self, s1_data, s2_data):
        """Returns the deterministic latent code z_mean (n_samples, ls)."""
        x = torch.tensor(
            np.concatenate((np.asarray(s1_data), np.asarray(s2_data)), axis=1),
            dtype=torch.float32, device=self.device)
        self.net.eval()
        with torch.no_grad():
            z_mean, _ = self.net.encode(x)
        return z_mean.cpu().numpy()

    def save_encoder(self, path, force_path=True):
        """Saves ONE self-contained checkpoint file: state_dict plus the
        architecture hyperparameters needed to rebuild `_CNCVAENet`, all as
        top-level keys in a single dict - matches ../CNC-VAE/models/cncvae.py's
        save_encoder() format exactly, so existing SHAP/analysis notebooks
        written against CNC-VAE's checkpoints
        (`ckpt = torch.load(path); net = _CNCVAENet(ckpt['input_size'], ckpt['ds'],
        ckpt['ls'], ckpt['dropout'], ckpt['act']); net.load_state_dict(ckpt['state_dict'])`)
        work unchanged against CNC-DEC's checkpoints too.
        """
        directory = os.path.dirname(path) or '.'
        if force_path and not os.path.exists(directory):
            os.makedirs(directory)

        a = self.args
        torch.save({
            'state_dict': self.net.state_dict(),
            'input_size': a.input_size,
            'ds': a.ds,
            'ls': a.ls,
            'dropout': a.dropout,
            'act': a.act,
        }, path)


class _Args:
    """Trivial namespace object - just turns constructor kwargs into
    attributes, e.g. _Args(ds=64).ds == 64. Used so CNCVAE can store its
    hyperparameters as self.args.<name> (same access pattern as CNC-VAE's
    argparse.Namespace) without requiring an actual argparse.Namespace to be
    built externally."""
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)


def load_cncvae_model(path):
    """Load a saved CNC-VAE encoder/decoder (weights saved via CNCVAE.save_encoder).
    Rebuilds a full CNCVAE wrapper (not just the raw _CNCVAENet) so the
    returned object's .predict()/.train() etc. all work immediately."""
    ckpt = torch.load(path, map_location='cpu')
    model = CNCVAE(input_size=ckpt['input_size'], ds=ckpt['ds'], ls=ckpt['ls'],
                   act=ckpt['act'], dropout=ckpt['dropout'])
    model.build_model()
    model.net.load_state_dict(ckpt['state_dict'])
    model.net.eval()
    return model
