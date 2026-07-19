"""
PyTorch port of the DAM-IC DEC repo's scripts/xvae.py (X-shaped VAE for mixed
numerical + binary data types). Ported to PyTorch for the same reason
CNC-VAE was: TensorFlow/Keras has DLL-loading problems on this Windows machine.

Set 1 (s1) must always be numerical (here: RNA-seq CPM expression) and is
reconstructed with MSE / a linear output head. Set 2 (s2) must always be
binary/one-hot categorical (here: MetaSheet_1_4.csv clinical/genotype dummy
columns) and is reconstructed with binary cross-entropy / a sigmoid output
head. The two branches are encoded separately, concatenated, and passed
through a shared bottleneck before the latent z_mean / z_log_sigma heads -
same topology as the original Keras xvae.build_model().
"""
import os

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from models.common import kl_regu, mmd, sampling

_ACTIVATIONS = {'elu': nn.ELU, 'relu': nn.ReLU, 'tanh': nn.Tanh, 'sigmoid': nn.Sigmoid}


class _XVAENet(nn.Module):
    """The "X" shape: two separate encoder branches (one per data type) that
    merge into a shared bottleneck before the latent heads, then a shared
    decoder trunk that splits back into two type-specific output branches -
    contrast with CNC-VAE, which just concatenates everything into one branch
    upfront. This lets each data type have its own appropriately-scaled
    hidden representation before being combined."""

    def __init__(self, s1_input_size, s2_input_size, ds1, ds2, ds12, ls, dropout, act='elu'):
        super().__init__()
        act_cls = _ACTIVATIONS.get(act, nn.ELU)

        # Two independent branches, one per input type, each own Dense+BatchNorm+activation...
        self.enc1 = nn.Sequential(nn.Linear(s1_input_size, ds1), nn.BatchNorm1d(ds1), act_cls())
        self.enc2 = nn.Sequential(nn.Linear(s2_input_size, ds2), nn.BatchNorm1d(ds2), act_cls())
        # ...concatenated and passed through a shared bottleneck layer before
        # the two latent heads - this is where the "X" shape's branches merge.
        self.enc12 = nn.Sequential(nn.Linear(ds1 + ds2, ds12), nn.BatchNorm1d(ds12), act_cls())

        self.fc_mu = nn.Linear(ds12, ls)
        self.fc_logvar = nn.Linear(ds12, ls)
        # matches Keras z_log_sigma Dense(kernel_initializer='zeros')
        nn.init.zeros_(self.fc_logvar.weight)
        nn.init.zeros_(self.fc_logvar.bias)

        # Mirror of the encoder: shared trunk from the latent code, then split
        # back into two type-specific output branches.
        self.dec_shared = nn.Sequential(
            nn.Linear(ls, ds12), nn.BatchNorm1d(ds12), act_cls(), nn.Dropout(dropout)
        )
        self.dec1 = nn.Sequential(nn.Linear(ds12, ds1), nn.BatchNorm1d(ds1), act_cls())
        self.dec2 = nn.Sequential(nn.Linear(ds12, ds2), nn.BatchNorm1d(ds2), act_cls())
        self.out1 = nn.Linear(ds1, s1_input_size)  # linear activation (numeric)
        self.out2 = nn.Linear(ds2, s2_input_size)  # logits; sigmoid applied via BCEWithLogits

    def encode(self, x1, x2):
        """x1 = numeric input (RNA-seq), x2 = binary input (clinical dummies).
        Returns (z_mean, z_log_sigma) - the deterministic distribution
        parameters, not a sample."""
        h1 = self.enc1(x1)
        h2 = self.enc2(x2)
        h = self.enc12(torch.cat([h1, h2], dim=-1))
        return self.fc_mu(h), self.fc_logvar(h)

    def decode(self, z):
        """Returns (s1_out, s2_logits): s1_out is the numeric reconstruction
        (already on the right scale, no activation needed), s2_logits are RAW
        LOGITS for the binary reconstruction - call torch.sigmoid() on them to
        get actual [0, 1] probabilities (kept as logits here because the loss
        function uses the numerically-stabler BCEWithLogits variant)."""
        h = self.dec_shared(z)
        s1_out = self.out1(self.dec1(h))
        s2_logits = self.out2(self.dec2(h))
        return s1_out, s2_logits

    def forward(self, x1, x2):
        z_mean, z_log_sigma = self.encode(x1, x2)
        z = sampling(z_mean, z_log_sigma)
        s1_out, s2_logits = self.decode(z)
        return s1_out, s2_logits, z_mean, z_log_sigma, z


class xvae:
    """X-shaped VAE. SET 1 must always be numerical, SET 2 must always be
    binary (categorical) - same convention as the original repo."""

    def __init__(self, s1_input_size, s2_input_size, ds1=48, ds2=None, ds12=None,
                 ls=32, weighted=True, act='elu', dropout=0.2, distance='kl',
                 beta=25, epochs=250, bs=64, save_model=False):
        self.args = _Args(
            s1_input_size=s1_input_size, s2_input_size=s2_input_size,
            ds1=ds1, ds2=ds1 if ds2 is None else ds2,
            ds12=ds1 if ds12 is None else ds12,
            ls=ls, act=act, dropout=dropout, distance=distance, beta=beta,
            epochs=epochs, bs=bs, weighted=weighted, save_model=save_model,
        )
        self.net = None
        self.optimizer = None
        self.device = torch.device('cpu')

    def build_model(self, seed=5192):
        torch.manual_seed(seed)
        np.random.seed(seed)
        a = self.args
        self.net = _XVAENet(a.s1_input_size, a.s2_input_size, a.ds1, a.ds2,
                            a.ds12, a.ls, a.dropout, a.act).to(self.device)
        self.optimizer = torch.optim.Adam(self.net.parameters(), lr=0.001, betas=(0.9, 0.999))

    def _loss(self, x1, x2, s1_out, s2_logits, z_mean, z_log_sigma, z):
        """Reconstruction loss = MSE for the numeric branch (s1) + binary
        cross-entropy for the binary branch (s2), each type getting the loss
        function that actually matches its data - unlike CNC-VAE, which uses
        one MSE loss over everything concatenated together.

        `weighted` controls how the two branches' losses are combined: if
        True (default), they're pooled into one shared per-feature average
        (a branch with many more features naturally gets more total weight);
        if False, each branch is averaged independently first, so a
        small-but-important branch (e.g. few clinical features) doesn't get
        drowned out by a much larger one (e.g. many genes).
        """
        a = self.args
        # sum-over-features per sample, matching Keras mean(...)*n_features
        s1_loss = F.mse_loss(s1_out, x1, reduction='none').sum(dim=1)
        s2_loss = F.binary_cross_entropy_with_logits(s2_logits, x2, reduction='none').sum(dim=1)

        if a.weighted:
            reconstruction_loss = (s1_loss + s2_loss) / (a.s1_input_size + a.s2_input_size)
        else:
            reconstruction_loss = s1_loss / a.s1_input_size + s2_loss / a.s2_input_size

        if a.distance == 'mmd':
            prior = torch.randn_like(z)
            distance = mmd(prior, z)
            vae_loss = reconstruction_loss.mean() + a.beta * distance
        elif a.distance == 'kl':
            distance = kl_regu(z_mean, z_log_sigma)
            vae_loss = (reconstruction_loss + a.beta * distance).mean()
        else:
            raise ValueError('{} not recognised as distance.'.format(a.distance))
        return vae_loss

    def train(self, s1_train, s2_train, s1_test=None, s2_test=None, seed=5192, verbose=True):
        torch.manual_seed(seed)
        np.random.seed(seed)
        a = self.args

        x1 = torch.tensor(np.asarray(s1_train), dtype=torch.float32, device=self.device)
        x2 = torch.tensor(np.asarray(s2_train), dtype=torch.float32, device=self.device)
        test_x1 = test_x2 = None
        if s1_test is not None and s2_test is not None:
            test_x1 = torch.tensor(np.asarray(s1_test), dtype=torch.float32, device=self.device)
            test_x2 = torch.tensor(np.asarray(s2_test), dtype=torch.float32, device=self.device)

        n = x1.shape[0]
        for epoch in range(a.epochs):
            self.net.train()
            perm = torch.randperm(n)
            epoch_loss, n_batches = 0.0, 0
            for i in range(0, n, a.bs):
                idx = perm[i:i + a.bs]
                if idx.numel() < 2:  # BatchNorm needs more than one sample
                    continue
                b1, b2 = x1[idx], x2[idx]
                self.optimizer.zero_grad()
                s1_out, s2_logits, z_mean, z_log_sigma, z = self.net(b1, b2)
                loss = self._loss(b1, b2, s1_out, s2_logits, z_mean, z_log_sigma, z)
                loss.backward()
                self.optimizer.step()
                epoch_loss += loss.item()
                n_batches += 1

            if verbose and (epoch + 1) % max(1, a.epochs // 20) == 0:
                if test_x1 is not None:
                    self.net.eval()
                    with torch.no_grad():
                        s1_out, s2_logits, z_mean, z_log_sigma, z = self.net(test_x1, test_x2)
                        val_loss = self._loss(test_x1, test_x2, s1_out, s2_logits,
                                              z_mean, z_log_sigma, z).item()
                    print('Epoch {}/{} - loss: {:.4f} - val_loss: {:.4f}'.format(
                        epoch + 1, a.epochs, epoch_loss / max(n_batches, 1), val_loss))
                else:
                    print('Epoch {}/{} - loss: {:.4f}'.format(
                        epoch + 1, a.epochs, epoch_loss / max(n_batches, 1)))

        if a.save_model:
            self.save_encoder('xvae_model/model_weights.pt')

    def predict(self, s1_data, s2_data, output='encoder'):
        """output='encoder' returns the deterministic latent code z_mean
        (n_samples, ls). output='decoder' additionally reconstructs s1/s2
        (s2 passed through sigmoid)."""
        x1 = torch.tensor(np.asarray(s1_data), dtype=torch.float32, device=self.device)
        x2 = torch.tensor(np.asarray(s2_data), dtype=torch.float32, device=self.device)
        self.net.eval()
        with torch.no_grad():
            z_mean, _ = self.net.encode(x1, x2)
            if output == 'encoder':
                return z_mean.cpu().numpy()
            elif output == 'decoder':
                s1_out, s2_logits = self.net.decode(z_mean)
                return s1_out.cpu().numpy(), torch.sigmoid(s2_logits).cpu().numpy()
            raise ValueError("output must be 'encoder' or 'decoder'")

    def save_encoder(self, path, force_path=True):
        """Saves ONE self-contained checkpoint file: state_dict plus the
        architecture hyperparameters needed to rebuild `_XVAENet`, all as
        top-level keys in a single dict - same convention as
        ../CNC-VAE/models/cncvae.py's save_encoder(), so it loads the same
        simple way:
        `ckpt = torch.load(path); net = _XVAENet(ckpt['s1_input_size'], ckpt['s2_input_size'],
        ckpt['ds1'], ckpt['ds2'], ckpt['ds12'], ckpt['ls'], ckpt['dropout'], ckpt['act']);
        net.load_state_dict(ckpt['state_dict'])`. See XVAEEncoderMu below for a
        single-tensor-input wrapper (SHAP etc. expect one input, not two).
        """
        directory = os.path.dirname(path) or '.'
        if force_path and not os.path.exists(directory):
            os.makedirs(directory)

        a = self.args
        torch.save({
            'state_dict': self.net.state_dict(),
            's1_input_size': a.s1_input_size,
            's2_input_size': a.s2_input_size,
            'ds1': a.ds1,
            'ds2': a.ds2,
            'ds12': a.ds12,
            'ls': a.ls,
            'dropout': a.dropout,
            'act': a.act,
        }, path)


class _Args:
    def __init__(self, **kwargs):
        for k, v in kwargs.items():
            setattr(self, k, v)


def load_xvae_model(path):
    """Load a saved xvae encoder/decoder (weights saved via xvae.save_encoder)."""
    ckpt = torch.load(path, map_location='cpu')
    model = xvae(
        s1_input_size=ckpt['s1_input_size'], s2_input_size=ckpt['s2_input_size'],
        ds1=ckpt['ds1'], ds2=ckpt['ds2'], ds12=ckpt['ds12'], ls=ckpt['ls'],
        dropout=ckpt['dropout'], act=ckpt['act'],
    )
    model.build_model()
    model.net.load_state_dict(ckpt['state_dict'])
    model.net.eval()
    return model


class XVAEEncoderMu(nn.Module):
    """Single-tensor-input wrapper around `_XVAENet`'s encoder mean output,
    for tools (SHAP, etc.) that expect a model taking one input, not xvae's
    native two (s1 numeric, s2 binary). Splits the input into s1 (first
    `s1_input_size` columns) and s2 (remaining columns) before calling
    `net.encode(s1, s2)` - so `X = np.hstack([X_numeric, X_binary])` (same
    concatenation order used everywhere else in this project: run_xdec.py,
    models/dec.py) is the input this wrapper expects.

    Usage:
        ckpt = torch.load('encoder_xvae.pt', map_location='cpu')
        net = _XVAENet(ckpt['s1_input_size'], ckpt['s2_input_size'], ckpt['ds1'],
                       ckpt['ds2'], ckpt['ds12'], ckpt['ls'], ckpt['dropout'], ckpt['act'])
        net.load_state_dict(ckpt['state_dict']); net.eval()
        mu_model = XVAEEncoderMu(net, ckpt['s1_input_size'])
        explainer = shap.DeepExplainer(mu_model, background)  # background: torch.Tensor
    """

    def __init__(self, net, s1_input_size):
        super().__init__()
        self.net = net
        self.s1_input_size = s1_input_size

    def forward(self, x):
        x1 = x[:, :self.s1_input_size]
        x2 = x[:, self.s1_input_size:]
        z_mean, _ = self.net.encode(x1, x2)
        return z_mean
