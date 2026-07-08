import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F

from models.common import kl_regu, mmd, sampling

_ACTIVATIONS = {'elu': nn.ELU, 'relu': nn.ReLU, 'tanh': nn.Tanh, 'sigmoid': nn.Sigmoid}


class _CNCVAENet(nn.Module):
    """Encoder: Linear -> BatchNorm -> activation -> {z_mean, z_log_sigma} -> reparameterize.
    Decoder: Linear -> BatchNorm -> activation -> Dropout -> Linear (no output activation,
    matching the original CancerAI-CL CNC-VAE Clin+mRNA integration branch - MSE
    reconstruction loss, not bounded to [0, 1] like the Clin+CNA/sigmoid branch).
    """

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
    """CNC-VAE for Clin+mRNA integration (AMD cohort has no CNA arm, unlike METABRIC).

    PyTorch port of the original tf.keras implementation - same architecture and
    training objective, switched to PyTorch to avoid Windows TensorFlow DLL
    loading issues and to match the reference notebooks' own framework.
    """

    def __init__(self, args):
        self.args = args
        self.net = None
        self.optimizer = None
        self.device = torch.device('cpu')

    def build_model(self):
        torch.manual_seed(42)
        np.random.seed(42)
        self.net = _CNCVAENet(
            self.args.input_size, self.args.ds, self.args.ls, self.args.dropout, self.args.act
        ).to(self.device)
        self.optimizer = torch.optim.Adam(self.net.parameters(), lr=0.001, betas=(0.9, 0.999))

    def _loss(self, x, x_hat, z_mean, z_log_sigma, z):
        recon = F.mse_loss(x_hat, x, reduction='mean')
        if self.args.distance == 'mmd':
            prior = torch.randn_like(z)
            distance = mmd(prior, z)
        else:
            distance = kl_regu(z_mean, z_log_sigma).mean()
        return recon + self.args.beta * distance

    def train(self, s1_train, s2_train, s1_test=None, s2_test=None):
        train_x = torch.tensor(
            np.concatenate((s1_train, s2_train), axis=-1), dtype=torch.float32, device=self.device
        )
        test_x = None
        if s1_test is not None and s2_test is not None:
            test_x = torch.tensor(
                np.concatenate((s1_test, s2_test), axis=-1), dtype=torch.float32, device=self.device
            )

        n = train_x.shape[0]
        for epoch in range(self.args.epochs):
            self.net.train()
            perm = torch.randperm(n)
            epoch_loss, n_batches = 0.0, 0
            for i in range(0, n, self.args.bs):
                idx = perm[i:i + self.args.bs]
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

            if test_x is not None:
                self.net.eval()
                with torch.no_grad():
                    x_hat, z_mean, z_log_sigma, z = self.net(test_x)
                    val_loss = self._loss(test_x, x_hat, z_mean, z_log_sigma, z).item()
                print('Epoch {}/{} - loss: {:.4f} - val_loss: {:.4f}'.format(
                    epoch + 1, self.args.epochs, epoch_loss / max(n_batches, 1), val_loss))
            else:
                print('Epoch {}/{} - loss: {:.4f}'.format(
                    epoch + 1, self.args.epochs, epoch_loss / max(n_batches, 1)))

        if self.args.save_model:
            torch.save(self.net.state_dict(), self.args.model_out)

    def predict(self, s1_data, s2_data):
        x = torch.tensor(
            np.concatenate((s1_data, s2_data), axis=1), dtype=torch.float32, device=self.device
        )
        self.net.eval()
        with torch.no_grad():
            z_mean, _ = self.net.encode(x)
        return z_mean.cpu().numpy()

    def save_encoder(self, path):
        """Save the encoder (architecture args + weights) for standalone reuse, e.g. SHAP.

        To reload:
            ckpt = torch.load(path, map_location='cpu')
            net = _CNCVAENet(ckpt['input_size'], ckpt['ds'], ckpt['ls'], ckpt['dropout'], ckpt['act'])
            net.load_state_dict(ckpt['state_dict'])
            net.eval()
        """
        torch.save({
            'state_dict': self.net.state_dict(),
            'input_size': self.args.input_size,
            'ds': self.args.ds,
            'ls': self.args.ls,
            'dropout': self.args.dropout,
            'act': self.args.act,
        }, path)
