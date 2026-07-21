# The CNC-VAE model itself: a single-hidden-layer variational autoencoder that
# takes the CONCATENATED clinical + gene-expression vector as input ("CNC" =
# Concatenated iNputs Cnc-vae), compresses it down to a small latent embedding,
# and reconstructs the input from that embedding. The trained encoder's latent
# output is what downstream scripts use as a compact patient representation for
# classification, t-SNE, SHAP, etc.
#
# _CNCVAENet is the raw PyTorch network (what you'd load for inference elsewhere,
# e.g. in a SHAP script). CNCVAE is a thin wrapper around it that also owns the
# optimizer and the training loop, so run_cncvae.py/compare_representations.py
# can just do: cncvae = CNCVAE(args); cncvae.build_model(); cncvae.train(...).
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
        """
        input_size = number of concatenated clinical + gene features (the model's input/output width)
        ds         = hidden ("dense") layer width, shared by encoder and decoder
        ls         = latent size - how many numbers each patient gets compressed down to
        dropout    = dropout rate applied in the decoder only (regularization)
        act        = activation function name, looked up in _ACTIVATIONS above
        """
        super().__init__()
        act_cls = _ACTIVATIONS.get(act, nn.ELU)

        # Encoder body: input -> one hidden layer -> BatchNorm -> activation.
        # Output of this is fed into two separate heads (fc_mu, fc_logvar) below,
        # not directly used as the latent code itself.
        self.enc_body = nn.Sequential(
            nn.Linear(input_size, ds),
            nn.BatchNorm1d(ds),
            act_cls(),
        )
        # Two linear heads turn the hidden representation into the mean and
        # log-variance of a Gaussian distribution over the latent space - this
        # is what makes it a *variational* autoencoder rather than a plain one.
        self.fc_mu = nn.Linear(ds, ls)
        self.fc_logvar = nn.Linear(ds, ls)
        # Zero-initializing the log-variance head means every sample starts out
        # with unit variance (log(1)=0) - a stable starting point for training.
        nn.init.zeros_(self.fc_logvar.weight)
        nn.init.zeros_(self.fc_logvar.bias)

        # Decoder: latent vector -> hidden layer -> BatchNorm -> activation ->
        # Dropout -> back out to input_size. No output activation (e.g. no
        # sigmoid) since the reconstruction target (concatenated clinical+gene
        # features) isn't bounded to [0, 1], and the loss is plain MSE.
        self.decoder = nn.Sequential(
            nn.Linear(ls, ds),
            nn.BatchNorm1d(ds),
            act_cls(),
            nn.Dropout(dropout),
            nn.Linear(ds, input_size),
        )

    def encode(self, x):
        """Input -> (z_mean, z_log_sigma). This is the function to call for
        inference/downstream use (e.g. SHAP) - it's deterministic given x,
        unlike forward() which also draws a random sample."""
        h = self.enc_body(x)
        return self.fc_mu(h), self.fc_logvar(h)

    def forward(self, x):
        """Full VAE forward pass used during training: encode -> sample a latent
        vector z via the reparameterization trick -> decode z back to a
        reconstruction. Returns everything the loss function needs."""
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
        """Construct the network and optimizer. Must be called before train()/predict().
        Reads all its hyperparameters off self.args (ds, ls, dropout, act, input_size),
        which the calling script (run_cncvae.py etc.) is expected to have set."""
        torch.manual_seed(42)  # fixed seed -> reproducible weight initialization
        np.random.seed(42)
        self.net = _CNCVAENet(
            self.args.input_size, self.args.ds, self.args.ls, self.args.dropout, self.args.act
        ).to(self.device)
        self.optimizer = torch.optim.Adam(self.net.parameters(), lr=0.001, betas=(0.9, 0.999))

    def _loss(self, x, x_hat, z_mean, z_log_sigma, z):
        """VAE loss = reconstruction error + beta * (distance from the prior).
        distance is either KL divergence (--distance kl) or MMD (--distance mmd);
        beta controls how strongly the latent space is pulled toward the prior
        vs. how much weight goes on reconstructing the input accurately."""
        recon = F.binary_cross_entropy_with_logits(x_hat, x, reduction='mean')
        if self.args.distance == 'mmd':
            prior = torch.randn_like(z)  # samples from the target N(0, I) distribution
            distance = mmd(prior, z)
        else:
            distance = kl_regu(z_mean, z_log_sigma).mean()
        return recon + self.args.beta * distance

    def train(self, s1_train, s2_train, s1_test=None, s2_test=None):
        """Train the VAE with mini-batch gradient descent.

        s1_* and s2_* are the two input sources (clinical, gene expression) as
        separate numpy arrays - they're concatenated here into one input tensor,
        since CNC-VAE ("Concatenated iNputs") trains on the combined vector
        rather than treating the two sources specially. Passing s1_test/s2_test
        additionally reports a validation loss each epoch (no effect on training
        itself - the model never sees the test set during optimization).
        """
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
            self.net.train()  # enables dropout + batch-statistics-based BatchNorm
            perm = torch.randperm(n)  # shuffle sample order every epoch
            epoch_loss, n_batches = 0.0, 0
            for i in range(0, n, self.args.bs):
                idx = perm[i:i + self.args.bs]
                if idx.numel() < 2:  # BatchNorm needs more than one sample
                    continue
                batch = train_x[idx]
                self.optimizer.zero_grad()
                x_hat, z_mean, z_log_sigma, z = self.net(batch)
                loss = self._loss(batch, x_hat, z_mean, z_log_sigma, z)
                loss.backward()       # compute gradients
                self.optimizer.step()  # apply gradients (Adam update)
                epoch_loss += loss.item()
                n_batches += 1

            if test_x is not None:
                # Evaluate on the held-out set in inference mode (no dropout,
                # BatchNorm uses running statistics instead of this batch's own).
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
            # Weights-only save of the full VAE (encoder+decoder+optimizer state
            # is NOT included). For a standalone, reloadable encoder (e.g. for
            # SHAP), use save_encoder() instead.
            torch.save(self.net.state_dict(), self.args.model_out)

    def predict(self, s1_data, s2_data):
        """Encode (clinical, gene expression) into the latent embedding.

        Returns z_mean only (the deterministic center of the latent
        distribution) rather than a random sample - this is what downstream
        classifiers/plots/SHAP should use, since it's reproducible and doesn't
        inject sampling noise into results.
        """
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
