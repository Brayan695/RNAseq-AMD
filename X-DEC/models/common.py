# Same math building blocks as ../CNC-VAE/models/common.py (MMD and KL
# regularizers for pulling the VAE's latent space toward a standard normal
# prior, plus the reparameterization trick for differentiable latent
# sampling) - X-DEC's xvae uses the same regularizers, just on a dual-branch
# (numeric + binary) architecture instead of a single concatenated input.
import torch


def mmd(x, y):
    """Maximum Mean Discrepancy between samples x and y: near 0 when x and y
    look like they came from the same distribution, grows as they differ.
    Gaussian/RBF-kernel bandwidth 2*dim - the standard InfoVAE convention,
    matching CNC-VAE/models/common.py."""
    def kernel(a, b):
        dim = a.size(1)
        a = a.unsqueeze(1)
        b = b.unsqueeze(0)
        dist = (a - b).pow(2).sum(2)
        return torch.exp(-dist / (2.0 * dim))
    return kernel(x, x).mean() + kernel(y, y).mean() - 2 * kernel(x, y).mean()


def kl_regu(z_mean, z_log_sigma):
    """Per-sample KL( N(z_mean, exp(z_log_sigma)) || N(0, I) ) - the standard
    VAE regularization term. Caller averages over the batch (matches the
    original Keras kl_regu, which also left the batch-mean to the caller)."""
    kl = 1 + z_log_sigma - z_mean.pow(2) - z_log_sigma.exp()
    kl = torch.sum(kl, dim=-1)
    return -0.5 * kl


def sampling(z_mean, z_log_sigma):
    """Reparameterization trick: sample z ~ N(z_mean, exp(z_log_sigma)) in a
    way that's still differentiable w.r.t. z_mean/z_log_sigma."""
    std = torch.exp(0.5 * z_log_sigma)
    eps = torch.randn_like(std)
    return z_mean + eps * std
