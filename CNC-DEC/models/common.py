import torch


def mmd(x, y):
    """Gaussian/RBF-kernel MMD between samples x and y, bandwidth 2*dim - the
    standard InfoVAE convention, matching CNC-VAE/models/common.py."""
    def kernel(a, b):
        dim = a.size(1)
        a = a.unsqueeze(1)
        b = b.unsqueeze(0)
        dist = (a - b).pow(2).sum(2)
        return torch.exp(-dist / (2.0 * dim))
    return kernel(x, x).mean() + kernel(y, y).mean() - 2 * kernel(x, y).mean()


def kl_regu(z_mean, z_log_sigma):
    """Per-sample KL( N(z_mean, exp(z_log_sigma)) || N(0, I) ). Caller averages
    over the batch (matches the original Keras kl_regu, which also left the
    batch-mean to the caller)."""
    kl = 1 + z_log_sigma - z_mean.pow(2) - z_log_sigma.exp()
    kl = torch.sum(kl, dim=-1)
    return -0.5 * kl


def sampling(z_mean, z_log_sigma):
    """Reparameterization trick: sample from an isotropic unit Gaussian."""
    std = torch.exp(0.5 * z_log_sigma)
    eps = torch.randn_like(std)
    return z_mean + eps * std
