# Small math building blocks shared by the VAE model: the two alternative
# "keep the latent space close to a standard normal prior" regularizers (MMD or
# KL divergence - pick one via --distance), plus the reparameterization trick
# used to sample a latent vector z in a way that's still differentiable.
import torch


def mmd(x, y):
    """Gaussian/RBF-kernel MMD (Maximum Mean Discrepancy) between samples x and y,
    bandwidth 2*dim - the standard InfoVAE convention, matching the CNC_VAE_AMD
    reference notebooks.

    Intuition: mmd(x, y) is near 0 when x and y look like they came from the same
    distribution, and grows as they look more different. Used as an alternative
    to KL divergence for pulling the encoder's latent output z toward samples
    drawn from a standard normal distribution ("prior").
    """
    def kernel(a, b):
        # Pairwise Gaussian-kernel similarity between every row of a and every
        # row of b: kernel(a,b)[i,j] = exp(-||a_i - b_j||^2 / bandwidth).
        dim = a.size(1)
        a = a.unsqueeze(1)  # (n_a, 1, dim)
        b = b.unsqueeze(0)  # (1, n_b, dim)
        dist = (a - b).pow(2).sum(2)  # (n_a, n_b) squared distance between every pair
        return torch.exp(-dist / (2.0 * dim))
    # MMD^2 estimator: mean within-x similarity + mean within-y similarity
    # - 2 * mean cross similarity. Zero if x and y are identically distributed.
    return kernel(x, x).mean() + kernel(y, y).mean() - 2 * kernel(x, y).mean()


def kl_regu(z_mean, z_log_sigma):
    """Per-sample KL( N(z_mean, exp(z_log_sigma)) || N(0, I) ). Caller averages
    over the batch (matches the original Keras kl_regu, which also left the
    batch-mean to the caller).

    This is the standard VAE regularization term: it penalizes the encoder for
    producing a latent distribution (mean=z_mean, variance=exp(z_log_sigma)) that
    strays far from a standard normal N(0, I). Closed-form KL divergence between
    two Gaussians - no sampling/approximation needed.
    """
    kl = 1 + z_log_sigma - z_mean.pow(2) - z_log_sigma.exp()
    kl = torch.sum(kl, dim=-1)  # sum over latent dimensions -> one value per sample
    return -0.5 * kl


def sampling(z_mean, z_log_sigma):
    """Reparameterization trick: sample from an isotropic unit Gaussian.

    Draws z ~ N(z_mean, exp(z_log_sigma)) in a way that gradients can still flow
    back through z_mean/z_log_sigma during training - done by sampling noise
    (eps) independently and combining it with z_mean/std via simple arithmetic,
    rather than sampling z directly (which wouldn't be differentiable).
    """
    std = torch.exp(0.5 * z_log_sigma)
    eps = torch.randn_like(std)  # random noise, same shape as std, not connected to the graph
    return z_mean + eps * std
