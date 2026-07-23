# Binary classifier head trained on top of a frozen, already-trained encoder's
# latent embedding (z) - the encoder itself is never touched here; only this
# head's own weights get optimized. Architecture is picked by latent_dim,
# matching whichever --ls/neurons_e the loaded encoder was trained with:
#
#   latent_dim=16:  Linear(16, 4) -> ReLU -> Linear(4, 1) -> Sigmoid
#   latent_dim=8:   Linear(8, 4)  -> ReLU -> Linear(4, 1) -> Sigmoid
#   latent_dim=4:   Linear(4, 1)  -> Sigmoid                (no hidden layer)
#
# Outputs probabilities directly (Sigmoid is the last layer, not left as raw
# logits) - train with nn.BCELoss() to match.
import torch.nn as nn

_ARCHITECTURES = {
    16: lambda: nn.Sequential(nn.Linear(16, 4), nn.ReLU(), nn.Linear(4, 1), nn.Sigmoid()),
    8: lambda: nn.Sequential(nn.Linear(8, 4), nn.ReLU(), nn.Linear(4, 1), nn.Sigmoid()),
    4: lambda: nn.Sequential(nn.Linear(4, 1), nn.Sigmoid()),
}


class ClassifierHead(nn.Module):
    def __init__(self, latent_dim):
        super().__init__()
        if latent_dim not in _ARCHITECTURES:
            raise ValueError(
                'No classifier-head architecture defined for latent_dim={} - only 16, 8, '
                'and 4 are supported (matching the --ls/neurons_e values used throughout '
                'this project). Train an encoder with one of those latent sizes, or add a '
                'new entry to _ARCHITECTURES above if you need a different size.'.format(
                    latent_dim)
            )
        self.latent_dim = latent_dim
        self.net = _ARCHITECTURES[latent_dim]()

    def forward(self, z):
        """z: (n_samples, latent_dim) -> (n_samples,) predicted probabilities in [0, 1]."""
        return self.net(z).squeeze(-1)
