import os
import numpy as np


def normalizeRNA(X):
    """Min-max normalize a numeric matrix to [0, 1] (column-wise)."""
    mins = X.min(axis=0)
    maxs = X.max(axis=0)
    # Genes that are constant would otherwise divide by zero.
    ranges = np.where(maxs - mins == 0, 1.0, maxs - mins)
    return (X - mins) / ranges


def save_results(savedir, savefile, **arrays):
    save_path = os.path.join(savedir, savefile)
    np.savez(save_path, **arrays)
