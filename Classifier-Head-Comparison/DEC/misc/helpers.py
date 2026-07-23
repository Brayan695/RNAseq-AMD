# Small, framework-agnostic preprocessing/IO helpers shared by the plain-DEC scripts.
import os
import numpy as np


def normalizeRNA(X):
    """Min-max normalize a numeric matrix to [0, 1] (column-wise).

    Only called on the gene-expression block in run_dec.py, NOT on the
    clinical block - see the note in misc/dataset.py's _build_curated_clinical
    about 'age' consequently staying unscaled in the final feature matrix.
    """
    mins = X.min(axis=0)
    maxs = X.max(axis=0)
    # Genes that are constant would otherwise divide by zero.
    ranges = np.where(maxs - mins == 0, 1.0, maxs - mins)
    return (X - mins) / ranges


def save_results(savedir, savefile, **arrays):
    """Save any number of named arrays (e.g. z=latent_embedding, y_pred=cluster_labels,
    centroids=cluster_centers) to a single .npz file under savedir/savefile."""
    save_path = os.path.join(savedir, savefile)
    np.savez(save_path, **arrays)
