import os
import numpy as np


def normalizeRNA(*args):
    """Min-max normalize RNA-seq data to [0, 1].

    Pass one array to normalize it on its own range, or two (train, test) to
    normalize both together on the range computed from their union.
    """
    if len(args) > 1:
        normalize_data = np.concatenate((args[0], args[1]), axis=0)
    else:
        normalize_data = args[0]

    mins = normalize_data.min(axis=0)
    maxs = normalize_data.max(axis=0)
    # Genes that are constant would otherwise divide by zero.
    ranges = np.where(maxs - mins == 0, 1.0, maxs - mins)

    if len(args) > 1:
        normalize_data = (normalize_data - mins) / ranges
        return normalize_data[:args[0].shape[0]], normalize_data[args[0].shape[0]:]
    else:
        return (normalize_data - mins) / ranges


def save_results(savedir, savefile, **arrays):
    save_path = os.path.join(savedir, savefile)
    np.savez(save_path, **arrays)
