# Small, framework-agnostic preprocessing/IO helpers shared by the training scripts.
import os
import numpy as np


def normalizeRNA(*args):
    """Min-max normalize RNA-seq data to [0, 1].

    Pass one array to normalize it on its own range, or two (train, test) to
    normalize both together on the range computed from their union - this
    matches the original CNC-VAE preprocessing.
    """
    if len(args) > 1:
        # Stack train+test so both are scaled by the SAME min/max per gene -
        # scaling them separately would put train and test on different scales.
        normalize_data = np.concatenate((args[0], args[1]), axis=0)
    else:
        normalize_data = args[0]

    mins = normalize_data.min(axis=0)
    maxs = normalize_data.max(axis=0)
    # Genes that are constant within a fold would otherwise divide by zero.
    ranges = np.where(maxs - mins == 0, 1.0, maxs - mins)

    if len(args) > 1:
        normalize_data = (normalize_data - mins) / ranges
        # Split back apart into (train, test) using the original row counts.
        return normalize_data[:args[0].shape[0]], normalize_data[args[0].shape[0]:]
    else:
        return (normalize_data - mins) / ranges


def save_embedding(savedir, savefile, *args):
    """Save one (whole-dataset) or two (train, test) latent embeddings to a
    single .npz file under savedir/savefile."""
    save_path = os.path.join(savedir, savefile)
    if len(args) > 1:
        np.savez(save_path, emb_train=args[0], emb_test=args[1])
    else:
        np.savez(save_path, emb_train=args[0])
