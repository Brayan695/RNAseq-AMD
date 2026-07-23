"""
Port of the DAM-IC DEC repo's clustering_functions.compute_cluster_stability /
sample_cluster_stability / map_cluster_to_overlap (scripts/clustering_functions.py),
adapted for CNC-DEC's single concatenated-input CNCVAE.

Retrains a fresh CNCVAE + DEC under repeated K-fold resampling, maps each
resampled run's cluster labels back onto the reference (whole-cohort)
clustering by maximum sample overlap, and reports per-sample stability: the
fraction of resamples in which a sample landed back in "its" reference
cluster.

map_cluster_to_overlap here is a generalized, non-hardcoded-to-6-clusters
rewrite of the original (which indexed `[0,1,2,3,4,5]` regardless of
n_clusters) - functionally the same greedy max-overlap matching, sized to
whatever n_clusters is configured.
"""
import time

import numpy as np
import pandas as pd
from sklearn.model_selection import RepeatedKFold

from models.cncvae import CNCVAE
from models.dec import DEC


def map_cluster_to_overlap(y_pred, reference_labels, train_index, n_clusters):
    """Maps predicted cluster labels (from a resampled fit) onto reference
    cluster labels by maximum sample overlap (greedy one-to-one assignment)."""
    y_pred = np.asarray(y_pred)
    ref = np.asarray(reference_labels)[train_index]

    overlap = np.zeros((n_clusters, n_clusters))
    for new_cluster in range(n_clusters):
        mask = y_pred == new_cluster
        if mask.any():
            counts = np.bincount(ref[mask], minlength=n_clusters)
            overlap[:, new_cluster] = counts

    mapper = overlap.astype(float)
    mapper[mapper == 0] = np.nan

    y_pred_mapped = np.full(len(y_pred), np.nan)
    while not np.all(np.isnan(mapper)):
        flat_idx = np.nanargmax(mapper)
        original_cluster, new_cluster = np.unravel_index(flat_idx, mapper.shape)
        y_pred_mapped[y_pred == new_cluster] = original_cluster
        mapper[original_cluster, :] = np.nan
        mapper[:, new_cluster] = np.nan

    # Clusters with zero overlap with any reference cluster (all-NaN column)
    # keep their original label rather than being dropped.
    unmapped = np.isnan(y_pred_mapped)
    y_pred_mapped[unmapped] = y_pred[unmapped]
    return y_pred_mapped


def sample_cluster_stability(y_pred_mat, reference_labels):
    """Per-sample fraction of resampling iterations where the (overlap-mapped)
    predicted cluster matched the reference (whole-cohort) cluster label.
    NaN entries (sample held out of that iteration's training fold) are
    excluded from the per-sample average rather than counted as mismatches."""
    reference = pd.Series(reference_labels, index=y_pred_mat.index)
    matches = y_pred_mat.eq(reference, axis=0)
    matches = matches.astype(float)
    matches[y_pred_mat.isna()] = np.nan
    return matches.mean(axis=1)


def compute_cluster_stability(s1, s2, n_clusters, reference_labels, cncvae_kwargs=None,
                              dec_kwargs=None, k=5, rep=2, seed=5192, verbose=True):
    """Retrains CNCVAE + DEC under repeated K-fold resampling and reports
    per-sample cluster stability against the reference (whole-cohort) fit.

    Parameters
    ----------
    s1, s2 : ndarray
        The two input blocks (concatenated internally by CNCVAE - order
        doesn't matter for this model, unlike X-DEC).
    n_clusters : int
    reference_labels : ndarray
        Cluster labels from the whole-cohort DEC fit, used both as the
        stability reference and to map each resample's clusters onto it.
    cncvae_kwargs : dict, optional
        Keyword args forwarded to `CNCVAE(...)` (ds, ls, epochs, ...).
    dec_kwargs : dict, optional
        Keyword args forwarded to `DEC.fit(...)` (maxiter, batch_size, ...).
    k, rep : int
        Repeated K-fold parameters. Original repo defaults to k=10, rep=5
        (tuned for cohorts of thousands); default here is scaled down for
        n~150-166 and CPU-only training.

    Returns
    -------
    stability : pandas.Series
        Per-sample stability score in [0, 1].
    y_pred_mat : pandas.DataFrame
        Raw per-iteration (overlap-mapped) cluster predictions, samples x
        iterations, NaN where a sample was held out of that iteration.
    """
    cncvae_kwargs = cncvae_kwargs or {}
    dec_kwargs = dict(dec_kwargs or {})
    dec_kwargs.pop('y', None)
    dec_kwargs['verbose'] = False

    n = s1.shape[0]
    y_pred_mat = pd.DataFrame(index=np.arange(n), columns=range(k * rep), dtype=float)

    rkf = RepeatedKFold(n_splits=k, n_repeats=rep, random_state=seed)
    for i, (train_index, _) in enumerate(rkf.split(s1)):
        t0 = time.time()
        s1_train, s2_train = s1[train_index], s2[train_index]

        model = CNCVAE(input_size=s1_train.shape[1] + s2_train.shape[1], **cncvae_kwargs)
        model.build_model(seed=seed + i)
        model.train(s1_train, s2_train, seed=seed + i, verbose=False)

        dec = DEC(model, n_clusters=n_clusters, seed=seed + i)
        y_pred_it, _, _, _ = dec.fit(s1_train, s2_train, y=None, **dec_kwargs)

        y_pred_mat.iloc[train_index, i] = map_cluster_to_overlap(
            y_pred_it, reference_labels, train_index, n_clusters)

        if verbose:
            print('stability iteration {}/{} took {:.1f}s'.format(i + 1, k * rep, time.time() - t0))

    stability = sample_cluster_stability(y_pred_mat, reference_labels)
    return stability, y_pred_mat
