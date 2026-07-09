"""
Shared utilities for compare_clusters.py and stability_compare.py.

- calculate_metrics: same Hungarian-matched accuracy / ARI / NMI port used
  internally by ../X-DEC, ../CNC-DEC, and ../DEC (scripts/metrics.py in the
  DAM-IC repo).
- jaccard_cluster_agreement: generalized (non-hardcoded-to-6-clusters) port
  of the DAM-IC repo's clustering_functions.jaccard_similarity - computes
  per-cluster-pair Jaccard scores between two clusterings of the SAME
  samples and greedily maps one clustering's labels onto the other's by
  maximum Jaccard overlap.
"""
import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


def cluster_accuracy(y_true, y_pred):
    y_true = np.asarray(y_true).astype(np.int64)
    y_pred = np.asarray(y_pred).astype(np.int64)
    assert y_pred.size == y_true.size
    d = max(y_pred.max(), y_true.max()) + 1
    w = np.zeros((d, d), dtype=np.int64)
    for i in range(y_pred.size):
        w[y_pred[i], y_true[i]] += 1
    row_ind, col_ind = linear_sum_assignment(w.max() - w)
    return sum(w[i, j] for i, j in zip(row_ind, col_ind)) / y_pred.size


def calculate_metrics(y_true, y_pred):
    acc = np.round(cluster_accuracy(y_true, y_pred), 5)
    nmi = np.round(normalized_mutual_info_score(y_true, y_pred), 5)
    ari = np.round(adjusted_rand_score(y_true, y_pred), 5)
    return acc, ari, nmi


def jaccard_cluster_agreement(labels_a, labels_b, n_clusters):
    """Computes the Jaccard similarity between every pair of clusters from
    two different clusterings of the same samples, then greedily maps
    clusters in `labels_b` onto clusters in `labels_a` by maximum Jaccard
    score.

    Returns
    -------
    jac : ndarray, shape (n_clusters, n_clusters)
        jac[i, j] = Jaccard(cluster i of a, cluster j of b).
    mapping : dict
        {b_cluster_label: a_cluster_label}.
    labels_b_mapped : ndarray
        `labels_b` relabelled onto `labels_a`'s cluster numbering.
    """
    labels_a = np.asarray(labels_a)
    labels_b = np.asarray(labels_b)

    jac = np.zeros((n_clusters, n_clusters))
    for i in range(n_clusters):
        a_mask = labels_a == i
        for j in range(n_clusters):
            b_mask = labels_b == j
            union = np.sum(a_mask | b_mask)
            jac[i, j] = np.sum(a_mask & b_mask) / union if union else 0.0

    mapper = jac.copy()
    mapping = {}
    for _ in range(n_clusters):
        if np.all(mapper <= 0):
            break
        i, j = np.unravel_index(np.argmax(mapper), mapper.shape)
        mapping[int(j)] = int(i)
        mapper[i, :] = -1
        mapper[:, j] = -1

    labels_b_mapped = np.array([mapping.get(int(lbl), int(lbl)) for lbl in labels_b])
    return jac, mapping, labels_b_mapped
