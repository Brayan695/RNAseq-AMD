"""PyTorch/numpy port of the DAM-IC DEC repo's scripts/metrics.py clustering
metrics (accuracy via Hungarian matching, NMI, ARI)."""
import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score


def cluster_accuracy(y_true, y_pred):
    """Best-match clustering accuracy: finds the cluster-to-label assignment
    that maximises overlap (Hungarian algorithm), same as the original
    scripts/metrics.py accuracy()."""
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
