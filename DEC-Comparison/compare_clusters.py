"""
Compares X-DEC's and CNC-DEC's clustering results on the same AMD cohort:

  1. Clustering quality: ARI / NMI / accuracy of each model's clusters vs.
     the known mgs_level label (calculate_metrics, same Hungarian-matched
     accuracy as ../X-DEC and ../CNC-DEC use internally).
  2. Cross-model agreement: how much do the two models' cluster assignments
     agree with EACH OTHER, via a generalized (non-hardcoded-cluster-count)
     port of the DAM-IC repo's clustering_functions.jaccard_similarity -
     the same tool the original repo uses to compare clusterings, applied
     here to two models on one cohort instead of one model on two cohorts
     (see ../X-DEC/models/stability.py's docstring for why
     clustering_functions.compare_clusters() itself doesn't apply).
  3. Latent-space visualizations (PCA) for both models, colored by their own
     predicted cluster and by the true mgs_level.
  4. Cross-tab heatmaps: each model's clusters vs. mgs_level, and X-DEC's
     clusters vs. CNC-DEC's clusters.

Inputs are the `mgs_level.npz` files each model's run_*.py script already
saves (z, y_pred, y_proba, centroids, y, sample_id, classes).
"""
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

from cluster_compare_utils import calculate_metrics, jaccard_cluster_agreement

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--xdec_npz', required=True,
                     help='Path to X-DEC\'s saved results, e.g. '
                          '../X-DEC/results/XDEC_Clin+mRNA_integration/xdec_LS_16_.../mgs_level.npz')
parser.add_argument('--cncdec_npz', required=True,
                     help='Path to CNC-DEC\'s saved results, e.g. '
                          '../CNC-DEC/results/CNCDEC_Clin+mRNA_integration/cncdec_LS_16_.../mgs_level.npz')
parser.add_argument('--out', default='comparison', help='Output directory for plots/tables')


def load_npz(path):
    d = np.load(path, allow_pickle=True)
    return {k: d[k] for k in d.files}


def align_by_sample_id(a, b):
    """Reorders both result dicts onto their common sample_id set, in a
    consistent order, so per-sample comparisons line up."""
    common = np.intersect1d(a['sample_id'], b['sample_id'])
    if len(common) == 0:
        raise ValueError('X-DEC and CNC-DEC results have no overlapping sample_id - '
                         'were they trained on the same cohort/clinical_mode?')
    if len(common) < len(a['sample_id']) or len(common) < len(b['sample_id']):
        print('WARNING: only {} of {} (X-DEC) / {} (CNC-DEC) samples overlap - '
              'comparing on the intersection only.'.format(
                  len(common), len(a['sample_id']), len(b['sample_id'])))

    def reindex(d, ids):
        # get_indexer(ids) gives, for each id in `ids`, its row position in
        # d['sample_id'] - reordering every per-sample array by that index
        # list puts both dicts' arrays in the exact same sample order.
        order = pd.Index(d['sample_id']).get_indexer(ids)
        out = dict(d)
        for key in ('z', 'y_pred', 'y_proba', 'y', 'sample_id'):
            out[key] = d[key][order]
        return out

    return reindex(a, common), reindex(b, common)


def pca2d(z, seed=5192):
    return PCA(n_components=2, random_state=seed).fit_transform(z)


def scatter_panel(ax, proj, labels, label_names, title):
    for val in np.unique(labels):
        mask = labels == val
        name = label_names[val] if label_names is not None else str(val)
        ax.scatter(proj[mask, 0], proj[mask, 1], label=name, alpha=0.75, s=24)
    ax.set_title(title)
    ax.legend(fontsize=8)


def crosstab_heatmap(ax, row_labels, col_labels, row_names, col_names, title):
    ct = pd.crosstab(pd.Series(row_labels), pd.Series(col_labels))
    im = ax.imshow(ct.values, cmap='Blues')
    ax.set_xticks(range(len(ct.columns)))
    ax.set_xticklabels([col_names[c] if col_names is not None else c for c in ct.columns])
    ax.set_yticks(range(len(ct.index)))
    ax.set_yticklabels([row_names[r] if row_names is not None else r for r in ct.index])
    for i in range(ct.shape[0]):
        for j in range(ct.shape[1]):
            ax.text(j, i, ct.values[i, j], ha='center', va='center')
    ax.set_title(title)
    return ct


def main():
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    xdec = load_npz(args.xdec_npz)
    cncdec = load_npz(args.cncdec_npz)
    xdec, cncdec = align_by_sample_id(xdec, cncdec)
    n = len(xdec['sample_id'])
    print('Comparing on {} common samples.'.format(n))

    y = xdec['y']  # mgs_level encoded labels - identical for both after alignment
    classes = xdec['classes']

    # ----- 1. clustering quality vs. mgs_level -----
    xdec_acc, xdec_ari, xdec_nmi = calculate_metrics(y, xdec['y_pred'])
    cncdec_acc, cncdec_ari, cncdec_nmi = calculate_metrics(y, cncdec['y_pred'])
    print('X-DEC   vs mgs_level: acc={:.5f} ari={:.5f} nmi={:.5f}'.format(xdec_acc, xdec_ari, xdec_nmi))
    print('CNC-DEC vs mgs_level: acc={:.5f} ari={:.5f} nmi={:.5f}'.format(cncdec_acc, cncdec_ari, cncdec_nmi))

    # ----- 2. cross-model agreement -----
    n_clusters = max(xdec['y_pred'].max(), cncdec['y_pred'].max()) + 1
    jac_mat, mapping, cncdec_mapped = jaccard_cluster_agreement(
        xdec['y_pred'], cncdec['y_pred'], n_clusters)
    agreement = float(np.mean(xdec['y_pred'] == cncdec_mapped))
    cross_acc, cross_ari, cross_nmi = calculate_metrics(xdec['y_pred'], cncdec['y_pred'])
    print('X-DEC vs CNC-DEC agreement (best-match mapped): {:.5f}'.format(agreement))
    print('X-DEC vs CNC-DEC (label-invariant): ari={:.5f} nmi={:.5f}'.format(cross_ari, cross_nmi))

    summary = pd.DataFrame([
        {'model': 'X-DEC', 'acc_vs_mgs_level': xdec_acc, 'ari_vs_mgs_level': xdec_ari,
         'nmi_vs_mgs_level': xdec_nmi},
        {'model': 'CNC-DEC', 'acc_vs_mgs_level': cncdec_acc, 'ari_vs_mgs_level': cncdec_ari,
         'nmi_vs_mgs_level': cncdec_nmi},
        {'model': 'X-DEC vs CNC-DEC', 'acc_vs_mgs_level': agreement,
         'ari_vs_mgs_level': cross_ari, 'nmi_vs_mgs_level': cross_nmi},
    ]).set_index('model')
    summary.to_csv(os.path.join(args.out, 'comparison_summary.csv'))
    print(summary)

    # ----- 3. latent-space PCA panels -----
    xdec_proj = pca2d(xdec['z'])
    cncdec_proj = pca2d(cncdec['z'])

    fig, axes = plt.subplots(2, 2, figsize=(11, 10))
    scatter_panel(axes[0, 0], xdec_proj, xdec['y_pred'],
                 {i: 'cluster {}'.format(i) for i in range(n_clusters)}, 'X-DEC latent - predicted cluster')
    scatter_panel(axes[0, 1], xdec_proj, y, dict(enumerate(classes)), 'X-DEC latent - true mgs_level')
    scatter_panel(axes[1, 0], cncdec_proj, cncdec['y_pred'],
                 {i: 'cluster {}'.format(i) for i in range(n_clusters)}, 'CNC-DEC latent - predicted cluster')
    scatter_panel(axes[1, 1], cncdec_proj, y, dict(enumerate(classes)), 'CNC-DEC latent - true mgs_level')
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, 'latent_space_comparison.png'), bbox_inches='tight')
    print('Saved', os.path.join(args.out, 'latent_space_comparison.png'))

    # ----- 4. cross-tab heatmaps -----
    class_names = {i: c for i, c in enumerate(classes)}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    crosstab_heatmap(axes[0], xdec['y_pred'], y, None, class_names, 'X-DEC cluster vs. mgs_level')
    crosstab_heatmap(axes[1], cncdec['y_pred'], y, None, class_names, 'CNC-DEC cluster vs. mgs_level')
    crosstab_heatmap(axes[2], xdec['y_pred'], cncdec['y_pred'], None, None, 'X-DEC cluster vs. CNC-DEC cluster')
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, 'crosstab_heatmaps.png'), bbox_inches='tight')
    print('Saved', os.path.join(args.out, 'crosstab_heatmaps.png'))


if __name__ == '__main__':
    main()
