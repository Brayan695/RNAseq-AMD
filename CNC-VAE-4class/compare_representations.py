"""Reproduce the comparison table / t-SNE / embedding export from the
CNC_VAE_AMD_V5_4class reference notebook, using this project's tf.keras CNC-VAE
instead of the notebook's PyTorch model, so the two can be compared side by side.

Data: the 81 genes used in aak100_cpmdat.csv, pulled from gene_input.csv (matched
by expression VALUE on the 166 shared samples, since gene_input.csv only has gene
SYMBOLS while aak100_cpmdat.csv uses ENSG ids) for all 453 donors, joined with
clinical/genotype covariates from MetaSheet_Processed.csv. Target is the four-class
MGS1..MGS4 AMD grade (0..3).

Evaluation protocol (matches the notebook):
  1. Train one CNC-VAE on the whole cohort (no held-out fold) using the 'curated'
     clinical feature set (age/sex/genotype + prevalence-filtered, leakage-excluded
     oc_/mh_ history flags).
  2. Extract the latent mean z_mean for every sample (transductive - the same data
     used to fit the VAE/PCA is then scored via 5-fold CV below). This mirrors the
     notebook's protocol; it is more optimistic than a fully held-out fold and is
     not a claim of generalization performance, just a like-for-like comparison.
  3. Score three representations - CNC-VAE latent, PCA (same latent size), and the
     raw concatenated features - with 5-fold StratifiedKFold cross_val_score for
     NB / SVM / RF, reporting accuracy, balanced_accuracy and (one-vs-rest macro)
     roc_auc.

Note: exact numbers won't match the PyTorch notebook bit-for-bit (different
framework, different weight initialization scheme), but should land in the same
ballpark under the same architecture/hyperparameters/evaluation protocol.
"""
import argparse
import os
import warnings

warnings.filterwarnings('ignore')

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.manifold import TSNE
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC

from misc.dataset import DEFAULT_AAK_FILE, DEFAULT_CLIN_FILE, DEFAULT_GENE_INPUT_FILE, get_data
from models.cncvae import CNCVAE

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--ds', type=int, default=256, help='Dense layer size (notebook default: 256)')
parser.add_argument('--ls', type=int, default=16, help='Latent size (notebook default: 16)')
parser.add_argument('--dropout', type=float, default=0.2)
parser.add_argument('--distance', type=str, default='mmd', choices=['kl', 'mmd'])
parser.add_argument('--target_ratio', type=float, default=1.0,
                     help='Weighted-distance : reconstruction ratio at epoch 0 - beta is '
                          'auto-resolved to hit this ratio, matching the notebook\'s '
                          '"init" balance_mode')
parser.add_argument('--epochs', type=int, default=150)
parser.add_argument('--bs', type=int, default=16, help='Batch size (notebook default: 16)')
parser.add_argument('--gene_input_file', type=str, default=DEFAULT_GENE_INPUT_FILE)
parser.add_argument('--aak_file', type=str, default=DEFAULT_AAK_FILE)
parser.add_argument('--clin_file', type=str, default=DEFAULT_CLIN_FILE)
parser.add_argument('--clinical_mode', type=str, default='curated', choices=['full', 'curated'])
parser.add_argument('--seed', type=int, default=0, help='Notebook uses SEED=0')
parser.add_argument('--out', type=str, default='comparison')

GRADE_COLORS = {0: '#3b6fb6', 1: '#4daf4a', 2: '#ff7f00', 3: '#c0392b'}


def _mmd_numpy(x, y):
    """Same 2*dim-bandwidth RBF/MMD as models.common.mmd, computed on plain arrays
    for the pre-training beta-resolution probe (no need for a TF graph here)."""
    def kernel(a, b):
        dim = a.shape[1]
        dist = np.sum((a[:, None, :] - b[None, :, :]) ** 2, axis=2)
        return np.exp(-dist / (2.0 * dim))
    return float(kernel(x, x).mean() + kernel(y, y).mean() - 2 * kernel(x, y).mean())


def _kl_numpy(z_mean, z_log_sigma):
    """Same KL-divergence formula as models.common.kl_regu, computed on plain
    numpy arrays (no TF graph needed for this one-off pre-training measurement)."""
    return float(np.mean(-0.5 * np.sum(1 + z_log_sigma - z_mean ** 2 - np.exp(z_log_sigma), axis=1)))


def resolve_beta(args, X, target_ratio):
    """Build a probe model (untrained, fixed-seed init) to measure the reconstruction
    and distance terms at epoch 0, then solve for the beta that makes the weighted
    distance term start at `target_ratio` * reconstruction - matching the notebook's
    resolve_beta()/balance_mode='init'.

    Why: reconstruction loss and the KL/MMD distance term live on very
    different numeric scales, so a single fixed beta that works for one
    dataset/config can make one term totally dominate on another. This solves
    for beta automatically instead of requiring manual per-run tuning.
    """
    # Build a throwaway model with the same architecture (beta itself isn't
    # used here) purely to see what it outputs right after weight init, i.e.
    # before any training has happened ("epoch 0").
    probe_args = argparse.Namespace(**vars(args))
    probe_args.beta = 0.0
    probe_args.act = 'elu'
    probe_args.input_size = X.shape[1]
    probe_args.save_model = False

    probe = CNCVAE(probe_args)
    probe.build_model()

    z_mean, z_log_sigma, z = probe.encoder.predict(X, batch_size=args.bs, verbose=0)
    x_hat = probe.vae.predict(X, batch_size=args.bs, verbose=0)
    recon0 = float(np.mean((X - x_hat) ** 2))

    if args.distance == 'mmd':
        rng = np.random.default_rng(args.seed)
        prior = rng.standard_normal(z.shape).astype(np.float32)
        distance0 = _mmd_numpy(prior, z)
    else:
        distance0 = _kl_numpy(z_mean, z_log_sigma)

    # Solve for beta such that beta * distance0 == target_ratio * recon0.
    beta_eff = target_ratio * recon0 / (distance0 + 1e-12)
    print('init-balance: recon_0={:.4f}, {}_0={:.4f}  ->  beta_eff={:.2f} '
          '(weighted {} starts at {:g}x reconstruction)'.format(
              recon0, args.distance.upper(), distance0, beta_eff, args.distance.upper(), target_ratio))
    return beta_eff


def classifiers(n_features, seed):
    """The 3 simple classifiers used to score each representation - kept
    simple/fast since this is a sanity check of separability, not a model to deploy."""
    return {
        'NB': GaussianNB(),
        'SVM': SVC(kernel='rbf', C=1.5, gamma=1.0 / n_features, probability=True, random_state=seed),
        'RF': RandomForestClassifier(n_estimators=50, max_features=0.5, random_state=seed),
    }


def cv_scores(features, target, label, seed):
    """Score one representation with 5-fold CV across all 3 classifiers.
    `target` here is the 4-class (0..3) MGS grade, so AUC uses one-vs-rest
    macro averaging rather than plain binary roc_auc."""
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    n_classes = len(np.unique(target))
    # AUC falls back to one-vs-rest macro AUC for the four-class grade target.
    auc_scoring = 'roc_auc' if n_classes == 2 else 'roc_auc_ovr'
    rows = []
    for name, clf in classifiers(features.shape[1], seed).items():
        acc = cross_val_score(clf, features, target, cv=cv, scoring='accuracy').mean()
        bal = cross_val_score(clf, features, target, cv=cv, scoring='balanced_accuracy').mean()
        auc = cross_val_score(clf, features, target, cv=cv, scoring=auc_scoring).mean()
        rows.append({'representation': label, 'model': name,
                     'accuracy': acc, 'balanced_acc': bal, 'roc_auc': auc})
    return rows


def tsne2d(features, seed):
    """2D t-SNE projection for visualization only - perplexity scaled to
    sample count since a fixed perplexity behaves poorly on small datasets."""
    perp = min(30, max(5, features.shape[0] // 5))
    return TSNE(n_components=2, perplexity=perp, init='pca', random_state=seed).fit_transform(features)


def main():
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    d = get_data(args.gene_input_file, args.aak_file, args.clin_file, clinical_mode=args.clinical_mode)
    y = d['y']
    print('rna:', d['rnanp'].shape, ' clin ({}):'.format(args.clinical_mode), d['clin'].shape,
          ' classes:', dict(zip(d['label_classes'], np.bincount(y))))

    # ----- mRNA source: log1p + min-max (matches the notebook) -----
    mrna_log = np.log1p(d['rnanp'].astype(np.float64))
    mrna_min, mrna_max = mrna_log.min(axis=0), mrna_log.max(axis=0)
    mrna_range = np.where(mrna_max - mrna_min == 0, 1.0, mrna_max - mrna_min)
    X_mrna = ((mrna_log - mrna_min) / mrna_range).astype(np.float32)

    # ----- clinical source: already numeric/curated by get_data() -----
    X_clin = d['clin'].astype(np.float32)

    X = np.hstack([X_mrna, X_clin]).astype(np.float32)
    print('concatenated input:', X.shape, '(mRNA {} + clinical {})'.format(X_mrna.shape[1], X_clin.shape[1]))

    # ----- train CNC-VAE on the whole cohort with an auto-balanced beta -----
    beta_eff = resolve_beta(args, X, args.target_ratio)

    args.beta = beta_eff
    args.act = 'elu'
    args.input_size = X.shape[1]
    args.save_model = False

    cncvae = CNCVAE(args)
    cncvae.build_model()
    cncvae.train(X_clin, X_mrna)  # concatenated internally as [clin, mrna] == X above
    Z = cncvae.predict(X_clin, X_mrna)
    print('latent embedding:', Z.shape)

    # ----- baselines -----
    Z_pca = PCA(n_components=args.ls, random_state=args.seed).fit_transform(X)

    # ----- 5-fold CV comparison table -----
    results = []
    results += cv_scores(Z, y, 'CNC-VAE latent', args.seed)
    results += cv_scores(Z_pca, y, 'PCA ({})'.format(args.ls), args.seed)
    results += cv_scores(X, y, 'Raw concatenated', args.seed)

    results_df = pd.DataFrame(results).set_index(['representation', 'model']).round(3)
    print(results_df)
    results_df.to_csv(os.path.join(args.out, 'comparison_table.csv'))

    # ----- 3-panel t-SNE -----
    panels = {'Raw concatenated': X, 'PCA ({})'.format(args.ls): Z_pca, 'CNC-VAE latent': Z}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    for ax, (title, feat) in zip(axes, panels.items()):
        emb = tsne2d(feat, args.seed)
        for cls_idx, cls_name in enumerate(d['label_classes']):
            mask = y == cls_idx
            ax.scatter(emb[mask, 0], emb[mask, 1], s=18, alpha=0.75,
                       label=cls_name, color=GRADE_COLORS.get(cls_idx))
        ax.set_title(title)
        ax.set_xlabel('t-SNE 1')
        ax.set_ylabel('t-SNE 2')
        ax.legend()
    fig.suptitle('Clinical + mRNA integration, colored by AMD grade', y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, 'tsne_comparison.png'), bbox_inches='tight')
    print('Saved', os.path.join(args.out, 'tsne_comparison.png'))

    # ----- save the latent embedding -----
    # mgs_class = integer label (0..3), mgs_grade = human-readable string
    # ('MGS1'..'MGS4') recovered via label_classes[y].
    out_df = pd.DataFrame({'sample_id': d['sample_id'], 'mgs_class': y})
    out_df['mgs_grade'] = d['label_classes'][y]
    for j in range(Z.shape[1]):
        out_df['z{}'.format(j)] = Z[:, j]
    out_path = os.path.join(args.out, 'cncvae_latent_embedding.csv')
    out_df.to_csv(out_path, index=False)
    print('Saved', out_path, out_df.shape)


if __name__ == '__main__':
    main()
