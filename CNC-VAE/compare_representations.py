"""Reproduce the comparison table / t-SNE / embedding export from the CNC_VAE_AMD_V3
reference notebook, using this project's PyTorch CNC-VAE port, so the two can be
compared side by side.

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
     NB / SVM / RF, reporting accuracy, balanced_accuracy and roc_auc.

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
import torch
import torch.nn.functional as F
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.manifold import TSNE
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC

from misc.dataset import DEFAULT_CLIN_FILE, DEFAULT_RNA_FILE, get_data
from models.common import kl_regu, mmd
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
parser.add_argument('--rna_file', type=str, default=DEFAULT_RNA_FILE)
parser.add_argument('--clin_file', type=str, default=DEFAULT_CLIN_FILE)
parser.add_argument('--clinical_mode', type=str, default='curated', choices=['full', 'curated'])
parser.add_argument('--seed', type=int, default=0, help='Notebook uses SEED=0')
parser.add_argument('--out', type=str, default='comparison')
parser.add_argument('--save_model', action='store_true',
                     help='Save the trained encoder (encoder_cncvae.pt) to --out, e.g. for SHAP')

CLASS_COLORS = {0: '#3b6fb6', 1: '#c0392b'}


def resolve_beta(args, X, target_ratio):
    """Build a probe model (untrained, fixed-seed init) to measure the reconstruction
    and distance terms at epoch 0, then solve for the beta that makes the weighted
    distance term start at `target_ratio` * reconstruction - matching the notebook's
    resolve_beta()/balance_mode='init'.

    Why this matters: reconstruction loss and the KL/MMD distance term are on
    very different numeric scales, so a fixed beta that works for one dataset
    can make the other term completely dominate the loss on a different one.
    This picks beta automatically so the two terms start out roughly balanced,
    rather than requiring manual tuning per dataset/config.
    """
    # Build a throwaway model with the same architecture but beta=0 (not used
    # for this probe - we're only interested in what it outputs before any
    # training has happened, i.e. right after weight initialization).
    probe_args = argparse.Namespace(**vars(args))
    probe_args.beta = 0.0
    probe_args.act = 'elu'
    probe_args.input_size = X.shape[1]
    probe_args.save_model = False

    probe = CNCVAE(probe_args)
    probe.build_model()

    # One forward pass over the whole dataset with the freshly-initialized
    # (untrained) weights, to measure both loss terms at "epoch 0".
    probe.net.eval()
    x_t = torch.tensor(X, dtype=torch.float32)
    with torch.no_grad():
        x_hat, z_mean, z_log_sigma, z = probe.net(x_t)
    recon0 = float(F.mse_loss(x_hat, x_t, reduction='mean').item())

    if args.distance == 'mmd':
        rng = np.random.default_rng(args.seed)
        prior = torch.tensor(rng.standard_normal(z.shape).astype(np.float32))
        distance0 = float(mmd(prior, z).item())
    else:
        distance0 = float(kl_regu(z_mean, z_log_sigma).mean().item())

    # Solve for beta such that beta * distance0 == target_ratio * recon0, i.e.
    # the weighted distance term starts at target_ratio times the reconstruction
    # term (target_ratio=1.0 means "start balanced").
    beta_eff = target_ratio * recon0 / (distance0 + 1e-12)
    print('init-balance: recon_0={:.4f}, {}_0={:.4f}  ->  beta_eff={:.2f} '
          '(weighted {} starts at {:g}x reconstruction)'.format(
              recon0, args.distance.upper(), distance0, beta_eff, args.distance.upper(), target_ratio))
    return beta_eff


def classifiers(n_features, seed):
    """The 3 simple classifiers used to score each representation (raw features,
    PCA, CNC-VAE latent) - kept simple and fast since this is just a sanity
    check of how separable the classes are, not a model to deploy."""
    return {
        'NB': GaussianNB(),
        'SVM': SVC(kernel='rbf', C=1.5, gamma=1.0 / n_features, probability=True, random_state=seed),
        'RF': RandomForestClassifier(n_estimators=50, max_features=0.5, random_state=seed),
    }


def cv_scores(features, target, label, seed):
    """Score one representation (a 2D feature matrix) with 5-fold CV across all
    3 classifiers, returning one result row per classifier with accuracy,
    balanced accuracy, and ROC AUC."""
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    rows = []
    for name, clf in classifiers(features.shape[1], seed).items():
        acc = cross_val_score(clf, features, target, cv=cv, scoring='accuracy').mean()
        bal = cross_val_score(clf, features, target, cv=cv, scoring='balanced_accuracy').mean()
        auc = cross_val_score(clf, features, target, cv=cv, scoring='roc_auc').mean()
        rows.append({'representation': label, 'model': name,
                     'accuracy': acc, 'balanced_acc': bal, 'roc_auc': auc})
    return rows


def tsne2d(features, seed):
    """2D t-SNE projection for visualization only (not used anywhere else in
    the pipeline) - perplexity is scaled to the sample count since t-SNE
    behaves poorly with a fixed perplexity on very small datasets."""
    perp = min(30, max(5, features.shape[0] // 5))
    return TSNE(n_components=2, perplexity=perp, init='pca', random_state=seed).fit_transform(features)


def main():
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    d = get_data(args.rna_file, args.clin_file, clinical_mode=args.clinical_mode)
    y = d['y']
    print('rna:', d['rnanp'].shape, ' clin ({}):'.format(args.clinical_mode), d['clin'].shape,
          ' classes:', dict(zip(d['label_classes'], np.bincount(y))))

    # ----- mRNA source: log1p + min-max (matches the notebook) -----
    # log1p compresses the long right tail typical of RNA-seq CPM values (a few
    # highly-expressed genes would otherwise dominate the reconstruction loss).
    mrna_log = np.log1p(d['rnanp'].astype(np.float64))
    mrna_min, mrna_max = mrna_log.min(axis=0), mrna_log.max(axis=0)
    mrna_range = np.where(mrna_max - mrna_min == 0, 1.0, mrna_max - mrna_min)
    X_mrna = ((mrna_log - mrna_min) / mrna_range).astype(np.float32)

    # ----- clinical source: already numeric/curated by get_data() -----
    X_clin = d['clin'].astype(np.float32)

    # This is the actual model input: clinical + gene-expression features
    # concatenated into one vector per sample ("CNC" = Concatenated iNputs).
    X = np.hstack([X_mrna, X_clin]).astype(np.float32)
    print('concatenated input:', X.shape, '(mRNA {} + clinical {})'.format(X_mrna.shape[1], X_clin.shape[1]))

    # ----- train CNC-VAE on the whole cohort with an auto-balanced beta -----
    beta_eff = resolve_beta(args, X, args.target_ratio)

    save_model = args.save_model
    args.beta = beta_eff
    args.act = 'elu'
    args.input_size = X.shape[1]
    args.save_model = False  # suppress CNCVAE.train()'s own vae.save_weights (no model_out set here);
                              # the encoder alone is saved separately below when requested.

    cncvae = CNCVAE(args)
    cncvae.build_model()
    cncvae.train(X_clin, X_mrna)  # concatenated internally as [clin, mrna] == X above
    Z = cncvae.predict(X_clin, X_mrna)  # (n_samples, ls) latent embedding
    print('latent embedding:', Z.shape)

    if save_model:
        encoder_path = os.path.join(args.out, 'encoder_cncvae.pt')
        cncvae.save_encoder(encoder_path)
        print('Saved encoder to', encoder_path)

    # ----- baselines -----
    # PCA at the same output dimensionality as the VAE's latent size - the
    # question this whole script answers is "does the VAE's learned embedding
    # actually beat a much simpler linear compression (PCA), or the raw
    # features with no compression at all?"
    Z_pca = PCA(n_components=args.ls, random_state=args.seed).fit_transform(X)

    # ----- 5-fold CV comparison table -----
    # Same downstream classifiers, same CV splits, applied to all three
    # representations (CNC-VAE latent / PCA / raw) so the comparison is
    # apples-to-apples.
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
                       label=cls_name, color=CLASS_COLORS.get(cls_idx))
        ax.set_title(title)
        ax.set_xlabel('t-SNE 1')
        ax.set_ylabel('t-SNE 2')
        ax.legend()
    fig.suptitle('Clinical + mRNA integration, colored by AMD grade', y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, 'tsne_comparison.png'), bbox_inches='tight')
    print('Saved', os.path.join(args.out, 'tsne_comparison.png'))

    # ----- save the latent embedding -----
    out_df = pd.DataFrame({'sample_id': d['sample_id'], 'mgs_bin': y})
    out_df['mgs_grade'] = d['label_classes'][y]
    for j in range(Z.shape[1]):
        out_df['z{}'.format(j)] = Z[:, j]
    out_path = os.path.join(args.out, 'cncvae_latent_embedding.csv')
    out_df.to_csv(out_path, index=False)
    print('Saved', out_path, out_df.shape)


if __name__ == '__main__':
    main()
