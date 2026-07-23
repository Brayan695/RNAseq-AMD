# Loads a previously-trained X-DEC (xvae) encoder (no retraining) and trains
# a small binary classifier head on top of its FROZEN latent embedding, to
# predict mgs_level (MGS1 vs MGS4) directly from the compressed representation.
#
# The classifier head's architecture (models/classifier_head.py) is picked
# automatically from the loaded encoder's own latent size (--ls the checkpoint
# was trained with): 16->4->1, 8->4->1, or plain 4->1. Only the head's weights
# are ever trained - the encoder is used purely for inference (ae.predict()).
#
# Since there are only ~166 samples total, a single train/test split would be
# too noisy to trust (see utils/train_evaluate.py's own Bayesian-search CV fix
# elsewhere in this project for exactly this lesson) - so the head is trained
# and evaluated via stratified 5-fold CV, with metrics computed on the
# concatenated out-of-fold predictions (every sample gets exactly one
# held-out prediction across the whole CV) rather than averaging per-fold
# metrics, which would be noisy with ~33 samples per fold.
import argparse
import json
import os

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import torch
import torch.nn as nn
from sklearn.metrics import (
    average_precision_score, confusion_matrix, precision_recall_curve, roc_auc_score, roc_curve,
)
from sklearn.model_selection import StratifiedKFold

from misc.dataset import DEFAULT_CLIN_FILE, DEFAULT_RNA_FILE, get_data
from misc.helpers import normalizeRNA
from models.xvae import load_xvae_model
from models.classifier_head import ClassifierHead

parser = argparse.ArgumentParser(
    description='Train a binary classifier head on top of a frozen, already-trained X-DEC (xvae) '
                "encoder's latent embedding, evaluated via stratified k-fold CV "
                '(PR-AUC / ROC-AUC / confusion matrix).'
)
parser.add_argument('--encoder', required=True, help='Path to the saved encoder .pt file')
parser.add_argument('--clinical_mode', type=str, default='curated', choices=['full', 'curated'],
                     help='Must match what the checkpoint was originally trained with')
parser.add_argument('--rna_file', type=str, default=DEFAULT_RNA_FILE)
parser.add_argument('--clin_file', type=str, default=DEFAULT_CLIN_FILE)
parser.add_argument('--label_col', type=str, default='mgs_level')
parser.add_argument('--n_splits', type=int, default=5, help='Number of stratified CV folds')
parser.add_argument('--epochs', type=int, default=200, help='Classifier-head training epochs per fold')
parser.add_argument('--bs', type=int, default=32, help='Classifier-head training batch size')
parser.add_argument('--lr', type=float, default=1e-3)
parser.add_argument('--seed', type=int, default=5192)
parser.add_argument('--out', type=str, default='classifier_results')


def train_head(z_train, y_train, latent_dim, epochs, bs, lr, seed):
    """Trains a FRESH ClassifierHead from scratch on one fold's train split -
    a new head every fold, so no information leaks across folds."""
    torch.manual_seed(seed)
    head = ClassifierHead(latent_dim)
    optimizer = torch.optim.Adam(head.parameters(), lr=lr)
    loss_fn = nn.BCELoss()

    z_t = torch.tensor(z_train, dtype=torch.float32)
    y_t = torch.tensor(y_train, dtype=torch.float32)
    n = z_t.shape[0]

    head.train()
    for _ in range(epochs):
        perm = torch.randperm(n)
        for i in range(0, n, bs):
            idx = perm[i:i + bs]
            optimizer.zero_grad()
            proba = head(z_t[idx])
            loss = loss_fn(proba, y_t[idx])
            loss.backward()
            optimizer.step()
    return head


def main():
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    ae = load_xvae_model(args.encoder)  # build_model() + load_state_dict() + eval() already done
    latent_dim = ae.args.ls
    print('Loaded encoder: ds1={} ds2={} ds12={} ls={} s1_input_size={} s2_input_size={}'.format(
        ae.args.ds1, ae.args.ds2, ae.args.ds12, latent_dim, ae.args.s1_input_size, ae.args.s2_input_size))

    d = get_data(args.rna_file, args.clin_file, args.label_col, clinical_mode=args.clinical_mode)
    x_num = normalizeRNA(d['rnanp'])   # SET 1: numerical RNA-seq CPM
    x_bin = d['clin']                  # SET 2: binary clinical/genotype dummies

    if x_num.shape[1] != ae.args.s1_input_size or x_bin.shape[1] != ae.args.s2_input_size:
        raise ValueError(
            'Input widths (numeric={}, binary={}) do not match checkpoint '
            '(s1_input_size={}, s2_input_size={}) - check --clinical_mode matches '
            'how this encoder was trained'.format(
                x_num.shape[1], x_bin.shape[1], ae.args.s1_input_size, ae.args.s2_input_size)
        )

    z = ae.predict(x_num, x_bin, output='encoder')  # frozen - predict() already wraps eval()+no_grad
    y = d['y']
    print('latent embedding:', z.shape, ' class balance:', dict(zip(*np.unique(y, return_counts=True))))

    skf = StratifiedKFold(n_splits=args.n_splits, shuffle=True, random_state=args.seed)
    oof_proba = np.zeros(len(y))
    for fold, (train_idx, test_idx) in enumerate(skf.split(z, y)):
        head = train_head(z[train_idx], y[train_idx], latent_dim, args.epochs, args.bs, args.lr,
                           args.seed + fold)
        head.eval()
        with torch.no_grad():
            oof_proba[test_idx] = head(torch.tensor(z[test_idx], dtype=torch.float32)).numpy()
        print('Fold {}/{} done'.format(fold + 1, args.n_splits))

    y_pred = (oof_proba >= 0.5).astype(int)
    pr_auc = average_precision_score(y, oof_proba)
    roc_auc = roc_auc_score(y, oof_proba)
    cm = confusion_matrix(y, y_pred)
    print('PR-AUC={:.4f}  ROC-AUC={:.4f}'.format(pr_auc, roc_auc))
    print('Confusion matrix:\n', cm)

    precision, recall, _ = precision_recall_curve(y, oof_proba)
    fpr, tpr, _ = roc_curve(y, oof_proba)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    axes[0].plot(recall, precision)
    axes[0].set_xlabel('Recall'); axes[0].set_ylabel('Precision')
    axes[0].set_title('PR curve (AP={:.3f})'.format(pr_auc))

    axes[1].plot(fpr, tpr)
    axes[1].plot([0, 1], [0, 1], '--', color='gray')
    axes[1].set_xlabel('False positive rate'); axes[1].set_ylabel('True positive rate')
    axes[1].set_title('ROC curve (AUC={:.3f})'.format(roc_auc))

    im = axes[2].imshow(cm, cmap='Blues')
    classes = d['label_classes']
    axes[2].set_xticks(range(len(classes))); axes[2].set_xticklabels(classes)
    axes[2].set_yticks(range(len(classes))); axes[2].set_yticklabels(classes)
    axes[2].set_xlabel('Predicted'); axes[2].set_ylabel('True')
    for i in range(cm.shape[0]):
        for j in range(cm.shape[1]):
            axes[2].text(j, i, cm[i, j], ha='center', va='center')
    axes[2].set_title('Confusion matrix')
    fig.colorbar(im, ax=axes[2], fraction=0.046)

    fig.suptitle('X-DEC classifier head (latent_dim={}) - {}-fold CV, out-of-fold predictions'.format(
        latent_dim, args.n_splits))
    fig.tight_layout()
    plot_path = os.path.join(args.out, 'xdec_classifier_curves.png')
    fig.savefig(plot_path, bbox_inches='tight')
    print('Saved', plot_path)

    result_path = os.path.join(args.out, 'xdec_classifier_result.json')
    with open(result_path, 'w') as f:
        json.dump({
            'model': 'X-DEC', 'encoder': args.encoder, 'latent_dim': int(latent_dim),
            'n_splits': args.n_splits, 'pr_auc': float(pr_auc), 'roc_auc': float(roc_auc),
            'confusion_matrix': cm.tolist(), 'n_samples': int(len(y)),
        }, f, indent=2)
    print('Saved', result_path)


if __name__ == '__main__':
    main()
