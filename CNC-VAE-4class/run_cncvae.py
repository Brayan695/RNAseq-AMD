# Command-line entry point for training the 4-class CNC-VAE: same idea as
# ../../CNC-VAE/run_cncvae.py (whole-cohort mode via --fold 0, or one
# stratified CV fold via --fold N), but trains on the 81-gene subset pulled
# from gene_input.csv (all 453 donors, all 4 MGS grades) instead of the
# smaller 166-sample MGS1/MGS4-only aak100_cpmdat.csv. For the notebook-
# matching workflow (auto-tuned beta, baseline comparison, t-SNE), see
# compare_representations.py instead.
import argparse
import os

import numpy as np

from misc.dataset import AMDDataset, DEFAULT_AAK_FILE, DEFAULT_CLIN_FILE, DEFAULT_GENE_INPUT_FILE
from misc.helpers import normalizeRNA, save_embedding
from models.cncvae import CNCVAE

# Preset (hidden-layer-size -> other hyperparameters) bundles, selected via --ds.
configs = {
    16: {'ds': 16, 'act': 'elu', 'epochs': 150, 'bs': 64, 'dropout': 0.2},
    32: {'ds': 32, 'act': 'elu', 'epochs': 150, 'bs': 64, 'dropout': 0.2},
    64: {'ds': 64, 'act': 'elu', 'epochs': 150, 'bs': 64, 'dropout': 0.2},
    128: {'ds': 128, 'act': 'elu', 'epochs': 150, 'bs': 64, 'dropout': 0.2},
    256: {'ds': 256, 'act': 'elu', 'epochs': 150, 'bs': 16, 'dropout': 0.2},
}

parser = argparse.ArgumentParser(
    description='Train a CNC-VAE integrating the 81-gene AMD RNA-seq subset (matched from '
                 'gene_input.csv against the ENSG ids in aak100_cpmdat.csv, covering all 453 '
                 'donors) with clinical/genotype covariates (MetaSheet_Processed.csv), for the '
                 'four-class MGS1..MGS4 AMD grade task. Mirrors the CNC_VAE_AMD_V5_4class notebook.'
)
parser.add_argument('--ds', help='Intermediate dense layer size: 16, 32, 64, 128 or 256', type=int,
                     choices=[16, 32, 64, 128, 256], required=True)
parser.add_argument('--ls', help='Latent dimension size', type=int, required=True)
parser.add_argument('--beta', help='Beta-VAE regularization weight', type=float, default=1.0)
parser.add_argument('--distance', help='Regularizer: kl or mmd', type=str, default='kl')
parser.add_argument('--fold', help="Fold number (1..numfolds) to train/validate on, "
                                    "or '0' to train on the whole cohort", type=str, default='0')
parser.add_argument('--numfolds', help='Number of stratified CV folds (ignored when --fold 0)',
                     type=int, default=5)
parser.add_argument('--gene_input_file', help='Path to gene_input.csv (all 453 donors, gene SYMBOLS)',
                     type=str, default=DEFAULT_GENE_INPUT_FILE)
parser.add_argument('--aak_file', help='Path to aak100_cpmdat.csv - used only to pick which 81 genes '
                                        'to take from gene_input_file (matched by value on shared samples)',
                     type=str, default=DEFAULT_AAK_FILE)
parser.add_argument('--clin_file', help='Path to the clinical/genotype metadata CSV',
                     type=str, default=DEFAULT_CLIN_FILE)
parser.add_argument('--clinical_mode', help="'full' (every clinical column) or 'curated' "
                                             "(age/sex/genotype + prevalence-filtered, leakage-excluded "
                                             "history flags - matches the reference notebook)",
                     type=str, default='curated', choices=['full', 'curated'])
parser.add_argument('--label_col', help='Phenotype column (used for CV stratification and saved '
                                         'alongside the embeddings for downstream analysis)',
                     type=str, default='mgs_level')
parser.add_argument('--save_model', help='Save the trained VAE weights', action='store_true')
parser.add_argument('--writedir', help='Output directory - default is ./results', type=str, default='')
parser.add_argument('--seed', help='Random seed for the stratified CV split', type=int, default=42)

if __name__ == '__main__':
    args = parser.parse_args()
    config = configs[args.ds]
    for key, val in config.items():
        setattr(args, key, val)

    dataset = AMDDataset(
        args.gene_input_file, args.aak_file, args.clin_file, args.label_col,
        n_splits=args.numfolds, seed=args.seed, clinical_mode=args.clinical_mode,
    )

    root_dir = args.writedir or 'results'
    emb_save_dir = os.path.join(
        root_dir, 'CNCVAE_Clin+mRNA_integration',
        'cncvae_LS_{}_DS_{}_{}_beta_{}'.format(args.ls, args.ds, args.distance, args.beta),
    )
    os.makedirs(emb_save_dir, exist_ok=True)

    if args.fold == '0':
        # ----- Whole-cohort mode: train on every sample, no held-out test set -----
        print('TRAINING on the complete 453-donor AMD cohort (4-class)')

        whole = dataset.whole
        s1_train = whole['clin']
        # log1p compresses RNA-seq's long right tail before scaling to [0, 1].
        s2_train = normalizeRNA(np.log1p(whole['rnanp']))

        args.input_size = s1_train.shape[1] + s2_train.shape[1]
        args.model_out = os.path.join(emb_save_dir, 'vae_cncvae.weights.h5')

        cncvae = CNCVAE(args)
        cncvae.build_model()
        cncvae.train(s1_train, s2_train)
        emb_train = cncvae.predict(s1_train, s2_train)

        save_embedding(emb_save_dir, args.label_col + '.npz', emb_train)
        np.savez(
            os.path.join(emb_save_dir, args.label_col + '_labels.npz'),
            y=whole['y'], sample_id=whole['sample_id'], classes=whole['label_classes'],
        )
    else:
        # ----- Fold mode: train on one stratified CV split, evaluate held-out fold -----
        fold = int(args.fold)
        print('TRAINING on fold {}'.format(fold))

        train, test = dataset.fold(fold)
        s1_train, s1_test = train['clin'], test['clin']
        s2_train, s2_test = normalizeRNA(np.log1p(train['rnanp']), np.log1p(test['rnanp']))

        args.input_size = s1_train.shape[1] + s2_train.shape[1]
        args.model_out = os.path.join(emb_save_dir, 'vae_cncvae_fold{}.weights.h5'.format(fold))

        cncvae = CNCVAE(args)
        cncvae.build_model()
        cncvae.train(s1_train, s2_train, s1_test, s2_test)

        emb_train = cncvae.predict(s1_train, s2_train)
        emb_test = cncvae.predict(s1_test, s2_test)

        save_embedding(emb_save_dir, args.label_col + str(fold) + '.npz', emb_train, emb_test)
        np.savez(
            os.path.join(emb_save_dir, args.label_col + str(fold) + '_labels.npz'),
            y_train=train['y'], y_test=test['y'],
            sample_id_train=train['sample_id'], sample_id_test=test['sample_id'],
        )
