import argparse
import os

import numpy as np

from misc.dataset import AMDDataset, DEFAULT_CLIN_FILE, DEFAULT_RNA_FILE
from misc.helpers import normalizeRNA, save_embedding
from models.cncvae import CNCVAE

configs = {
    16: {'ds': 16, 'act': 'elu', 'epochs': 150, 'bs': 64, 'dropout': 0.2},
    32: {'ds': 32, 'act': 'elu', 'epochs': 150, 'bs': 64, 'dropout': 0.2},
    64: {'ds': 64, 'act': 'elu', 'epochs': 150, 'bs': 64, 'dropout': 0.2},
    128: {'ds': 128, 'act': 'elu', 'epochs': 150, 'bs': 64, 'dropout': 0.2},
}

parser = argparse.ArgumentParser(
    description='Train a CNC-VAE integrating AMD RNA-seq CPM expression with clinical/genotype '
                 'covariates (MetaSheet_1_4.csv), matched by sample_id.'
)
parser.add_argument('--ds', help='Intermediate dense layer size: 16, 32, 64 or 128', type=int,
                     choices=[16, 32, 64, 128], required=True)
parser.add_argument('--ls', help='Latent dimension size', type=int, required=True)
parser.add_argument('--beta', help='Beta-VAE regularization weight', type=float, default=1.0)
parser.add_argument('--distance', help='Regularizer: kl or mmd', type=str, default='kl')
parser.add_argument('--fold', help="Fold number (1..numfolds) to train/validate on, "
                                    "or '0' to train on the whole cohort", type=str, default='0')
parser.add_argument('--numfolds', help='Number of stratified CV folds (ignored when --fold 0)',
                     type=int, default=5)
parser.add_argument('--rna_file', help='Path to the gene-expression/CPM CSV '
                                        '(swap this in for a full-transcriptome file later)',
                     type=str, default=DEFAULT_RNA_FILE)
parser.add_argument('--clin_file', help='Path to the clinical/genotype metadata CSV',
                     type=str, default=DEFAULT_CLIN_FILE)
parser.add_argument('--clinical_mode', help="'full' (every clinical column) or 'curated' "
                                             "(age/sex/genotype + prevalence-filtered, leakage-excluded "
                                             "history flags - matches the CNC_VAE_AMD_V3 reference notebook)",
                     type=str, default='full', choices=['full', 'curated'])
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
        args.rna_file, args.clin_file, args.label_col, n_splits=args.numfolds, seed=args.seed,
        clinical_mode=args.clinical_mode,
    )

    root_dir = args.writedir or 'results'
    emb_save_dir = os.path.join(
        root_dir, 'CNCVAE_Clin+mRNA_integration',
        'cncvae_LS_{}_DS_{}_{}_beta_{}'.format(args.ls, args.ds, args.distance, args.beta),
    )
    os.makedirs(emb_save_dir, exist_ok=True)

    if args.fold == '0':
        print('TRAINING on the complete AMD cohort')

        whole = dataset.whole
        s1_train = whole['clin']
        s2_train = normalizeRNA(whole['rnanp'])

        args.input_size = s1_train.shape[1] + s2_train.shape[1]
        args.model_out = os.path.join(emb_save_dir, 'vae_cncvae.h5')

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
        fold = int(args.fold)
        print('TRAINING on fold {}'.format(fold))

        train, test = dataset.fold(fold)
        s1_train, s1_test = train['clin'], test['clin']
        s2_train, s2_test = normalizeRNA(train['rnanp'], test['rnanp'])

        args.input_size = s1_train.shape[1] + s2_train.shape[1]
        args.model_out = os.path.join(emb_save_dir, 'vae_cncvae_fold{}.h5'.format(fold))

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
