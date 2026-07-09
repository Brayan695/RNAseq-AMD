import argparse
import os

import numpy as np
import pandas as pd

from misc.dataset import get_data, DEFAULT_CLIN_FILE, DEFAULT_RNA_FILE
from misc.helpers import normalizeRNA, save_results
from models.cncvae import CNCVAE
from models.dec import DEC
from models.metrics import calculate_metrics

parser = argparse.ArgumentParser(
    description='Train a CNC-DEC model (CNC-VAE + Deep Embedded Clustering) on the AMD cohort. '
                 'CNC-VAE is the single concatenated-input variational encoder from ../CNC-VAE '
                 '(KL/MMD-regularized, MSE reconstruction over RNA-seq + clinical features '
                 'concatenated together - no type-specific branches, unlike X-DEC). A DEC '
                 'self-training clustering head (port of DAM-IC DECLayer.py + '
                 'clustering_functions.py) is fine-tuned jointly with the encoder on top of its '
                 'latent space.'
)
parser.add_argument('--ds', help='Encoder/decoder hidden layer size', type=int, default=64)
parser.add_argument('--ls', help='Latent dimension size', type=int, default=16)
parser.add_argument('--act', help='Activation function for hidden layers', type=str, default='elu')
parser.add_argument('--dropout', help='Decoder dropout rate', type=float, default=0.2)
parser.add_argument('--distance', help='VAE regularizer: kl or mmd', type=str, default='kl')
parser.add_argument('--beta', help='Beta-VAE regularization weight', type=float, default=1.0)
parser.add_argument('--epochs', help='CNC-VAE pretraining epochs', type=int, default=150)
parser.add_argument('--bs', help='CNC-VAE pretraining batch size', type=int, default=32)
parser.add_argument('--n_clusters', help='Number of clusters to identify', type=int, default=2)
parser.add_argument('--dec_maxiter', help='Max DEC self-training iterations', type=int, default=2000)
parser.add_argument('--dec_update_interval', help='Iterations between target-distribution updates',
                     type=int, default=50)
parser.add_argument('--dec_batch_size', help='DEC self-training batch size', type=int, default=64)
parser.add_argument('--dec_tol', help='Stop when the fraction of samples changing cluster '
                                       'assignment drops below this', type=float, default=0.01)
parser.add_argument('--n_init', help='Number of K-means initializations for cluster centers',
                     type=int, default=20)
parser.add_argument('--stability', help='Also run a repeated K-fold cluster-stability check '
                                         '(retrains the model k*rep times - slow)',
                     action='store_true')
parser.add_argument('--stability_k', help='Stability check: number of K-fold splits', type=int, default=5)
parser.add_argument('--stability_rep', help='Stability check: number of repeats', type=int, default=2)
parser.add_argument('--rna_file', help='Path to the gene-expression/CPM CSV', type=str,
                     default=DEFAULT_RNA_FILE)
parser.add_argument('--clin_file', help='Path to the clinical/genotype metadata CSV', type=str,
                     default=DEFAULT_CLIN_FILE)
parser.add_argument('--clinical_mode', help="'curated' (age/sex/genotype + prevalence-filtered, "
                                             "leakage-excluded history flags - recommended for "
                                             "clustering on this cohort size) or 'full' (every "
                                             "clinical column)", type=str, default='curated',
                     choices=['full', 'curated'])
parser.add_argument('--label_col', help='Phenotype column (used only for reporting ARI/NMI/'
                                         'accuracy against the DEC clusters, never for training)',
                     type=str, default='mgs_level')
parser.add_argument('--save_model', help='Save the trained CNC-VAE encoder (for the SHAP workflow)',
                     action='store_true')
parser.add_argument('--writedir', help='Output directory - default is ./results', type=str, default='')
parser.add_argument('--seed', help='Random seed', type=int, default=5192)


def run_cncdec(x_num, x_bin, args, y=None):
    """Pretrains the CNC-VAE on the concatenated input, then runs DEC
    self-training clustering on top of its latent space."""
    input_size = x_num.shape[1] + x_bin.shape[1]

    autoencoder = CNCVAE(
        input_size=input_size, ds=args.ds, ls=args.ls, act=args.act, dropout=args.dropout,
        distance=args.distance, beta=args.beta, epochs=args.epochs, bs=args.bs,
    )
    autoencoder.build_model(seed=args.seed)
    autoencoder.train(x_num, x_bin, seed=args.seed)

    dec = DEC(autoencoder, n_clusters=args.n_clusters, seed=args.seed)
    y_pred, y_proba, centroids, z = dec.fit(
        x_num, x_bin, y=y, maxiter=args.dec_maxiter, batch_size=args.dec_batch_size,
        update_interval=args.dec_update_interval, tol=args.dec_tol, n_init=args.n_init,
    )
    return autoencoder, dec, y_pred, y_proba, centroids, z


if __name__ == '__main__':
    args = parser.parse_args()

    data = get_data(args.rna_file, args.clin_file, args.label_col, clinical_mode=args.clinical_mode)

    x_num = normalizeRNA(data['rnanp'])  # RNA-seq CPM, min-max scaled
    x_bin = data['clin']                 # binary clinical/genotype dummies
    print('CNC-DEC input: {} samples, {} genes + {} clinical = {} concatenated features'.format(
        x_num.shape[0], x_num.shape[1], x_bin.shape[1], x_num.shape[1] + x_bin.shape[1]))

    root_dir = args.writedir or 'results'
    run_dir = os.path.join(
        root_dir, 'CNCDEC_Clin+mRNA_integration',
        'cncdec_LS_{}_DS_{}_{}_beta_{}_k{}'.format(args.ls, args.ds, args.distance, args.beta,
                                                    args.n_clusters),
    )
    os.makedirs(run_dir, exist_ok=True)

    autoencoder, dec, y_pred, y_proba, centroids, z = run_cncdec(x_num, x_bin, args, y=data['y'])

    if data['y'] is not None:
        acc, ari, nmi = calculate_metrics(data['y'], y_pred)
        print('Final clustering vs. {}: acc={:.5f}, ari={:.5f}, nmi={:.5f}'.format(
            args.label_col, acc, ari, nmi))

    if args.save_model:
        encoder_path = os.path.join(run_dir, 'encoder_cncvae.pt')
        autoencoder.save_encoder(encoder_path)
        print('Saved encoder to', encoder_path)

    save_results(
        run_dir, args.label_col + '.npz',
        z=z, y_pred=y_pred, y_proba=y_proba, centroids=centroids,
        y=data['y'], sample_id=data['sample_id'], classes=data['label_classes'],
    )
    print('Saved results to', os.path.join(run_dir, args.label_col + '.npz'))

    # ----- latent embedding CSV, drop-in shape for the existing XGBoost pipeline -----
    # (Bayesian Hyperparameter Tuning/utils/load_preprocess.py expects 'sample_id' and
    # a 'mgs_level' column containing the literal MGS1/MGS4 strings.)
    emb_df = pd.DataFrame({'sample_id': data['sample_id']})
    for j in range(z.shape[1]):
        emb_df['z{}'.format(j)] = z[:, j]
    emb_df['mgs_level'] = data['label_classes'][data['y']]
    emb_path = os.path.join(run_dir, 'cncdec_latent_embedding.csv')
    emb_df.to_csv(emb_path, index=False)
    print('Saved latent embedding to', emb_path, emb_df.shape)

    # ----- optional cluster-stability check -----
    if args.stability:
        from models.stability import compute_cluster_stability

        print('Running cluster-stability check (k={}, rep={}, {} total resamples)...'.format(
            args.stability_k, args.stability_rep, args.stability_k * args.stability_rep))
        cncvae_kwargs = dict(ds=args.ds, ls=args.ls, act=args.act, dropout=args.dropout,
                             distance=args.distance, beta=args.beta, epochs=args.epochs, bs=args.bs)
        dec_kwargs = dict(maxiter=args.dec_maxiter, batch_size=args.dec_batch_size,
                          update_interval=args.dec_update_interval, tol=args.dec_tol,
                          n_init=args.n_init)
        stability, y_pred_mat = compute_cluster_stability(
            x_num, x_bin, args.n_clusters, reference_labels=y_pred,
            cncvae_kwargs=cncvae_kwargs, dec_kwargs=dec_kwargs,
            k=args.stability_k, rep=args.stability_rep, seed=args.seed,
        )
        stability_path = os.path.join(run_dir, 'cluster_stability.csv')
        pd.DataFrame({'sample_id': data['sample_id'], 'cluster': y_pred,
                      'stability': stability.to_numpy()}).to_csv(stability_path, index=False)
        print('Mean cluster stability: {:.3f}'.format(stability.mean()))
        print('Saved', stability_path)
