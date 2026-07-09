import argparse
import os

import numpy as np
import pandas as pd

from misc.dataset import get_data, DEFAULT_CLIN_FILE, DEFAULT_RNA_FILE
from misc.helpers import normalizeRNA, save_results
from models.xvae import xvae
from models.dec import XDEC
from models.metrics import calculate_metrics

parser = argparse.ArgumentParser(
    description='Train an X-DEC model (X-shaped VAE + Deep Embedded Clustering) on the AMD '
                 'cohort, integrating RNA-seq CPM expression (aak100_cpmdat.csv, numerical branch) '
                 'with clinical/genotype covariates (MetaSheet_1_4.csv, binary branch), matched by '
                 'sample_id. Port of the DAM-IC Deep-Embedded-Clustering-generalisability-and-'
                 'adaptation-for-mixed-data-types repo (scripts/xvae.py, DECLayer.py, '
                 'clustering_functions.py), rewritten in PyTorch to avoid TensorFlow/Keras '
                 'DLL-loading issues on Windows.'
)
parser.add_argument('--ds1', help='First hidden layer size for the numerical (RNA) branch',
                     type=int, default=48)
parser.add_argument('--ds2', help='First hidden layer size for the binary (clinical) branch. '
                                   'Defaults to --ds1 if omitted.', type=int, default=None)
parser.add_argument('--ds12', help='Hidden layer size that joins both branches. Defaults to --ds1.',
                     type=int, default=None)
parser.add_argument('--ls', help='Latent dimension size', type=int, default=16)
parser.add_argument('--act', help='Activation function for hidden layers', type=str, default='elu')
parser.add_argument('--dropout', help='Decoder dropout rate', type=float, default=0.2)
parser.add_argument('--distance', help='VAE regularizer: kl or mmd', type=str, default='kl')
parser.add_argument('--beta', help='Beta-VAE regularization weight', type=float, default=1.0)
parser.add_argument('--unweighted', help='Disable weighting reconstruction loss by branch size',
                     action='store_true')
parser.add_argument('--epochs', help='X-VAE pretraining epochs', type=int, default=150)
parser.add_argument('--bs', help='X-VAE pretraining batch size', type=int, default=32)
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
parser.add_argument('--save_model', help='Save the trained X-VAE encoder (for the SHAP workflow)',
                     action='store_true')
parser.add_argument('--writedir', help='Output directory - default is ./results', type=str, default='')
parser.add_argument('--seed', help='Random seed', type=int, default=5192)


def run_xdec(x_num, x_bin, args, y=None):
    """x_num: numerical RNA-seq matrix (xvae SET 1). x_bin: binary/one-hot
    clinical matrix (xvae SET 2). Pretrains the X-VAE, then runs X-DEC
    self-training clustering on top of its latent space."""
    args.s1_input_size = x_num.shape[1]
    args.s2_input_size = x_bin.shape[1]

    if args.s2_input_size > x_num.shape[0]:
        print('WARNING: clinical branch has {} columns but only {} samples. Consider '
              '--clinical_mode curated to reduce dimensionality.'.format(
                  args.s2_input_size, x_num.shape[0]))

    autoencoder = xvae(
        s1_input_size=args.s1_input_size, s2_input_size=args.s2_input_size,
        ds1=args.ds1, ds2=args.ds2, ds12=args.ds12, ls=args.ls, act=args.act,
        dropout=args.dropout, distance=args.distance, beta=args.beta,
        weighted=not args.unweighted, epochs=args.epochs, bs=args.bs,
    )
    autoencoder.build_model(seed=args.seed)
    autoencoder.train(x_num, x_bin, seed=args.seed)

    xdec = XDEC(autoencoder, n_clusters=args.n_clusters, seed=args.seed)
    y_pred, y_proba, centroids, z = xdec.fit(
        x_num, x_bin, y=y, maxiter=args.dec_maxiter, batch_size=args.dec_batch_size,
        update_interval=args.dec_update_interval, tol=args.dec_tol, n_init=args.n_init,
    )
    return autoencoder, xdec, y_pred, y_proba, centroids, z


if __name__ == '__main__':
    args = parser.parse_args()

    data = get_data(args.rna_file, args.clin_file, args.label_col, clinical_mode=args.clinical_mode)

    x_num = normalizeRNA(data['rnanp'])  # SET 1: numerical RNA-seq CPM
    x_bin = data['clin']                 # SET 2: binary clinical/genotype dummies
    print('X-DEC input: {} samples, {} genes (numeric) + {} clinical (binary)'.format(
        x_num.shape[0], x_num.shape[1], x_bin.shape[1]))

    root_dir = args.writedir or 'results'
    run_dir = os.path.join(
        root_dir, 'XDEC_Clin+mRNA_integration',
        'xdec_LS_{}_DS1_{}_DS2_{}_DS12_{}_{}_beta_{}_k{}'.format(
            args.ls, args.ds1, args.ds2 or args.ds1, args.ds12 or args.ds1,
            args.distance, args.beta, args.n_clusters),
    )
    os.makedirs(run_dir, exist_ok=True)

    autoencoder, xdec, y_pred, y_proba, centroids, z = run_xdec(x_num, x_bin, args, y=data['y'])

    if data['y'] is not None:
        acc, ari, nmi = calculate_metrics(data['y'], y_pred)
        print('Final clustering vs. {}: acc={:.5f}, ari={:.5f}, nmi={:.5f}'.format(
            args.label_col, acc, ari, nmi))

    if args.save_model:
        encoder_path = os.path.join(run_dir, 'encoder_xvae.pt')
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
    emb_path = os.path.join(run_dir, 'xdec_latent_embedding.csv')
    emb_df.to_csv(emb_path, index=False)
    print('Saved latent embedding to', emb_path, emb_df.shape)

    # ----- optional cluster-stability check -----
    if args.stability:
        from models.stability import compute_cluster_stability

        print('Running cluster-stability check (k={}, rep={}, {} total resamples)...'.format(
            args.stability_k, args.stability_rep, args.stability_k * args.stability_rep))
        xvae_kwargs = dict(ds1=args.ds1, ds2=args.ds2, ds12=args.ds12, ls=args.ls, act=args.act,
                           dropout=args.dropout, distance=args.distance, beta=args.beta,
                           weighted=not args.unweighted, epochs=args.epochs, bs=args.bs)
        dec_kwargs = dict(maxiter=args.dec_maxiter, batch_size=args.dec_batch_size,
                          update_interval=args.dec_update_interval, tol=args.dec_tol,
                          n_init=args.n_init)
        stability, y_pred_mat = compute_cluster_stability(
            x_num, x_bin, args.n_clusters, reference_labels=y_pred,
            xvae_kwargs=xvae_kwargs, dec_kwargs=dec_kwargs,
            k=args.stability_k, rep=args.stability_rep, seed=args.seed,
        )
        stability_path = os.path.join(run_dir, 'cluster_stability.csv')
        pd.DataFrame({'sample_id': data['sample_id'], 'cluster': y_pred,
                      'stability': stability.to_numpy()}).to_csv(stability_path, index=False)
        print('Mean cluster stability: {:.3f}'.format(stability.mean()))
        print('Saved', stability_path)
