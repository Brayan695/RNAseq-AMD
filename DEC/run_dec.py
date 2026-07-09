import argparse
import os

import numpy as np

from misc.dataset import get_data, DEFAULT_CLIN_FILE, DEFAULT_RNA_FILE
from misc.helpers import normalizeRNA, save_results
from models.autoencoder import MLPAutoencoder
from models.dec import DEC
from models.metrics import calculate_metrics

parser = argparse.ArgumentParser(
    description='Train a plain Deep Embedded Clustering (DEC) model on the AMD cohort. Port of '
                 'the DAM-IC Deep-Embedded-Clustering-generalisability-and-adaptation-for-mixed-'
                 'data-types repo (scripts/DECLayer.py + clustering_functions.cluster_mlp_'
                 'autoencoder/dec_cluster), rewritten in PyTorch to avoid TensorFlow/Keras '
                 'DLL-loading issues on Windows. This is the plain single-input DEC (one MLP '
                 'autoencoder), not the X-shaped/mixed-type variant - RNA-seq CPM expression '
                 '(aak100_cpmdat.csv) and binary clinical/genotype dummies (MetaSheet_1_4.csv) '
                 'are matched by sample_id and concatenated into one feature matrix.'
)
parser.add_argument('--neurons_h', help='Hidden layer size', type=int, default=64)
parser.add_argument('--neurons_e', help='Embedding (latent) layer size', type=int, default=8)
parser.add_argument('--epochs', help='Autoencoder pretraining epochs', type=int, default=500)
parser.add_argument('--bs', help='Autoencoder pretraining batch size', type=int, default=64)
parser.add_argument('--n_clusters', help='Number of clusters to identify', type=int, default=2)
parser.add_argument('--dec_maxiter', help='Max DEC self-training iterations', type=int, default=2000)
parser.add_argument('--dec_update_interval', help='Iterations between target-distribution updates',
                     type=int, default=50)
parser.add_argument('--dec_batch_size', help='DEC self-training batch size', type=int, default=64)
parser.add_argument('--dec_tol', help='Stop when the fraction of samples changing cluster '
                                       'assignment drops below this', type=float, default=0.01)
parser.add_argument('--n_init', help='Number of K-means initializations for cluster centers',
                     type=int, default=20)
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
parser.add_argument('--save_model', help='Save the trained autoencoder + DEC cluster centers',
                     action='store_true')
parser.add_argument('--writedir', help='Output directory - default is ./results', type=str, default='')
parser.add_argument('--seed', help='Random seed', type=int, default=5192)

if __name__ == '__main__':
    args = parser.parse_args()

    data = get_data(args.rna_file, args.clin_file, args.label_col, clinical_mode=args.clinical_mode)

    # Single combined feature matrix: RNA-seq min-max scaled to [0, 1] (matches
    # the sigmoid decoder output), clinical dummies already 0/1.
    X = np.hstack([normalizeRNA(data['rnanp']), data['clin']]).astype(np.float32)
    print('Combined feature matrix: {} samples x {} features ({} genes + {} clinical)'.format(
        X.shape[0], X.shape[1], data['rnanp'].shape[1], data['clin'].shape[1]))

    root_dir = args.writedir or 'results'
    run_dir = os.path.join(
        root_dir, 'DEC',
        'dec_H{}_E{}_k{}'.format(args.neurons_h, args.neurons_e, args.n_clusters),
    )
    os.makedirs(run_dir, exist_ok=True)

    autoencoder = MLPAutoencoder(
        n_features=X.shape[1], neurons_h=args.neurons_h, neurons_e=args.neurons_e,
        epochs=args.epochs, bs=args.bs,
    )
    autoencoder.build_model(seed=args.seed)
    autoencoder.train(X, seed=args.seed)

    dec = DEC(autoencoder, n_clusters=args.n_clusters, seed=args.seed)
    y_pred, y_proba, centroids, z = dec.fit(
        X, y=data['y'], maxiter=args.dec_maxiter, batch_size=args.dec_batch_size,
        update_interval=args.dec_update_interval, tol=args.dec_tol, n_init=args.n_init,
    )

    if data['y'] is not None:
        acc, ari, nmi = calculate_metrics(data['y'], y_pred)
        print('Final clustering vs. {}: acc={:.5f}, ari={:.5f}, nmi={:.5f}'.format(
            args.label_col, acc, ari, nmi))

    if args.save_model:
        encoder_path = os.path.join(run_dir, 'encoder_dec.pt')
        autoencoder.save_encoder(encoder_path)
        print('Saved encoder to', encoder_path)

    save_results(
        run_dir, args.label_col + '.npz',
        z=z, y_pred=y_pred, y_proba=y_proba, centroids=centroids,
        y=data['y'], sample_id=data['sample_id'], classes=data['label_classes'],
    )
    print('Saved results to', os.path.join(run_dir, args.label_col + '.npz'))
