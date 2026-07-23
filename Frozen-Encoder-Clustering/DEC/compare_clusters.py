# Loads a previously-trained DEC (plain MLP) autoencoder (no retraining - see
# models/dec.py:DEC, the frozen-encoder variant in this copy of DEC) and
# clusters on top of its FROZEN latent embedding, then reports AMI and ARI
# against the known mgs_level label.
#
# Companion to ../CNC-DEC/compare_clusters.py and ../X-DEC/compare_clusters.py -
# each writes its own JSON result to ../results/, then run
# ../summarize_comparison.py to see all three side by side. Run this script
# from INSIDE this folder (it imports models/misc as local packages, same as
# every other script in the original DEC project).
import argparse
import json
import os

import numpy as np
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

from misc.dataset import DEFAULT_CLIN_FILE, DEFAULT_RNA_FILE, get_data
from misc.helpers import normalizeRNA
from models.autoencoder import load_autoencoder
from models.dec import DEC

parser = argparse.ArgumentParser(
    description='Cluster a previously-trained (frozen, not retrained) DEC autoencoder '
                'and report AMI/ARI against the known mgs_level label.'
)
parser.add_argument('--encoder', required=True,
                     help='Path to the saved encoder .pt file (its .architecture.csv sidecar '
                          'must sit right next to it)')
parser.add_argument('--clinical_mode', type=str, default='curated', choices=['full', 'curated'],
                     help='Must match what the checkpoint was originally trained with')
parser.add_argument('--rna_file', type=str, default=DEFAULT_RNA_FILE)
parser.add_argument('--clin_file', type=str, default=DEFAULT_CLIN_FILE)
parser.add_argument('--label_col', type=str, default='mgs_level')
parser.add_argument('--n_clusters', type=int, default=2)
parser.add_argument('--seed', type=int, default=5192)
parser.add_argument('--out', type=str, default='../results/dec_frozen_cluster.json')


def main():
    args = parser.parse_args()

    ae = load_autoencoder(args.encoder)  # build_model() + load_state_dict() + eval() already done
    a = ae.args
    print('Loaded encoder: neurons_h={} neurons_e={} l1_coef={} n_features={}'.format(
        a.neurons_h, a.neurons_e, a.l1_coef, a.n_features))

    d = get_data(args.rna_file, args.clin_file, args.label_col, clinical_mode=args.clinical_mode)
    X = np.hstack([normalizeRNA(d['rnanp']), d['clin']]).astype(np.float32)

    if X.shape[1] != a.n_features:
        raise ValueError(
            'Input width {} does not match checkpoint n_features {} - '
            'check --clinical_mode matches how this encoder was trained'.format(
                X.shape[1], a.n_features)
        )

    # DEC here is the frozen-encoder variant (see models/dec.py) - it still
    # runs the real K-means-init + Q/P self-training algorithm, it just never
    # updates the encoder's weights, only the DECHead's cluster centers.
    dec = DEC(ae, n_clusters=args.n_clusters, seed=args.seed)
    y_pred, y_proba, centroids, z = dec.fit(X, y=d['y'])

    ami = adjusted_mutual_info_score(d['y'], y_pred)
    ari = adjusted_rand_score(d['y'], y_pred)
    print('AMI={:.5f}  ARI={:.5f}'.format(ami, ari))

    os.makedirs(os.path.dirname(args.out) or '.', exist_ok=True)
    with open(args.out, 'w') as f:
        json.dump({
            'model': 'DEC', 'encoder': args.encoder, 'ls': a.neurons_e, 'n_clusters': args.n_clusters,
            'ami': ami, 'ari': ari, 'n_samples': int(len(d['y'])),
        }, f, indent=2)
    print('Saved', args.out)


if __name__ == '__main__':
    main()
