# Loads a previously-trained CNC-VAE encoder (no retraining - see
# models/dec.py:DEC, the frozen-encoder variant in this copy of CNC-DEC) and
# clusters on top of its FROZEN latent embedding, then reports AMI and ARI
# against the known mgs_level label.
#
# Companion to ../DEC/compare_clusters.py and ../X-DEC/compare_clusters.py -
# each writes its own JSON result to ../results/, then run
# ../summarize_comparison.py to see all three side by side. Run this script
# from INSIDE this folder (it imports models/misc as local packages, same as
# every other script in the original CNC-DEC project).
import argparse
import json
import os

from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

from misc.dataset import DEFAULT_CLIN_FILE, DEFAULT_RNA_FILE, get_data
from misc.helpers import normalizeRNA
from models.cncvae import load_cncvae_model
from models.dec import DEC

parser = argparse.ArgumentParser(
    description='Cluster a previously-trained (frozen, not retrained) CNC-VAE encoder '
                'and report AMI/ARI against the known mgs_level label.'
)
parser.add_argument('--encoder', required=True, help='Path to the saved encoder .pt file')
parser.add_argument('--clinical_mode', type=str, default='curated', choices=['full', 'curated'],
                     help='Must match what the checkpoint was originally trained with')
parser.add_argument('--rna_file', type=str, default=DEFAULT_RNA_FILE)
parser.add_argument('--clin_file', type=str, default=DEFAULT_CLIN_FILE)
parser.add_argument('--label_col', type=str, default='mgs_level')
parser.add_argument('--n_clusters', type=int, default=2)
parser.add_argument('--seed', type=int, default=5192)
parser.add_argument('--out', type=str, default='../results/cncdec_frozen_cluster.json')


def main():
    args = parser.parse_args()

    ae = load_cncvae_model(args.encoder)  # build_model() + load_state_dict() + eval() already done
    a = ae.args
    print('Loaded encoder: ds={} ls={} dropout={} act={} input_size={}'.format(
        a.ds, a.ls, a.dropout, a.act, a.input_size))

    d = get_data(args.rna_file, args.clin_file, args.label_col, clinical_mode=args.clinical_mode)
    x_num = normalizeRNA(d['rnanp'])
    x_bin = d['clin']

    if x_num.shape[1] + x_bin.shape[1] != a.input_size:
        raise ValueError(
            'Input width {} does not match checkpoint input_size {} - '
            'check --clinical_mode matches how this encoder was trained'.format(
                x_num.shape[1] + x_bin.shape[1], a.input_size)
        )

    # DEC here is the frozen-encoder variant (see models/dec.py) - it still
    # runs the real K-means-init + Q/P self-training algorithm, it just never
    # updates the encoder's weights, only the DECHead's cluster centers.
    dec = DEC(ae, n_clusters=args.n_clusters, seed=args.seed)
    y_pred, y_proba, centroids, z = dec.fit(x_num, x_bin, y=d['y'])

    ami = adjusted_mutual_info_score(d['y'], y_pred)
    ari = adjusted_rand_score(d['y'], y_pred)
    print('AMI={:.5f}  ARI={:.5f}'.format(ami, ari))

    os.makedirs(os.path.dirname(args.out) or '.', exist_ok=True)
    with open(args.out, 'w') as f:
        json.dump({
            'model': 'CNC-DEC', 'encoder': args.encoder, 'ls': a.ls, 'n_clusters': args.n_clusters,
            'ami': ami, 'ari': ari, 'n_samples': int(len(d['y'])),
        }, f, indent=2)
    print('Saved', args.out)


if __name__ == '__main__':
    main()
