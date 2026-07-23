# Reloads a standalone autoencoder saved by MLPAutoencoder.save_encoder()
# (e.g. via run_dec.py --save_model) and re-runs it over the AMD cohort to get
# the latent embedding back out, without retraining.
#
# Unlike the VAE-based variants (CNC-VAE/CNC-DEC/X-DEC), this autoencoder's
# checkpoint is TWO files: the weights (path) plus a sidecar architecture CSV
# (path + '.architecture.csv') - models/autoencoder.py's load_autoencoder(path)
# already knows how to read both together into a ready-to-use MLPAutoencoder,
# so this script is a thin CLI around that, matching the other three folders'
# load_encoder.py in shape/output.
#
# Also unlike the VAE variants, MLPAutoencoder.predict() takes ONE combined
# array (RNA + clinical already concatenated), not two separate arguments -
# see build_combined() below, which matches run_dec.py's own preprocessing
# exactly (plain min-max RNA normalization, no log1p variant in this folder).
import argparse
import os

import numpy as np
import pandas as pd

from misc.dataset import DEFAULT_CLIN_FILE, DEFAULT_RNA_FILE, get_data
from misc.helpers import normalizeRNA
from models.autoencoder import load_autoencoder

parser = argparse.ArgumentParser(
    description='Reload a saved DEC autoencoder and encode the AMD cohort into its latent space.'
)
parser.add_argument('--encoder', required=True,
                     help='Path to the saved encoder .pt file (its .architecture.csv sidecar '
                          'must sit right next to it)')
parser.add_argument('--clinical_mode', type=str, default='curated', choices=['full', 'curated'],
                     help="Must match what the checkpoint was trained with (run_dec.py's own "
                          "default is 'curated')")
parser.add_argument('--rna_file', type=str, default=DEFAULT_RNA_FILE)
parser.add_argument('--clin_file', type=str, default=DEFAULT_CLIN_FILE)
parser.add_argument('--label_col', type=str, default='mgs_level')
parser.add_argument('--out', type=str, default='', help='Output CSV path for the latent embedding '
                                                          '(default: alongside --encoder)')


def build_combined(d):
    """Matches run_dec.py's own preprocessing exactly: RNA min-max scaled to
    [0, 1], concatenated with the (already 0/1) clinical dummies."""
    return np.hstack([normalizeRNA(d['rnanp']), d['clin']]).astype(np.float32)


def main():
    args = parser.parse_args()

    model = load_autoencoder(args.encoder)
    a = model.args
    print('Loaded encoder: neurons_h={} neurons_e={} l1_coef={} n_features={}'.format(
        a.neurons_h, a.neurons_e, a.l1_coef, a.n_features))

    d = get_data(args.rna_file, args.clin_file, args.label_col, clinical_mode=args.clinical_mode)
    X = build_combined(d)

    if X.shape[1] != a.n_features:
        raise ValueError(
            'Input width {} does not match checkpoint n_features {} - '
            'check --clinical_mode matches how this encoder was trained'.format(
                X.shape[1], a.n_features)
        )

    z = model.predict(X)
    print('latent embedding:', z.shape)

    out_path = args.out or os.path.join(
        os.path.dirname(args.encoder) or '.', 'reencoded_latent_embedding.csv'
    )
    out_df = pd.DataFrame({'sample_id': d['sample_id'], args.label_col: d['label_classes'][d['y']]})
    for j in range(z.shape[1]):
        out_df['z{}'.format(j)] = z[:, j]
    out_df.to_csv(out_path, index=False)
    print('Saved', out_path, out_df.shape)


if __name__ == '__main__':
    main()
