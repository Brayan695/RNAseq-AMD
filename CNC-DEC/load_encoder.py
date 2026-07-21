# Reloads a standalone encoder saved by CNCVAE.save_encoder() (encoder_cncvae*.pt,
# e.g. via run_cncdec.py --save_model or CNC-DEC Runner.ipynb) and re-runs it over
# the AMD cohort to get the latent embedding back out, without retraining.
#
# Unlike CNC-VAE's load_encoder.py, this folder's models/cncvae.py already ships
# a load_cncvae_model(path) helper that rebuilds the full CNCVAE wrapper (weights
# + architecture) in one call - this script is a thin CLI around that, matching
# CNC-VAE's load_encoder.py in shape/output so the two are easy to compare.
#
# CNC-DEC only ever normalizes RNA with plain min-max (misc.helpers.normalizeRNA) -
# unlike CNC-VAE, there's no log1p-normalized variant script in this folder, so
# there's no --normalization flag to get wrong here.
import argparse
import os

import pandas as pd

from misc.dataset import DEFAULT_CLIN_FILE, DEFAULT_RNA_FILE, get_data
from misc.helpers import normalizeRNA
from models.cncvae import load_cncvae_model

parser = argparse.ArgumentParser(
    description='Reload a saved CNC-DEC encoder and encode the AMD cohort into its latent space.'
)
parser.add_argument('--encoder', required=True, help='Path to the saved encoder .pt file')
parser.add_argument('--clinical_mode', type=str, default='curated', choices=['full', 'curated'],
                     help="Must match what the checkpoint was trained with (run_cncdec.py's own "
                          "default is 'curated')")
parser.add_argument('--rna_file', type=str, default=DEFAULT_RNA_FILE)
parser.add_argument('--clin_file', type=str, default=DEFAULT_CLIN_FILE)
parser.add_argument('--label_col', type=str, default='mgs_level')
parser.add_argument('--out', type=str, default='', help='Output CSV path for the latent embedding '
                                                          '(default: alongside --encoder)')


def main():
    args = parser.parse_args()

    cncvae = load_cncvae_model(args.encoder)
    print('Loaded encoder: ds={} ls={} dropout={} act={} input_size={}'.format(
        cncvae.args.ds, cncvae.args.ls, cncvae.args.dropout, cncvae.args.act, cncvae.args.input_size))

    d = get_data(args.rna_file, args.clin_file, args.label_col, clinical_mode=args.clinical_mode)
    X_clin = d['clin']
    X_rna = normalizeRNA(d['rnanp'])

    input_size = X_clin.shape[1] + X_rna.shape[1]
    if input_size != cncvae.args.input_size:
        raise ValueError(
            'Input width {} does not match checkpoint input_size {} - '
            'check --clinical_mode matches how this encoder was trained'.format(
                input_size, cncvae.args.input_size)
        )

    z = cncvae.predict(X_clin, X_rna)
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
