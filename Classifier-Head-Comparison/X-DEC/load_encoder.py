# Reloads a standalone encoder saved by xvae.save_encoder() (encoder_xvae*.pt,
# e.g. via run_xdec.py --save_model or XDEC Runner.ipynb) and re-runs it over
# the AMD cohort to get the latent embedding back out, without retraining.
#
# Unlike CNC-VAE's load_encoder.py, this folder's models/xvae.py already ships
# a load_xvae_model(path) helper that rebuilds the full xvae wrapper (weights +
# dual-branch architecture) in one call - this script is a thin CLI around
# that. Note the dual-input shape: xvae keeps clinical (s2, binary) and RNA
# (s1, numeric) as two SEPARATE arrays all the way through - never concatenated
# like CNC-VAE/CNC-DEC - so predict() takes two arguments, not one.
#
# X-DEC only ever normalizes RNA with plain min-max (misc.helpers.normalizeRNA) -
# there's no log1p-normalized variant script in this folder, so there's no
# --normalization flag to get wrong here.
import argparse
import os

import pandas as pd

from misc.dataset import DEFAULT_CLIN_FILE, DEFAULT_RNA_FILE, get_data
from misc.helpers import normalizeRNA
from models.xvae import load_xvae_model

parser = argparse.ArgumentParser(
    description='Reload a saved X-DEC encoder and encode the AMD cohort into its latent space.'
)
parser.add_argument('--encoder', required=True, help='Path to the saved encoder .pt file')
parser.add_argument('--clinical_mode', type=str, default='curated', choices=['full', 'curated'],
                     help="Must match what the checkpoint was trained with (run_xdec.py's own "
                          "default is 'curated')")
parser.add_argument('--rna_file', type=str, default=DEFAULT_RNA_FILE)
parser.add_argument('--clin_file', type=str, default=DEFAULT_CLIN_FILE)
parser.add_argument('--label_col', type=str, default='mgs_level')
parser.add_argument('--out', type=str, default='', help='Output CSV path for the latent embedding '
                                                          '(default: alongside --encoder)')


def main():
    args = parser.parse_args()

    model = load_xvae_model(args.encoder)
    a = model.args
    print('Loaded encoder: ds1={} ds2={} ds12={} ls={} dropout={} act={} '
          's1_input_size={} s2_input_size={}'.format(
              a.ds1, a.ds2, a.ds12, a.ls, a.dropout, a.act, a.s1_input_size, a.s2_input_size))

    d = get_data(args.rna_file, args.clin_file, args.label_col, clinical_mode=args.clinical_mode)
    x_num = normalizeRNA(d['rnanp'])   # SET 1: numerical RNA-seq CPM
    x_bin = d['clin']                  # SET 2: binary clinical/genotype dummies

    if x_num.shape[1] != a.s1_input_size or x_bin.shape[1] != a.s2_input_size:
        raise ValueError(
            'Input widths (numeric={}, binary={}) do not match checkpoint '
            '(s1_input_size={}, s2_input_size={}) - check --clinical_mode matches '
            'how this encoder was trained'.format(
                x_num.shape[1], x_bin.shape[1], a.s1_input_size, a.s2_input_size)
        )

    z = model.predict(x_num, x_bin, output='encoder')
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
