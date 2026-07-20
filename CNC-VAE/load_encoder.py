# Reload a standalone encoder saved by CNCVAE.save_encoder() (encoder_cncvae*.pt)
# and re-run it over the AMD cohort to get the latent embedding back out, without
# retraining. Unlike run_cncvae.py/compare_representations.py, this script never
# calls build_model()/train() - it only reconstructs the exact architecture the
# checkpoint was saved with and loads its weights.
#
# The checkpoint's own ds/ls/dropout/act/input_size are used (not passed on the
# command line), since changing any of those would mean the saved weights no
# longer fit the network shape. What you DO need to get right is how the input
# was preprocessed at training time - see --normalization below.
import argparse
import os

import numpy as np
import pandas as pd
import torch

from misc.dataset import DEFAULT_CLIN_FILE, DEFAULT_RNA_FILE, get_data
from misc.helpers import normalizeRNA
from models.cncvae import _CNCVAENet

parser = argparse.ArgumentParser(
    description='Reload a saved CNC-VAE encoder and encode the AMD cohort into its latent space.'
)
parser.add_argument('--encoder', required=True, help='Path to the saved encoder .pt file')
parser.add_argument('--clinical_mode', type=str, default='curated', choices=['full', 'curated'],
                     help="Must match what the checkpoint was trained with: 'curated' for "
                          "compare_representations.py, 'full' or 'curated' for run_cncvae.py "
                          "depending on how it was invoked")
parser.add_argument('--normalization', type=str, default='minmax', choices=['minmax', 'log1p_minmax'],
                     help="RNA preprocessing to match the checkpoint's training script: "
                          "'minmax' = run_cncvae.py (misc.helpers.normalizeRNA), "
                          "'log1p_minmax' = compare_representations.py (log1p then min-max)")
parser.add_argument('--rna_file', type=str, default=DEFAULT_RNA_FILE)
parser.add_argument('--clin_file', type=str, default=DEFAULT_CLIN_FILE)
parser.add_argument('--label_col', type=str, default='mgs_level')
parser.add_argument('--out', type=str, default='', help='Output CSV path for the latent embedding '
                                                          '(default: alongside --encoder)')


def build_input(args, d):
    """Reproduce the exact (clinical, RNA) preprocessing + concatenation order
    ([clin, mrna]) that both training scripts use, so the checkpoint sees input
    on the same scale it was trained on."""
    X_clin = d['clin'].astype(np.float32)

    if args.normalization == 'log1p_minmax':
        mrna_log = np.log1p(d['rnanp'].astype(np.float64))
        mrna_min, mrna_max = mrna_log.min(axis=0), mrna_log.max(axis=0)
        mrna_range = np.where(mrna_max - mrna_min == 0, 1.0, mrna_max - mrna_min)
        X_mrna = ((mrna_log - mrna_min) / mrna_range).astype(np.float32)
    else:
        X_mrna = normalizeRNA(d['rnanp'])

    return np.concatenate([X_clin, X_mrna], axis=1).astype(np.float32)


def main():
    args = parser.parse_args()

    ckpt = torch.load(args.encoder, map_location='cpu')
    net = _CNCVAENet(ckpt['input_size'], ckpt['ds'], ckpt['ls'], ckpt['dropout'], ckpt['act'])
    net.load_state_dict(ckpt['state_dict'])
    net.eval()
    print('Loaded encoder: ds={} ls={} dropout={} act={} input_size={}'.format(
        ckpt['ds'], ckpt['ls'], ckpt['dropout'], ckpt['act'], ckpt['input_size']))

    d = get_data(args.rna_file, args.clin_file, args.label_col, clinical_mode=args.clinical_mode)
    X = build_input(args, d)

    if X.shape[1] != ckpt['input_size']:
        # A mismatch here almost always means --clinical_mode doesn't match what
        # this checkpoint was trained with - fail loudly rather than silently
        # feeding the network a misaligned input vector.
        raise ValueError(
            'Input width {} does not match checkpoint input_size {} - '
            'check --clinical_mode matches how this encoder was trained'.format(
                X.shape[1], ckpt['input_size'])
        )

    with torch.no_grad():
        z_mean, _ = net.encode(torch.tensor(X, dtype=torch.float32))
    Z = z_mean.numpy()
    print('latent embedding:', Z.shape)

    out_path = args.out or os.path.join(
        os.path.dirname(args.encoder) or '.', 'reencoded_latent_embedding.csv'
    )
    out_df = pd.DataFrame({'sample_id': d['sample_id'], args.label_col: d['label_classes'][d['y']]})
    for j in range(Z.shape[1]):
        out_df['z{}'.format(j)] = Z[:, j]
    out_df.to_csv(out_path, index=False)
    print('Saved', out_path, out_df.shape)


if __name__ == '__main__':
    main()
