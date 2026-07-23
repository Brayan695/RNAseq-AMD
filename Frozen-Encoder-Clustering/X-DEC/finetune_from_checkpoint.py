# Warm-starts a saved X-DEC encoder checkpoint into a full xvae wrapper
# (build_model() + a fresh Adam optimizer), loads the checkpoint's weights
# into it, then calls .train() again - e.g. to continue X-VAE pretraining
# under a changed loss function (models/xvae.py's xvae._loss()) before DEC
# clustering is (re-)run on top of it.
#
# Unlike load_encoder.py (inference only, via the existing load_xvae_model()
# helper - no optimizer, no .train()), this rebuilds the full wrapper with your
# own choice of epochs/bs/beta/distance/weighted, since load_xvae_model() only
# takes architecture args. The optimizer itself is NOT restored from the
# checkpoint (save_encoder() never saved optimizer.state_dict()) - this is a
# warm start from the saved weights, not a byte-for-byte resume of the
# original run.
#
# Reminder on xvae's loss (models/xvae.py's _loss()): the numeric branch (s1,
# RNA) uses MSE; the binary branch (s2, clinical) already uses
# binary_cross_entropy_with_logits, unlike CNC-VAE/CNC-DEC which use one MSE
# loss over everything concatenated. If you're changing s1's reconstruction
# loss to something cross-entropy-based too, its inputs need to be in [0, 1]
# first - see the warning this script prints below.
#
# This only re-trains the xvae encoder (the pretraining stage) - it does not
# re-run DEC self-training. To re-cluster on top of the fine-tuned encoder
# afterward, pass the saved output of this script's --save_model into
# XDEC(model, n_clusters=...).fit(...) yourself (see run_xdec.py for the
# exact call shape).
import argparse
import os

import numpy as np
import pandas as pd
import torch

from misc.dataset import DEFAULT_CLIN_FILE, DEFAULT_RNA_FILE, get_data
from misc.helpers import normalizeRNA
from models.xvae import xvae

parser = argparse.ArgumentParser(
    description='Warm-start a saved X-DEC encoder checkpoint and continue training it '
                "(e.g. after editing xvae._loss() to use a different loss function)."
)
parser.add_argument('--encoder', required=True, help='Path to the saved encoder .pt file to warm-start from')
parser.add_argument('--clinical_mode', type=str, default='curated', choices=['full', 'curated'],
                     help='Must match what the checkpoint was originally trained with')
parser.add_argument('--rna_file', type=str, default=DEFAULT_RNA_FILE)
parser.add_argument('--clin_file', type=str, default=DEFAULT_CLIN_FILE)
parser.add_argument('--label_col', type=str, default='mgs_level')
parser.add_argument('--epochs', type=int, default=250, help='Additional epochs to fine-tune for')
parser.add_argument('--bs', type=int, default=64)
parser.add_argument('--beta', type=float, default=25)
parser.add_argument('--distance', type=str, default='kl',
                     help="Read by xvae._loss() if it still branches on kl/mmd - ignored if you've "
                          "already hardcoded it to KL-only")
parser.add_argument('--unweighted', action='store_true',
                     help='Disable weighting reconstruction loss by branch size (see xvae._loss())')
parser.add_argument('--seed', type=int, default=5192)
parser.add_argument('--save_model', action='store_true', help='Save the fine-tuned encoder afterward')
parser.add_argument('--out', type=str, default='finetuned',
                     help='Output directory for the fine-tuned encoder + latent embedding CSV')


def main():
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    ckpt = torch.load(args.encoder, map_location='cpu')
    print('Warm-starting from: ds1={} ds2={} ds12={} ls={} dropout={} act={} '
          's1_input_size={} s2_input_size={}'.format(
              ckpt['ds1'], ckpt['ds2'], ckpt['ds12'], ckpt['ls'], ckpt['dropout'], ckpt['act'],
              ckpt['s1_input_size'], ckpt['s2_input_size']))

    # Architecture comes from the checkpoint (must match exactly, or load_state_dict
    # below will fail with a shape mismatch); training hyperparameters (epochs/bs/
    # beta/distance/weighted) are THIS run's own - they were never part of the checkpoint.
    model = xvae(
        s1_input_size=ckpt['s1_input_size'], s2_input_size=ckpt['s2_input_size'],
        ds1=ckpt['ds1'], ds2=ckpt['ds2'], ds12=ckpt['ds12'], ls=ckpt['ls'],
        act=ckpt['act'], dropout=ckpt['dropout'],
        distance=args.distance, beta=args.beta, weighted=not args.unweighted,
        epochs=args.epochs, bs=args.bs,
    )
    model.build_model(seed=args.seed)  # fresh network + fresh Adam optimizer, matching the checkpoint's shape
    model.net.load_state_dict(ckpt['state_dict'])  # warm start: overwrite the fresh init with saved weights

    d = get_data(args.rna_file, args.clin_file, args.label_col, clinical_mode=args.clinical_mode)
    x_num = normalizeRNA(d['rnanp'])   # SET 1: numerical RNA-seq CPM
    x_bin = d['clin']                  # SET 2: binary clinical/genotype dummies

    if x_num.shape[1] != ckpt['s1_input_size'] or x_bin.shape[1] != ckpt['s2_input_size']:
        raise ValueError(
            'Input widths (numeric={}, binary={}) do not match checkpoint '
            '(s1_input_size={}, s2_input_size={}) - check --clinical_mode matches '
            'the original training run'.format(
                x_num.shape[1], x_bin.shape[1], ckpt['s1_input_size'], ckpt['s2_input_size'])
        )
    if x_num.min() < 0.0 or x_num.max() > 1.0:
        print('WARNING: numeric (RNA) branch contains values outside [0, 1] (min={:.3f}, max={:.3f}) - '
              "if xvae._loss()'s s1 reconstruction now uses cross-entropy-style loss, values outside "
              '[0, 1] will produce garbage/NaN loss.'.format(x_num.min(), x_num.max()))

    print('Fine-tuning for {} epochs...'.format(args.epochs))
    model.train(x_num, x_bin, seed=args.seed)

    z = model.predict(x_num, x_bin, output='encoder')
    print('latent embedding:', z.shape)

    if args.save_model:
        encoder_path = os.path.join(args.out, 'encoder_xvae_finetuned.pt')
        model.save_encoder(encoder_path)
        print('Saved fine-tuned encoder to', encoder_path)

    out_df = pd.DataFrame({'sample_id': d['sample_id'], args.label_col: d['label_classes'][d['y']]})
    for j in range(z.shape[1]):
        out_df['z{}'.format(j)] = z[:, j]
    out_path = os.path.join(args.out, 'finetuned_latent_embedding.csv')
    out_df.to_csv(out_path, index=False)
    print('Saved', out_path, out_df.shape)


if __name__ == '__main__':
    main()
