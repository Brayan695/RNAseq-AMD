# Warm-starts a saved CNC-DEC encoder checkpoint into a full CNCVAE wrapper
# (build_model() + a fresh Adam optimizer), loads the checkpoint's weights
# into it, then calls .train() again - e.g. to continue VAE pretraining under
# a changed loss function (models/cncvae.py's CNCVAE._loss()) before DEC
# clustering is (re-)run on top of it.
#
# Unlike load_encoder.py (inference only, via the existing load_cncvae_model()
# helper - no optimizer, no .train()), this rebuilds the full wrapper with your
# own choice of epochs/bs/beta/distance, since load_cncvae_model() only takes
# architecture args. The optimizer itself is NOT restored from the checkpoint
# (save_encoder() never saved optimizer.state_dict()) - this is a warm start
# from the saved weights, not a byte-for-byte resume of the original run.
#
# This only re-trains the CNCVAE encoder (the pretraining stage) - it does not
# re-run DEC self-training. To re-cluster on top of the fine-tuned encoder
# afterward, pass the saved output of this script's --save_model into
# DEC(cncvae, n_clusters=...).fit(...) yourself (see run_cncdec.py for the
# exact call shape).
import argparse
import os

import numpy as np
import pandas as pd
import torch

from misc.dataset import DEFAULT_CLIN_FILE, DEFAULT_RNA_FILE, get_data
from misc.helpers import normalizeRNA
from models.cncvae import CNCVAE

parser = argparse.ArgumentParser(
    description='Warm-start a saved CNC-DEC encoder checkpoint and continue training it '
                '(e.g. after editing CNCVAE._loss() to use a different loss function).'
)
parser.add_argument('--encoder', required=True, help='Path to the saved encoder .pt file to warm-start from')
parser.add_argument('--clinical_mode', type=str, default='curated', choices=['full', 'curated'],
                     help='Must match what the checkpoint was originally trained with')
parser.add_argument('--rna_file', type=str, default=DEFAULT_RNA_FILE)
parser.add_argument('--clin_file', type=str, default=DEFAULT_CLIN_FILE)
parser.add_argument('--label_col', type=str, default='mgs_level')
parser.add_argument('--epochs', type=int, default=150, help='Additional epochs to fine-tune for')
parser.add_argument('--bs', type=int, default=32)
parser.add_argument('--beta', type=float, default=1.0)
parser.add_argument('--distance', type=str, default='kl',
                     help="Read by CNCVAE._loss() if it still branches on kl/mmd - ignored if you've "
                          "already hardcoded it to KL-only")
parser.add_argument('--seed', type=int, default=5192)
parser.add_argument('--save_model', action='store_true', help='Save the fine-tuned encoder afterward')
parser.add_argument('--out', type=str, default='finetuned',
                     help='Output directory for the fine-tuned encoder + latent embedding CSV')


def main():
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    ckpt = torch.load(args.encoder, map_location='cpu')
    print('Warm-starting from: ds={} ls={} dropout={} act={} input_size={}'.format(
        ckpt['ds'], ckpt['ls'], ckpt['dropout'], ckpt['act'], ckpt['input_size']))

    # Architecture comes from the checkpoint (must match exactly, or load_state_dict
    # below will fail with a shape mismatch); training hyperparameters (epochs/bs/
    # beta/distance) are THIS run's own - they were never part of the checkpoint.
    cncvae = CNCVAE(
        input_size=ckpt['input_size'], ds=ckpt['ds'], ls=ckpt['ls'],
        act=ckpt['act'], dropout=ckpt['dropout'],
        distance=args.distance, beta=args.beta, epochs=args.epochs, bs=args.bs,
    )
    cncvae.build_model(seed=args.seed)  # fresh network + fresh Adam optimizer, matching the checkpoint's shape
    cncvae.net.load_state_dict(ckpt['state_dict'])  # warm start: overwrite the fresh init with saved weights

    d = get_data(args.rna_file, args.clin_file, args.label_col, clinical_mode=args.clinical_mode)
    X_clin = d['clin']
    X_rna = normalizeRNA(d['rnanp'])

    input_size = X_clin.shape[1] + X_rna.shape[1]
    if input_size != ckpt['input_size']:
        raise ValueError(
            'Input width {} does not match checkpoint input_size {} - '
            'check --clinical_mode matches the original training run'.format(
                input_size, ckpt['input_size'])
        )
    combined_min = min(X_clin.min(), X_rna.min())
    combined_max = max(X_clin.max(), X_rna.max())
    if combined_min < 0.0 or combined_max > 1.0:
        print('WARNING: input contains values outside [0, 1] (min={:.3f}, max={:.3f}) - '
              'if CNCVAE._loss() now uses cross-entropy reconstruction, values outside '
              '[0, 1] will produce garbage/NaN loss.'.format(combined_min, combined_max))

    print('Fine-tuning for {} epochs...'.format(args.epochs))
    cncvae.train(X_clin, X_rna, seed=args.seed)

    z = cncvae.predict(X_clin, X_rna)
    print('latent embedding:', z.shape)

    if args.save_model:
        encoder_path = os.path.join(args.out, 'encoder_cncvae_finetuned.pt')
        cncvae.save_encoder(encoder_path)
        print('Saved fine-tuned encoder to', encoder_path)

    out_df = pd.DataFrame({'sample_id': d['sample_id'], args.label_col: d['label_classes'][d['y']]})
    for j in range(z.shape[1]):
        out_df['z{}'.format(j)] = z[:, j]
    out_path = os.path.join(args.out, 'finetuned_latent_embedding.csv')
    out_df.to_csv(out_path, index=False)
    print('Saved', out_path, out_df.shape)


if __name__ == '__main__':
    main()
