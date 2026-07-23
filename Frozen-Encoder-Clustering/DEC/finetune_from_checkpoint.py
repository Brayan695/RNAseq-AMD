# Warm-starts a saved DEC autoencoder checkpoint and continues training it -
# e.g. to keep training for more epochs, or after editing
# MLPAutoencoder._loss() in models/autoencoder.py to use a different loss.
#
# Unlike load_encoder.py (inference only - no further training), this reuses
# the same load_autoencoder(path) helper to get a fully warm-started model
# (weights loaded, optimizer built), then overrides its training
# hyperparameters (epochs/bs) for this run before calling .train() again.
# The optimizer itself is NOT restored from the checkpoint (save_encoder()
# never saved optimizer.state_dict()) - this is a warm start from the saved
# weights, not a byte-for-byte resume of the original training run.
import argparse
import os

import numpy as np
import pandas as pd

from misc.dataset import DEFAULT_CLIN_FILE, DEFAULT_RNA_FILE, get_data
from misc.helpers import normalizeRNA
from models.autoencoder import load_autoencoder

parser = argparse.ArgumentParser(
    description='Warm-start a saved DEC autoencoder checkpoint and continue training it.'
)
parser.add_argument('--encoder', required=True,
                     help='Path to the saved encoder .pt file (its .architecture.csv sidecar '
                          'must sit right next to it)')
parser.add_argument('--clinical_mode', type=str, default='curated', choices=['full', 'curated'],
                     help='Must match what the checkpoint was originally trained with')
parser.add_argument('--rna_file', type=str, default=DEFAULT_RNA_FILE)
parser.add_argument('--clin_file', type=str, default=DEFAULT_CLIN_FILE)
parser.add_argument('--label_col', type=str, default='mgs_level')
parser.add_argument('--epochs', type=int, default=500, help='Additional epochs to fine-tune for')
parser.add_argument('--bs', type=int, default=64)
parser.add_argument('--seed', type=int, default=5192)
parser.add_argument('--save_model', action='store_true', help='Save the fine-tuned encoder afterward')
parser.add_argument('--out', type=str, default='finetuned',
                     help='Output directory for the fine-tuned encoder + latent embedding CSV')


def build_combined(d):
    """Matches run_dec.py's own preprocessing exactly: RNA min-max scaled to
    [0, 1], concatenated with the (already 0/1) clinical dummies."""
    return np.hstack([normalizeRNA(d['rnanp']), d['clin']]).astype(np.float32)


def main():
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    model = load_autoencoder(args.encoder)  # build_model() + load_state_dict() + eval() already done
    a = model.args
    print('Warm-starting from: neurons_h={} neurons_e={} l1_coef={} n_features={}'.format(
        a.neurons_h, a.neurons_e, a.l1_coef, a.n_features))

    # Training hyperparameters (epochs/bs) are THIS run's own - they were never
    # part of the checkpoint, so override whatever load_autoencoder() left on
    # model.args (its own defaults) with what was actually asked for here.
    a.epochs = args.epochs
    a.bs = args.bs

    d = get_data(args.rna_file, args.clin_file, args.label_col, clinical_mode=args.clinical_mode)
    X = build_combined(d)

    if X.shape[1] != a.n_features:
        raise ValueError(
            'Input width {} does not match checkpoint n_features {} - '
            'check --clinical_mode matches the original training run'.format(
                X.shape[1], a.n_features)
        )
    if X.min() < 0.0 or X.max() > 1.0:
        print('WARNING: input contains values outside [0, 1] (min={:.3f}, max={:.3f}) - '
              "if MLPAutoencoder's decoder/loss changes assume a [0, 1] target, values "
              'outside that range will produce garbage/NaN loss.'.format(X.min(), X.max()))

    print('Fine-tuning for {} epochs...'.format(args.epochs))
    model.train(X, seed=args.seed)

    z = model.predict(X)
    print('latent embedding:', z.shape)

    if args.save_model:
        encoder_path = os.path.join(args.out, 'encoder_dec_finetuned.pt')
        model.save_encoder(encoder_path)  # writes both the weights and the .architecture.csv sidecar
        print('Saved fine-tuned encoder to', encoder_path)

    out_df = pd.DataFrame({'sample_id': d['sample_id'], args.label_col: d['label_classes'][d['y']]})
    for j in range(z.shape[1]):
        out_df['z{}'.format(j)] = z[:, j]
    out_path = os.path.join(args.out, 'finetuned_latent_embedding.csv')
    out_df.to_csv(out_path, index=False)
    print('Saved', out_path, out_df.shape)


if __name__ == '__main__':
    main()
