"""
Compares X-DEC's and CNC-DEC's cluster-stability results (each produced by
running `run_xdec.py --stability` / `run_cncdec.py --stability` in their own
folders, which saves a `cluster_stability.csv` per model). This script only
loads and compares those already-computed results - it does not re-run any
training itself, since ../X-DEC/models and ../CNC-DEC/models can't both be
imported as `models` in one process (same package name, different folders).

Reports mean stability per model and a per-sample side-by-side comparison
(which samples are stable in one model but not the other).
"""
import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--xdec_stability', required=True,
                     help='Path to X-DEC\'s cluster_stability.csv (from run_xdec.py --stability)')
parser.add_argument('--cncdec_stability', required=True,
                     help='Path to CNC-DEC\'s cluster_stability.csv (from run_cncdec.py --stability)')
parser.add_argument('--out', default='comparison', help='Output directory for plots/tables')


def main():
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    xdec = pd.read_csv(args.xdec_stability).set_index('sample_id')
    cncdec = pd.read_csv(args.cncdec_stability).set_index('sample_id')

    common = xdec.index.intersection(cncdec.index)
    if len(common) == 0:
        raise ValueError('No overlapping sample_id between the two stability files.')
    if len(common) < len(xdec) or len(common) < len(cncdec):
        print('WARNING: only {} of {} (X-DEC) / {} (CNC-DEC) samples overlap.'.format(
            len(common), len(xdec), len(cncdec)))
    xdec = xdec.loc[common]
    cncdec = cncdec.loc[common]

    print('X-DEC   mean cluster stability: {:.3f} (median {:.3f})'.format(
        xdec['stability'].mean(), xdec['stability'].median()))
    print('CNC-DEC mean cluster stability: {:.3f} (median {:.3f})'.format(
        cncdec['stability'].mean(), cncdec['stability'].median()))

    combined = pd.DataFrame({
        'xdec_cluster': xdec['cluster'], 'xdec_stability': xdec['stability'],
        'cncdec_cluster': cncdec['cluster'], 'cncdec_stability': cncdec['stability'],
    }, index=common)
    combined['stability_gap'] = combined['xdec_stability'] - combined['cncdec_stability']
    out_path = os.path.join(args.out, 'stability_comparison.csv')
    combined.to_csv(out_path)
    print('Saved', out_path)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))

    axes[0].boxplot([xdec['stability'].to_numpy(), cncdec['stability'].to_numpy()],
                    labels=['X-DEC', 'CNC-DEC'])
    axes[0].set_ylabel('per-sample cluster stability')
    axes[0].set_title('Cluster stability distribution')

    axes[1].scatter(combined['xdec_stability'], combined['cncdec_stability'], alpha=0.6, s=20)
    lims = [0, 1]
    axes[1].plot(lims, lims, linestyle='--', color='gray', linewidth=1)
    axes[1].set_xlabel('X-DEC stability')
    axes[1].set_ylabel('CNC-DEC stability')
    axes[1].set_xlim(lims)
    axes[1].set_ylim(lims)
    axes[1].set_title('Per-sample stability: X-DEC vs. CNC-DEC')

    fig.tight_layout()
    fig_path = os.path.join(args.out, 'stability_comparison.png')
    fig.savefig(fig_path, bbox_inches='tight')
    print('Saved', fig_path)


if __name__ == '__main__':
    main()
