# Reads whatever *_frozen_cluster.json files exist in results/ (written by
# CNC-DEC/compare_clusters.py, DEC/compare_clusters.py, X-DEC/compare_clusters.py)
# and prints them side by side. Run this from THIS folder (Frozen-Encoder-
# Clustering/) after running one or more of those three scripts - it has no
# model-loading code of its own, so it doesn't hit the cross-package import
# collision the three per-model scripts avoid by each running as a separate
# process from inside its own folder.
import glob
import json

import pandas as pd

rows = []
for path in sorted(glob.glob('results/*_frozen_cluster.json')):
    with open(path) as f:
        rows.append(json.load(f))

if not rows:
    print('No results yet - run CNC-DEC/compare_clusters.py, DEC/compare_clusters.py, '
          'and/or X-DEC/compare_clusters.py first (each from inside its own folder), '
          'then re-run this script.')
else:
    df = pd.DataFrame(rows).set_index('model')
    print(df[['ls', 'n_clusters', 'n_samples', 'ami', 'ari', 'encoder']])
