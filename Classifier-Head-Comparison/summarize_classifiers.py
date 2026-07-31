# Reads whatever *_classifier_result.json files exist in each subfolder's
# classifier_results/ (written by CNC-VAE/classify_latent.py,
# CNC-DEC/classify_latent.py, DEC/classify_latent.py, X-DEC/classify_latent.py)
# and prints them side by side. Run this from THIS folder (Classifier-Head-
# Comparison/) after running one or more of those four scripts.
import glob
import json

import pandas as pd

rows = []
# glob('**', recursive=True) also catches results saved to subfolders (e.g.
# --out classifier_results/ls16), not just directly under classifier_results/.
for path in sorted(glob.glob('*/classifier_results/**/*_classifier_result.json', recursive=True)):
    with open(path) as f:
        rows.append(json.load(f))

if not rows:
    print('No results yet - run classify_latent.py in one or more of CNC-VAE/, CNC-DEC/, '
          'DEC/, X-DEC/ first (each from inside its own folder), then re-run this script.')
else:
    df = pd.DataFrame(rows).sort_values(['model', 'latent_dim'], ascending=[True, False])
    df = df.set_index(['model', 'latent_dim'])
    print(df[['clin_file', 'n_splits', 'n_samples', 'pr_auc', 'roc_auc', 'encoder']])
