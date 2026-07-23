# CNC-DEC (CNC-VAE + Deep Embedded Clustering)

Combines the CNC-VAE encoder (see [`../CNC-VAE`](../CNC-VAE)) with **Deep
Embedded Clustering (DEC)**: pretrains the VAE on reconstruction, then
fine-tunes the encoder *jointly* with a clustering head so the latent space
itself becomes more clusterable. Produces `--n_clusters` clusters directly
from the data — no label is used for training, only (optionally) for
reporting how well the discovered clusters line up with a known phenotype
afterward.

This is the **single-concatenated-input** variant (clinical + gene
expression concatenated into one vector, same as CNC-VAE). For the
dual-branch variant that keeps the two data types as separate model inputs,
see [`../X-DEC`](../X-DEC).

## How it works

Two-stage process, both run automatically by one call to `run_cncdec()`:

1. **Pretraining** — an ordinary CNC-VAE (see `../CNC-VAE/README.md` for the
   architecture diagram) is trained on reconstruction only, giving the
   encoder a reasonable starting latent space.
2. **DEC self-training** — K-means initializes cluster centers in that latent
   space; then, in a loop, the model computes soft cluster assignments `Q`
   (via a Student's-t kernel, same kernel t-SNE uses), sharpens that into a
   target distribution `P`, and takes a gradient step pulling `Q` toward `P`.
   This is repeated until few samples change cluster assignment between
   updates (`--dec_tol`) or `--dec_maxiter` is reached. Critically, **the
   encoder's weights are updated too** — not just the cluster centers — so
   the latent space itself is reshaped to be more clusterable, not just fit
   with a fixed embedding.

Model code:
- [`models/cncvae.py`](models/cncvae.py) — same architecture as `../CNC-VAE/models/cncvae.py`, just built from constructor kwargs (`CNCVAE(input_size=..., ds=..., ls=...)`) instead of an `argparse.Namespace`, so it's usable standalone.
- [`models/dec.py`](models/dec.py) — the `DECHead` clustering layer and the `DEC` class that runs the self-training loop.
- [`models/metrics.py`](models/metrics.py) — `calculate_metrics()`: clustering accuracy (via Hungarian matching), NMI, ARI — used only for reporting against a known label, never for training.

## How to run

From this folder (`CNC-DEC/`):

```bash
py run_cncdec.py --ds 64 --ls 16 --n_clusters 2
```

Or work interactively in [`CNC-DEC Runner.ipynb`](CNC-DEC%20Runner.ipynb),
which trains, saves results to `results/manual_run/`, then plots:
1. A 2-panel PCA scatter of the latent space (predicted clusters vs. true `mgs_level`)
2. A **cluster-vs-label crosstab heatmap** (counts of each cluster × each `mgs_level` value — the confusion-matrix-style view)
3. A histogram of cluster-assignment confidence (`max soft-assignment probability` per sample)

## What can be changed

All via command-line flags on `run_cncdec.py` — no code editing needed:

| Flag | Meaning |
|---|---|
| `--ls` | Latent dimension size (default 16) — the main knob to experiment with, same as CNC-VAE. |
| `--ds` | Encoder/decoder hidden layer width (default 64). |
| `--n_clusters` | Number of clusters DEC should find (default 2) — set this to match how many phenotype groups you're comparing against in the confusion matrix. |
| `--distance` / `--beta` | VAE pretraining regularizer (`kl` or `mmd`) and its weight. |
| `--epochs` / `--bs` | CNC-VAE pretraining epochs / batch size. |
| `--dec_maxiter` / `--dec_update_interval` / `--dec_batch_size` / `--dec_tol` | DEC self-training loop controls. |
| `--n_init` | Number of K-means restarts used to initialize cluster centers. |
| `--clinical_mode` | `'curated'` (default, recommended for this cohort size) or `'full'`. |
| `--label_col` | Phenotype column reported against the discovered clusters (never used in training). Default `mgs_level`. |
| `--save_model` | Save the trained encoder as `encoder_cncvae.pt` (same checkpoint format as CNC-VAE's `save_encoder`). |
| `--stability` | Also run a repeated K-fold cluster-stability check (retrains `k * rep` times — slow; see `models/stability.py`). |

Output directory name already encodes `--ls`/`--ds`/`--distance`/`--beta`/
`--n_clusters`, so different configurations don't overwrite each other:
`results/CNCDEC_Clin+mRNA_integration/cncdec_LS_{ls}_DS_{ds}_{distance}_beta_{beta}_k{n_clusters}/`.

## Packages needed

```bash
pip install -r requirements.txt
```
installs `torch>=2.0`, `numpy>=1.24`, `pandas>=1.5`, `scikit-learn>=1.2`,
`scipy>=1.10` (for the Hungarian-matching accuracy metric), `matplotlib>=3.7`
(needed for `CNC-DEC Runner.ipynb`'s plots).

Same data files as CNC-VAE: `Dataset/aak100_cpmdat.csv` and
`Dataset/MetaSheet_1_4.csv`, joined by `sample_id`.

## Outputs

`run_cncdec.py` writes to `results/CNCDEC_Clin+mRNA_integration/cncdec_LS_..._k{n}/`:
- `{label_col}.npz` — latent embedding `z`, `y_pred` (cluster labels), `y_proba` (soft assignments), `centroids`, plus `y`/`sample_id`/`classes` for interpretation
- `cncdec_latent_embedding.csv` — same embedding in the CSV shape the XGBoost pipeline expects
- `encoder_cncvae.pt` (with `--save_model`)
- `cluster_stability.csv` (with `--stability`)

Console output reports `acc`/`ari`/`nmi` against `--label_col` both
periodically during DEC training and at the end.
