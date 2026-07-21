# DEC (Plain Deep Embedded Clustering)

The simplest of the four clustering/embedding variants in this project: a
plain (non-variational) MLP autoencoder combined with **Deep Embedded
Clustering (DEC)**. No VAE regularization, no dual-branch input handling —
just reconstruction + clustering. Port of the DAM-IC
`Deep-Embedded-Clustering-generalisability-and-adaptation-for-mixed-data-types`
repo's plain single-input path, rewritten in PyTorch (TensorFlow/Keras has
DLL-loading problems on Windows).

If you want the VAE-regularized version, see [`../CNC-DEC`](../CNC-DEC)
(single concatenated input) or [`../X-DEC`](../X-DEC) (dual-branch,
numeric + binary inputs kept separate).

## How it works

```
[gene expression (min-max scaled) | clinical dummies]   (concatenated into one vector)
        │
        ▼
  Linear → ReLU              (hidden layer, width = --neurons_h)
        │
        ▼
  Linear → ReLU              (bottleneck/latent layer, width = --neurons_e — THIS is the embedding DEC clusters on)
        │
        ▼
  Linear → ReLU → Linear → Sigmoid   (decoder, reconstructs the input, scaled to [0,1])
```

Loss = MSE reconstruction + an L1 penalty on the first hidden layer's
activations (encourages a sparse intermediate representation). See
[`models/autoencoder.py`](models/autoencoder.py) (`_AENet`/`MLPAutoencoder`).

Same two-stage pipeline as CNC-DEC/X-DEC: pretrain the autoencoder on
reconstruction, then jointly fine-tune the encoder + a DEC clustering head
(`models/dec.py`) so the latent space becomes clusterable. See
`../CNC-DEC/README.md` for the DEC algorithm explanation (K-means init →
soft assignment `Q` → sharpened target `P` → KL-divergence gradient steps).

## How to run

From this folder (`DEC/`):

```bash
py run_dec.py --neurons_e 8 --n_clusters 2
```

Or work interactively in [`DEC Runner.ipynb`](DEC%20Runner.ipynb), which
trains and then plots the same 3-panel results as CNC-DEC/X-DEC:
1. PCA scatter of the latent space (predicted clusters vs. true `mgs_level`)
2. Cluster-vs-label crosstab heatmap
3. Cluster-assignment confidence histogram

## What can be changed

| Flag | Meaning |
|---|---|
| `--neurons_e` | **Embedding (latent) layer size** (default 8) — this is the equivalent of CNC-VAE's `--ls`. To try a different latent size, e.g. in the notebook: `MLPAutoencoder(n_features=X.shape[1], neurons_e=4)`. |
| `--neurons_h` | Hidden layer size (default 64). |
| `--epochs` / `--bs` | Autoencoder pretraining epochs (default 500) / batch size (default 64). |
| `--n_clusters` | Number of clusters to identify (default 2) — set to match the number of phenotype groups you're comparing against. |
| `--dec_maxiter` / `--dec_update_interval` / `--dec_batch_size` / `--dec_tol` | DEC self-training loop controls. |
| `--n_init` | K-means restarts for initializing cluster centers. |
| `--clinical_mode` | `'curated'` (default) or `'full'`. |
| `--label_col` | Phenotype column reported against discovered clusters (default `mgs_level`). |
| `--save_model` | Save the trained autoencoder (weights `.pt` + a sidecar `.architecture.csv` — see `models/autoencoder.py:save_encoder`). |

Output directory name encodes `--neurons_h`/`--neurons_e`/`--n_clusters`:
`results/DEC/dec_H{neurons_h}_E{neurons_e}_k{n_clusters}/`.

## Packages needed

```bash
pip install -r requirements.txt
```
installs `torch>=2.0`, `numpy>=1.24`, `pandas>=1.5`, `scikit-learn>=1.2`,
`scipy>=1.10`, `matplotlib>=3.7`.

Same data files as CNC-VAE: `Dataset/aak100_cpmdat.csv` and
`Dataset/MetaSheet_1_4.csv`, joined by `sample_id`.

## Outputs

`run_dec.py` writes to `results/DEC/dec_H{h}_E{e}_k{n}/`:
- `{label_col}.npz` — `z` (latent embedding), `y_pred` (cluster labels), `y_proba`, `centroids`, plus `y`/`sample_id`/`classes`
- `encoder_dec.pt` (with `--save_model`)

Console output reports `acc`/`ari`/`nmi` against `--label_col`, both
periodically during DEC training and at the end.

## Reloading a saved autoencoder

Unlike CNC-VAE/CNC-DEC/X-DEC (which bundle weights + architecture into one
`.pt` file), this autoencoder's `save_encoder()` writes weights to `path` and
the architecture hyperparameters to a **separate sidecar CSV**
(`path + '.architecture.csv'`). Use `models.autoencoder.load_autoencoder(path)`
to reload both together into a ready-to-use `MLPAutoencoder` wrapper. As with
the other variants, the latent size (`neurons_e`) is baked into the saved
weights — you cannot change it without retraining from scratch.
