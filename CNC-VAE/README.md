# CNC-VAE (Clin+mRNA Integration)

A PyTorch port of the CancerAI-CL **CNC-VAE** ("Concatenated iNputs VAE"): a
single-input variational autoencoder that concatenates clinical/genotype
features with gene-expression (RNA-seq CPM) features into one vector,
compresses it down to a small **latent embedding**, and reconstructs the
input from that embedding. The trained latent embedding is what downstream
scripts (classifiers, t-SNE, SHAP) use as a compact per-patient
representation.

Ported to PyTorch (rather than the original TensorFlow/Keras) to avoid
TensorFlow's DLL-loading problems on Windows.

## How it works

```
[clinical features | gene expression features]  (concatenated into one vector)
        │
        ▼
  Linear → BatchNorm → activation           (encoder hidden layer, width = --ds)
        │
   ┌────┴────┐
   ▼         ▼
 fc_mu   fc_logvar        (each Linear(ds → ls) — the latent mean/log-variance heads)
   │         │
   └────┬────┘
        ▼
  reparameterization trick → z   (the latent embedding, dimension = --ls)
        │
        ▼
  Linear → BatchNorm → activation → Dropout → Linear   (decoder, mirrors the encoder)
        │
        ▼
  reconstructed [clinical | gene expression] vector
```

Loss = MSE reconstruction error + `beta * distance`, where `distance` is
either KL divergence to a standard normal prior (`--distance kl`) or Maximum
Mean Discrepancy (`--distance mmd`). `beta` controls how much the latent
space is pulled toward the prior vs. how accurately the input is
reconstructed (beta-VAE).

Core model code: [`models/cncvae.py`](models/cncvae.py) — `_CNCVAENet` is the
raw PyTorch network, `CNCVAE` is a wrapper that owns training/prediction.

## The 4 scripts in this folder

| Script | Purpose |
|---|---|
| [`run_cncvae.py`](run_cncvae.py) | Command-line training: whole-cohort or stratified-CV-fold mode. Saves the latent embedding as `.npz`. |
| [`compare_representations.py`](compare_representations.py) | Reproduces the CNC_VAE_AMD_V3 reference notebook: auto-tunes beta, compares CNC-VAE latent vs. PCA vs. raw features with 5-fold CV classifiers, plots t-SNE, saves the embedding as `.csv`. |
| [`analyse_representations.py`](analyse_representations.py) | Reads the `.npz` files `run_cncvae.py` wrote and scores simple classifiers (logreg/nb/svm/rf) on them, optionally with a t-SNE plot. Does **not** retrain anything. |
| [`load_encoder.py`](load_encoder.py) | Reloads a saved encoder checkpoint (`encoder_cncvae*.pt`) and re-encodes data into its latent space **without retraining** — e.g. to regenerate an embedding CSV from an old model. |

## What can be changed

Everything below is a command-line flag — no code editing required for any
of it.

**Latent size (the main thing to experiment with):**
```bash
py run_cncvae.py --ds 64 --ls 8    # try 8-dimensional latent space
py run_cncvae.py --ds 64 --ls 4    # try 4-dimensional latent space
```
`--ls` sets the width of the bottleneck (`fc_mu`/`fc_logvar` in
`models/cncvae.py`). Different `--ls` values are safe to run back-to-back —
`run_cncvae.py`'s output directory name already encodes `--ls`/`--ds`/
`--distance`/`--beta` so runs don't overwrite each other. **`compare_representations.py`
does NOT do this** — its `--out` directory holds fixed filenames
(`comparison_table.csv`, `cncvae_latent_embedding.csv`, etc.), so give each
run of that script its own `--out` folder if you don't want to overwrite the
previous run's results.

**Other hyperparameters (all via flags):**
- `--ds` — hidden/dense layer width (encoder and decoder share this width). `run_cncvae.py` only accepts 16/32/64/128 (each preset also fixes epochs/batch-size/dropout/activation — see the `configs` dict at the top of the file); `compare_representations.py` takes any int and exposes epochs/bs/dropout as separate flags.
- `--beta` — beta-VAE regularization weight (or let `compare_representations.py` auto-resolve it via `--target_ratio`).
- `--distance` — `kl` or `mmd`.
- `--clinical_mode` — `'full'` (every clinical column) or `'curated'` (age/sex/genotype + prevalence-filtered, leakage-excluded history flags).
- `--fold` (`run_cncvae.py` only) — `0` for whole-cohort training, or `1..--numfolds` for one stratified CV fold (reports val_loss against the held-out fold).
- `--label_col` — phenotype column used for CV stratification / reporting (default `mgs_level`).
- `--save_model` — also save the encoder as a standalone, reloadable `.pt` checkpoint (see `load_encoder.py` below).

## Packages needed

```bash
pip install -r requirements.txt
```
which installs: `torch>=2.0`, `numpy>=1.24`, `pandas>=1.5`, `scikit-learn>=1.2`,
`matplotlib>=3.7`, `shap>=0.44`.

Data files expected (relative paths resolved automatically by
`misc/dataset.py`, two directories up into a sibling `Dataset/` folder):
- `Dataset/aak100_cpmdat.csv` — RNA-seq CPM expression matrix
- `Dataset/MetaSheet_1_4.csv` — clinical/genotype metadata

Both are joined by `sample_id`; only samples present in both files are used.

## How to run

From this folder (`CNC-VAE/`):

```bash
# Train on the whole cohort, latent size 16, hidden width 64
py run_cncvae.py --ds 64 --ls 16

# Train + evaluate one stratified CV fold (fold 1 of 5)
py run_cncvae.py --ds 64 --ls 16 --fold 1

# Reproduce the reference-notebook comparison table + t-SNE + embedding CSV
py compare_representations.py --ds 256 --ls 16 --epochs 150 --bs 16 --out comparison

# Score the embeddings run_cncvae.py --fold 0 produced
py analyse_representations.py --resdir results/CNCVAE_Clin+mRNA_integration/cncvae_LS_16_DS_64_kl_beta_1.0 --numfolds 0

# Reload a saved encoder and re-encode without retraining
py load_encoder.py --encoder comparison/encoder_cncvae.pt --clinical_mode curated --normalization log1p_minmax
```

## Outputs

- `run_cncvae.py` → `results/CNCVAE_Clin+mRNA_integration/cncvae_LS_{ls}_DS_{ds}_{distance}_beta_{beta}/` containing `{label_col}.npz` (or `{label_col}{fold}.npz` per fold), `{label_col}_labels.npz`, and (with `--save_model`) `encoder_cncvae.pt` / `vae_cncvae.pt`.
- `compare_representations.py` → `--out` folder containing `comparison_table.csv`, `tsne_comparison.png`, `cncvae_latent_embedding.csv`, and (with `--save_model`) `encoder_cncvae.pt`.
- `analyse_representations.py` → `--out` folder (default `analyses/`) containing a results CSV and optionally a t-SNE PNG.

## Reloading a saved encoder (`load_encoder.py`)

The **latent size cannot be changed on a saved model** — it's baked into the
weight matrix shapes at training time (`fc_mu`/`fc_logvar`: `ds→ls`,
decoder's first layer: `ls→ds`). To try a different `--ls`, retrain from
scratch; there's no shortcut via a checkpoint.

What a saved checkpoint (`encoder_cncvae*.pt`) IS good for is re-encoding
data **without retraining**, using the exact architecture it was saved with.
`load_encoder.py` does this — but you must tell it which preprocessing the
checkpoint was originally trained with via `--normalization`:
- `minmax` — if trained via `run_cncvae.py` (`misc.helpers.normalizeRNA`)
- `log1p_minmax` — if trained via `compare_representations.py` (log1p then min-max)

Getting this wrong won't error (the input width still matches), it will just
silently produce meaningless embeddings — so double check which script
produced the checkpoint you're loading.
