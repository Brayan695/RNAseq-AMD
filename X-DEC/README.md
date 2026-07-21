# X-DEC (X-shaped VAE + Deep Embedded Clustering)

The dual-branch variant: instead of concatenating clinical + gene-expression
features into one vector upfront (like CNC-VAE/CNC-DEC), X-DEC keeps them as
**two separate inputs** that each get their own encoder branch, matched to
their own data type and reconstruction loss, before merging into a shared
bottleneck. Port of the DAM-IC DEC repo's `xvae.py` / `DECLayer.py` /
`clustering_functions.py` (mixed numerical + binary data types path),
rewritten in PyTorch.

- **SET 1 (s1) must always be numerical** — here, RNA-seq CPM expression — reconstructed with MSE.
- **SET 2 (s2) must always be binary/one-hot** — here, clinical/genotype dummy columns — reconstructed with binary cross-entropy.

## How it works

```
 s1 (numeric: RNA-seq)        s2 (binary: clinical dummies)
        │                              │
        ▼                              ▼
 Linear→BatchNorm→act           Linear→BatchNorm→act      (separate branches, widths --ds1 / --ds2)
        │                              │
        └──────────────┬───────────────┘
                        ▼
              concat → Linear→BatchNorm→act    (shared bottleneck, width --ds12 — branches merge here)
                        │
                 ┌──────┴──────┐
                 ▼             ▼
              fc_mu       fc_logvar             (→ latent z, dimension --ls)
                        │
                        ▼
         shared decoder trunk → splits back into s1-branch / s2-branch
                        │
              ┌─────────┴─────────┐
              ▼                   ▼
       s1 reconstruction    s2 reconstruction (logits → sigmoid)
        (MSE loss)           (BCE-with-logits loss)
```

Model code: [`models/xvae.py`](models/xvae.py) (`_XVAENet`/`xvae` class),
math helpers in [`models/common.py`](models/common.py) (KL/MMD/reparameterization
— same as CNC-VAE's, just shared across a dual-branch architecture).

Same two-stage DEC pipeline as CNC-DEC/DEC (pretrain → joint fine-tune with
a clustering head): see [`models/dec.py`](models/dec.py) (`XDEC` class — same
Student's-t soft-assignment / target-sharpening algorithm as the other two
variants, just calling `xvae.net.encode(x1, x2)` instead of a single-input
`encode(x)`).

## How to run

From this folder (`X-DEC/`):

```bash
py run_xdec.py --ds1 48 --ls 16 --n_clusters 2
```

Or work interactively in [`XDEC Runner.ipynb`](XDEC%20Runner.ipynb), which
trains, saves results + an encoder checkpoint to `results/manual_run/`, then
plots the same 3-panel results as CNC-DEC/DEC (added so all three DEC
variants are visually comparable):
1. PCA scatter of the latent space (predicted clusters vs. true `mgs_level`)
2. Cluster-vs-label crosstab heatmap
3. Cluster-assignment confidence histogram

## What can be changed

| Flag | Meaning |
|---|---|
| `--ls` | Latent dimension size (default 16) — the main knob to experiment with. |
| `--ds1` | Hidden width of the numeric (RNA) branch (default 48). |
| `--ds2` | Hidden width of the binary (clinical) branch (defaults to `--ds1`). |
| `--ds12` | Hidden width of the shared bottleneck where both branches merge (defaults to `--ds1`). |
| `--n_clusters` | Number of clusters to identify (default 2). |
| `--distance` / `--beta` | VAE regularizer (`kl` or `mmd`) and its weight. |
| `--unweighted` | By default, the two branches' reconstruction losses are pooled into one shared per-feature average (a branch with more features gets more weight). Pass `--unweighted` to instead average each branch independently first, so a small-but-important branch (few clinical features) isn't drowned out by a much larger one (many genes). |
| `--epochs` / `--bs` | X-VAE pretraining epochs (default 150) / batch size (default 32). |
| `--dec_maxiter` / `--dec_update_interval` / `--dec_batch_size` / `--dec_tol` | DEC self-training loop controls. |
| `--n_init` | K-means restarts for initializing cluster centers. |
| `--clinical_mode` | `'curated'` (default) or `'full'`. |
| `--label_col` | Phenotype column reported against discovered clusters (default `mgs_level`). |
| `--save_model` | Save the trained encoder as `encoder_xvae.pt` (for the SHAP workflow — see `XVAEEncoderMu` in `models/xvae.py`). |
| `--stability` | Repeated K-fold cluster-stability check (slow — see `models/stability.py`). |

A one-time warning prints if the clinical branch has more columns than there
are samples (`--clinical_mode full` on this cohort size) — a nudge toward
`--clinical_mode curated` rather than a hard failure.

Output directory encodes `--ls`/`--ds1`/`--ds2`/`--ds12`/`--distance`/`--beta`/`--n_clusters`:
`results/XDEC_Clin+mRNA_integration/xdec_LS_{ls}_DS1_{ds1}_DS2_{ds2}_DS12_{ds12}_{distance}_beta_{beta}_k{n_clusters}/`.

## Packages needed

```bash
pip install -r requirements.txt
```
installs `torch>=2.0`, `numpy>=1.24`, `pandas>=1.5`, `scikit-learn>=1.2`,
`scipy>=1.10`, `matplotlib>=3.7` (needed for `XDEC Runner.ipynb`'s plots).

Same data files as CNC-VAE: `Dataset/aak100_cpmdat.csv` and
`Dataset/MetaSheet_1_4.csv`, joined by `sample_id`.

## Outputs

`run_xdec.py` writes to `results/XDEC_Clin+mRNA_integration/xdec_LS_..._k{n}/`:
- `{label_col}.npz` — `z`, `y_pred`, `y_proba`, `centroids`, plus `y`/`sample_id`/`classes`
- `xdec_latent_embedding.csv` — same embedding, drop-in shape for the XGBoost pipeline
- `encoder_xvae.pt` (with `--save_model`)
- `cluster_stability.csv` (with `--stability`)

Console output reports `acc`/`ari`/`nmi` against `--label_col`.

## Reloading a saved encoder

`encoder_xvae.pt` bundles weights + architecture (`s1_input_size`,
`s2_input_size`, `ds1`, `ds2`, `ds12`, `ls`, `dropout`, `act`) into one file —
reload with `models.xvae.load_xvae_model(path)`. As with the other variants,
the latent size is fixed at training time; there's no way to change `--ls`
on a saved checkpoint without retraining.

For tools that expect a single-tensor model input (e.g. SHAP), use
`XVAEEncoderMu` (in `models/xvae.py`) — it wraps the two-branch encoder to
accept one concatenated tensor `[s1 | s2]` and splits it internally.
