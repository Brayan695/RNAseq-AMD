# MEGENA_AMI — Gene Co-expression Network Analysis (AMI distance)

Builds a **gene co-expression network** from the AMD RNA-seq cohort using
**MEGENA** (Multiscale Embedded Gene Co-Expression Network Analysis), then
extracts per-module "eigengenes" and tests whether they differ between MGS
grades. "AMI" = **Adjusted Mutual Information** — the gene-gene similarity
measure used here to build the network (computed upstream in
`../Python Scripts/AMI.ipynb`), as opposed to Pearson correlation, JSD, or
Wasserstein distance used in the sibling `MEGENA_BD` / `MEGENA_SKL_JSD` /
`MEGENA_WD` / `MEGENA_WGCNA` folders (not covered by this README).

**⚠️ Before running anything:** all three scripts use **hardcoded absolute
paths** (`C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/...`). On a new
machine, find-and-replace that path prefix with wherever this repo lives
locally — the scripts will error immediately otherwise (file not found).

## The 3 scripts in this folder

| Script | Cohort | Distance file | What it adds |
|---|---|---|---|
| [`MEGENA_AMI_1_100.R`](MEGENA_AMI_1_100.R) | MGS1 only (mild) | `Dataset/AMI/ami_edges_control.csv` | Builds a network + modules for MGS1 alone, reports per-module eigengene summary stats (mean/SD/median/min/max). |
| [`MEGENA_AMI_4_100.R`](MEGENA_AMI_4_100.R) | MGS4 only (severe) | `Dataset/AMI/ami_edges_late.csv` | Identical to the above, for MGS4. |
| [`MEGENA_AMI_1_4_100.R`](MEGENA_AMI_1_4_100.R) | MGS1 + MGS4 combined | `Dataset/AMI/ami_edges.csv` | Builds one network on the combined cohort, **plus** tests each module's eigengene for a difference between MGS1 vs. MGS4 (Wilcoxon rank-sum + Student's t-test, BH-corrected). |

The `_1_100`/`_4_100` pair differ from each other in exactly 3 lines (input
file, plot title, `mgs_level` filter) — diffing them is the fastest way to
see what changes between a "single group" run and its counterpart.

## How it works (pipeline, same shape in all 3 scripts)

```r
# 1. Load the precomputed AMI gene-gene similarity edges + expression + gene metadata
distance = read.csv(".../Dataset/AMI/ami_edges*.csv")   # columns: from, to, weight
genes    = read.csv(".../Dataset/aak100_cpmdat.csv")     # samples x genes CPM matrix
info     = read.delim(".../Dataset/gene_info.tsv")       # ensembl_gene_id -> external_gene_name etc.

# 2. Build a Planar Filtered Network (PFN) - MEGENA's network-sparsification step
pfn = calculate.PFN(distance)
pfn_g = igraph::graph_from_data_frame(d = pfn, directed = FALSE)

# 3. Run MEGENA's multiscale clustering to find modules (+ hub genes) on that network
meg = do.MEGENA(g = pfn_g, mod.pval = 0.05, hub.pval = 0.05, min.size = 10, n.perm = 100)
module_summary = MEGENA.ModuleSummary(meg)
modules = module_summary$modules   # named list: module id -> vector of gene ids

# 4. Visualize interactively (nodes colored/grouped by module, hover = gene + module + degree)
visNetwork(nodes, edges, main = "...") %>% visIgraphLayout() %>% visOptions(...) %>% visLegend()

# 5. Score each module: modularity, conductance, density, transitivity, avg. intra-module correlation

# 6. Compute one "eigengene" per module = PC1 of that module's genes' expression (prcomp)
#    -> a single per-sample number summarizing that module's expression level

# 7. (1_4_100 only) Test each module eigengene for a difference between MGS1 vs MGS4
```

## What can be changed

All of these are plain variable edits at the top of the script (no
command-line flags — these are interactive R scripts, run top-to-bottom in
RStudio):

| What | Where | Notes |
|---|---|---|
| Input distance file | `distance = read.csv(...)` | Swap in a different `Dataset/AMI/ami_edges*.csv` variant (there are several precomputed: `_QNclass`, `_QNglobal`, `_QS`, etc. — see `Dataset/AMI/`) to compare normalization choices. |
| `mod.pval` / `hub.pval` | `do.MEGENA(...)` call | Significance thresholds for calling something a module / a hub gene. Lower = stricter (fewer, more confident modules/hubs). |
| `min.size` | `do.MEGENA(...)` call | Minimum genes per module — modules smaller than this are dropped. |
| `n.perm` | `do.MEGENA(...)` call | Number of permutations for the significance tests above — higher = more stable p-values, slower. |
| Which MGS grade | `MEGENA_AMI_1_100.R` / `_4_100.R`: `subset(genes, mgs_level == "MGS1")` | Change to any other grade present in `mgs_level` if the cohort has more than MGS1/MGS4. |
| Paired vs. unpaired testing | `MEGENA_AMI_1_4_100.R`, the commented-out `pair_col` block (~line 213) | By default there's no donor/pair ID, so the script automatically falls back to an **unpaired** Wilcoxon rank-sum test. If your data has a pairing column (e.g. same donor sampled twice), uncomment that block and set `pair_col` to enable a **paired** signed-rank test instead. |
| Plots | The commented-out blocks at the bottom of `_1_4_100.R` (correlation heatmap, per-module boxplots) | Left in as ready-to-uncomment extensions, not run by default. |

## Packages needed

```r
install.packages(c(
  "MEGENA",        # co-expression network module detection - see note below
  "igraph",        # graph object + graph metrics (modularity, conductance, density, transitivity)
  "visNetwork",     # interactive network visualization
  "RColorBrewer",   # color palettes for the plots
  "readr",          # fast CSV/TSV reading
  "dplyr",          # data wrangling (group_by/summarise, used in the paired-test path)
  "ggplot2",        # plotting (only used in the commented-out extension blocks)
  "pheatmap",       # heatmap plotting (only used in the commented-out extension blocks)
  "reshape2",       # long-format reshaping (only used in the commented-out extension blocks)
  "tidyr"           # pivot_wider, used via tidyr:: in the paired-test path
))
```

**MEGENA** is on CRAN (`install.packages("MEGENA")` should just work). If
that ever fails, the development version is on GitHub:
```r
install.packages("devtools")
devtools::install_github("songw01/MEGENA")
```
Source: [MEGENA @ CRAN/METACRAN](https://www.r-pkg.org/pkg/MEGENA), [songw01/MEGENA on GitHub](https://github.com/songw01/MEGENA).

## Data files needed

- `Dataset/aak100_cpmdat.csv` — RNA-seq CPM expression matrix, samples × genes, with an `mgs_level` phenotype column (`MGS1`/`MGS4`/...).
- `Dataset/gene_info.tsv` — Ensembl gene ID → gene symbol / genomic location lookup, used to label network nodes.
- `Dataset/AMI/ami_edges.csv`, `ami_edges_control.csv`, `ami_edges_late.csv` — precomputed gene-gene AMI similarity edge lists (`from, to, weight`), one per cohort variant. These are generated upstream by **`../Python Scripts/AMI.ipynb`** — if you need to regenerate them (new samples, different gene set, etc.), that notebook is the place to start; it's a separate Python workflow, not part of this R folder.

## How to run

Open the `.Rproj` in `Scripts/` (or just open the `.R` file in RStudio),
**update the hardcoded paths** at the top to match your machine, then
source top-to-bottom — either the whole file (`Ctrl+Alt+R` / *Source*) or
section-by-section if you want to inspect intermediate objects
(`pfn_g`, `meg`, `module_summary`, `ME_df`) along the way:

```r
# from an R console, with the working directory set anywhere (paths are absolute):
source("MEGENA_AMI_1_4_100.R")   # combined MGS1+MGS4 cohort + differential eigengene test
source("MEGENA_AMI_1_100.R")     # MGS1-only network
source("MEGENA_AMI_4_100.R")     # MGS4-only network
```

Each script is meant to be run **interactively**, not headlessly — the
`visNetwork(...)` call opens an interactive HTML widget in the RStudio
Viewer pane, and several `print(...)` calls are how you're meant to inspect
results (module metrics table, differential-eigengene results table), not
values that get saved to disk automatically.

## Outputs (in-session objects, not auto-saved to disk)

| Object | What it is |
|---|---|
| `pfn_g` | The `igraph` network object (Planar Filtered Network). |
| `meg` / `module_summary` | Raw and summarized MEGENA output; `module_summary$modules` is the named list of gene modules. |
| `module_metrics` | Data frame: per-module modularity/conductance/density/transitivity/avg. correlation. |
| `ME_df` | Per-sample eigengene matrix (one column per module, PC1 of that module's genes) + `Phenotype`. |
| `module_results` | (1_4_100 only) Per-module Wilcoxon + t-test results vs. MGS1/MGS4, with BH-adjusted FDR — this is the table to look at for "which modules differ between mild and severe AMD." |

If you want these saved to disk rather than just left in the R session,
add e.g. `write.csv(module_results, "module_results.csv")` — none of the
three scripts currently do this by default.
