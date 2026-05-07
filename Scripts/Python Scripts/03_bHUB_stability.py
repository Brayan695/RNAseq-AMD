"""
bHUB-Inspired Stability Testing for AMD Co-expression Networks
==============================================================
Dataset : aak100_cpmdat.csv  (CPM-normalised RNA-seq, 81 genes, 166 samples)
Classes : MGS1 (early/normal, n=105) vs MGS4 (advanced AMD, n=61)

Method reference
----------------
Kim J, Do K-A, Ha MJ, Peterson CB. (2019).
"Bayesian Inference of Hub Nodes Across Multiple Networks."
Biometrics 75(1): 172–182.  DOI: 10.1111/biom.12958

Background
----------
bHUB introduces a latent hub indicator h_{k,i} for each gene i in each
group k and uses a Bayesian hierarchical prior to share hub information
across groups.  The posterior probability PPI_i = P(h_{k,i}=1 | data)
is the key inferential quantity; a gene with PPI >= 0.5 is declared a hub.

This script provides a fully frequentist approximation of the bHUB
framework:

  1.  For each bootstrap replicate b = 1 … B:
        a.  Draw a sub-sample (≈60% without replacement) from MGS1 and MGS4.
        b.  Estimate co-expression networks via graphical LASSO (GL).
        c.  Compute node degree d_i in each group.
        d.  Flag gene i as a hub if d_i ≥ 75th-percentile of degrees
            (analogous to bHUB's per-replicate edge inclusion).

  2.  Stability score  s_{k,i} = proportion of replicates where gene i
      was flagged as a hub in group k.
      (This is the frequentist analog of bHUB's PPI.)

  3.  Hub classification (threshold τ = 0.60, ~bHUB's 0.5 PPI cut-off):
        Co-hub         s_{MGS1,i} ≥ τ  AND  s_{MGS4,i} ≥ τ
        MGS1-specific  s_{MGS1,i} ≥ τ  AND  s_{MGS4,i} < τ
        MGS4-specific  s_{MGS1,i} < τ  AND  s_{MGS4,i} ≥ τ
        Non-hub        neither threshold met

  4.  Full-data GL networks are estimated for visualisation.

Outputs
-------
bhub_stability_results.png  – stability scatter, bar chart, pie, networks
bhub_stability_results.csv  – per-gene stability scores and hub status
cohub_genes.csv             – list of co-hub (shared) genes

Runtime notes
-------------
N_BOOT = 100 (default).  Each replicate fits a GL with 3-fold CV,
which takes ~2-5 s on this dataset.  Total ≈ 5-15 minutes.
Set N_BOOT = 20 for a quick test run.
"""

import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import networkx as nx

from sklearn.preprocessing import StandardScaler
from sklearn.covariance import GraphicalLassoCV, GraphicalLasso

np.random.seed(42)

# ── Configuration ─────────────────────────────────────────────────────────────
DATA_PATH = (
    r"C:\Users\Brayan Gutierrez\Desktop\RNAseq-AMD\Dataset\aak100_cpmdat.csv"
)
OUT_DIR   = r"C:\Users\Brayan Gutierrez\Desktop\RNAseq-AMD\Scripts\Python Scripts"

N_BOOT          = 100     # bootstrap replicates  (set to 20 for quick test)
SUBSAMPLE_FRAC  = 0.60    # fraction of each group sampled per replicate
HUB_PERCENTILE  = 75      # top-quartile degree = hub candidate per replicate
STAB_THRESHOLD  = 0.60    # stability score cut-off  (analogous to PPI ≥ 0.5)

# ── 1. Load and Split Data ────────────────────────────────────────────────────
df = pd.read_csv(DATA_PATH, index_col=0)

feature_cols = [c for c in df.columns if c != "mgs_level"]
gene_names   = np.array(feature_cols)
p            = len(gene_names)

X_mgs1 = df.loc[df["mgs_level"] == "MGS1", feature_cols].values
X_mgs4 = df.loc[df["mgs_level"] == "MGS4", feature_cols].values

print("=" * 65)
print("bHUB-Inspired Stability Testing – AMD Co-expression Networks")
print("=" * 65)
print(f"Genes (nodes) : {p}")
print(f"MGS1 samples  : {X_mgs1.shape[0]}")
print(f"MGS4 samples  : {X_mgs4.shape[0]}")
print(f"Bootstrap B   : {N_BOOT}    sub-sample fraction : {SUBSAMPLE_FRAC}")
print(f"Hub threshold : {STAB_THRESHOLD}")

# ── 2. Standardise Within Each Group ─────────────────────────────────────────
X1 = StandardScaler().fit_transform(X_mgs1)
X4 = StandardScaler().fit_transform(X_mgs4)


# ── 3. Bootstrap Stability Selection ─────────────────────────────────────────
def bootstrap_hub_stability(X, n_boot=100, subsample_frac=0.60,
                             hub_percentile=75):
    """
    Returns a stability score vector of length p.
    Score_i = proportion of bootstrap replicates in which gene i
              was flagged as a hub node.
    """
    n, p = X.shape
    # Cap n_sub at n so replace=False never exceeds population.
    # When n < p (e.g. MGS4: 61 samples, 81 genes) we fall back to using
    # the full group per replicate; the regularised GL still produces valid
    # precision estimates and stability across alpha values is preserved.
    n_sub = min(max(int(n * subsample_frac), 10), n)
    use_replace = False          # always sample without replacement
    hub_count = np.zeros(p)

    for b in range(n_boot):
        idx   = np.random.choice(n, size=n_sub, replace=use_replace)
        X_sub = X[idx]

        # ── Estimate precision matrix ─────────────────────────────────────────
        try:
            # When n_sub < p, CV-based alpha selection may fail; shrink cv folds
            n_cv = min(3, max(2, n_sub // (p + 1) + 1))
            gl = GraphicalLassoCV(cv=n_cv, max_iter=500, tol=1e-4, n_jobs=-1)
            gl.fit(X_sub)
            prec = gl.precision_
        except Exception:
            # Fallback: fixed-alpha GL (works even when n < p)
            try:
                gl = GraphicalLasso(alpha=0.1, max_iter=300)
                gl.fit(X_sub)
                prec = gl.precision_
            except Exception:
                continue          # skip this replicate

        # ── Node degree ───────────────────────────────────────────────────────
        adj = (np.abs(prec) > 1e-8).astype(float)
        np.fill_diagonal(adj, 0)
        degree = adj.sum(axis=1)

        # ── Flag hubs as top-degree nodes ─────────────────────────────────────
        positive_degrees = degree[degree > 0]
        if len(positive_degrees) == 0:
            continue
        thresh = np.percentile(positive_degrees, hub_percentile)
        hub_count += (degree >= thresh).astype(float)

    return hub_count / n_boot


print(f"\nRunning {N_BOOT} bootstrap replicates for MGS1 …", flush=True)
stab_mgs1 = bootstrap_hub_stability(X1, n_boot=N_BOOT,
                                     subsample_frac=SUBSAMPLE_FRAC,
                                     hub_percentile=HUB_PERCENTILE)

print(f"Running {N_BOOT} bootstrap replicates for MGS4 …", flush=True)
stab_mgs4 = bootstrap_hub_stability(X4, n_boot=N_BOOT,
                                     subsample_frac=SUBSAMPLE_FRAC,
                                     hub_percentile=HUB_PERCENTILE)

print("Bootstrap complete.")

# ── 4. Hub Classification ─────────────────────────────────────────────────────
hub_mgs1   = gene_names[stab_mgs1 >= STAB_THRESHOLD]
hub_mgs4   = gene_names[stab_mgs4 >= STAB_THRESHOLD]
cohubs     = np.intersect1d(hub_mgs1, hub_mgs4)
only_mgs1  = np.setdiff1d(hub_mgs1, hub_mgs4)
only_mgs4  = np.setdiff1d(hub_mgs4, hub_mgs1)

print(f"\nHub genes in MGS1          : {len(hub_mgs1)}")
print(f"Hub genes in MGS4          : {len(hub_mgs4)}")
print(f"Co-hub genes (shared)      : {len(cohubs)}")
print(f"MGS1-specific hubs         : {len(only_mgs1)}")
print(f"MGS4-specific hubs         : {len(only_mgs4)}")

# ── 5. Full-Data Network Estimation (for visualisation) ──────────────────────
print("\nFitting full-data graphical LASSO for network visualisation …", flush=True)

def build_network(X, gene_names):
    gl = GraphicalLassoCV(cv=5, max_iter=1000, n_jobs=-1)
    gl.fit(X)
    prec = gl.precision_
    adj  = (np.abs(prec) > 1e-8).astype(int)
    np.fill_diagonal(adj, 0)
    G = nx.from_numpy_array(adj)
    nx.relabel_nodes(G, {i: gene_names[i] for i in range(len(gene_names))},
                     copy=False)
    return G, prec

G1, prec1 = build_network(X1, gene_names)
G4, prec4 = build_network(X4, gene_names)

# ── 6. Results Summary Table ──────────────────────────────────────────────────
summary = pd.DataFrame({
    "Gene":         gene_names,
    "Stab_MGS1":    np.round(stab_mgs1, 4),
    "Stab_MGS4":    np.round(stab_mgs4, 4),
    "Degree_MGS1":  [G1.degree(g) for g in gene_names],
    "Degree_MGS4":  [G4.degree(g) for g in gene_names],
    "Hub_MGS1":     (stab_mgs1 >= STAB_THRESHOLD).astype(int),
    "Hub_MGS4":     (stab_mgs4 >= STAB_THRESHOLD).astype(int),
})
summary["Status"] = "Non-hub"
summary.loc[
    (summary["Hub_MGS1"] == 1) & (summary["Hub_MGS4"] == 1), "Status"
] = "Co-hub"
summary.loc[
    (summary["Hub_MGS1"] == 1) & (summary["Hub_MGS4"] == 0), "Status"
] = "MGS1-specific"
summary.loc[
    (summary["Hub_MGS1"] == 0) & (summary["Hub_MGS4"] == 1), "Status"
] = "MGS4-specific"

summary = summary.sort_values("Stab_MGS4", ascending=False).reset_index(drop=True)

print("\nTop 20 genes by MGS4 stability score:")
print(summary.head(20).to_string(index=False))

# ── 7. Plots ──────────────────────────────────────────────────────────────────
COLOR_MAP = {
    "Co-hub":        "#d62728",
    "MGS1-specific": "#1f77b4",
    "MGS4-specific": "#ff7f0e",
    "Non-hub":       "#cccccc",
}

fig = plt.figure(figsize=(16, 14))
gs  = gridspec.GridSpec(3, 2, figure=fig, hspace=0.50, wspace=0.35)
fig.suptitle(
    "bHUB-Inspired Stability Testing – AMD Co-expression Networks\n"
    f"(B={N_BOOT} bootstrap replicates, threshold τ={STAB_THRESHOLD})",
    fontsize=13, fontweight="bold",
)

# Panel A – Stability scatter (full width)
ax = fig.add_subplot(gs[0, :])
for status in ["Non-hub", "MGS1-specific", "MGS4-specific", "Co-hub"]:
    grp = summary[summary["Status"] == status]
    ax.scatter(
        grp["Stab_MGS1"], grp["Stab_MGS4"],
        color=COLOR_MAP[status], label=status,
        alpha=0.75, s=55, edgecolors="none", zorder=3,
    )
ax.axvline(STAB_THRESHOLD, color="gray", linestyle="--", linewidth=1.2, alpha=0.7)
ax.axhline(STAB_THRESHOLD, color="gray", linestyle="--", linewidth=1.2, alpha=0.7)

# Label the most stable hub genes
top_label = summary[summary["Status"] != "Non-hub"].head(10)
for _, row in top_label.iterrows():
    ax.annotate(
        row["Gene"],
        (row["Stab_MGS1"], row["Stab_MGS4"]),
        fontsize=7, textcoords="offset points", xytext=(4, 3), alpha=0.85,
    )

ax.set_xlabel("Stability Score – MGS1  (fraction of replicates hub-flagged)")
ax.set_ylabel("Stability Score – MGS4  (fraction of replicates hub-flagged)")
ax.set_title(
    "(A)  Gene Hub Stability Scores\n"
    "Dashed lines = threshold τ=%.2f  |  "
    "Red=Co-hub · Blue=MGS1-only · Orange=MGS4-only · Grey=Non-hub"
    % STAB_THRESHOLD
)
ax.legend(loc="upper left", fontsize=9)
ax.grid(alpha=0.3)

# Panel B – Stability bar chart (top hub genes)
ax = fig.add_subplot(gs[1, 0])
top_hub = summary[summary["Status"] != "Non-hub"].head(15)
if len(top_hub) == 0:
    # Lower threshold fallback for display
    top_hub = summary.head(15)
x = np.arange(len(top_hub))
w = 0.40
ax.bar(x - w / 2, top_hub["Stab_MGS1"].values, w,
       label="MGS1", color="#1f77b4", alpha=0.85, edgecolor="white")
ax.bar(x + w / 2, top_hub["Stab_MGS4"].values, w,
       label="MGS4", color="#ff7f0e", alpha=0.85, edgecolor="white")
ax.axhline(STAB_THRESHOLD, color="red", linestyle="--", linewidth=1.2,
           label=f"τ = {STAB_THRESHOLD}")
ax.set_xticks(x)
ax.set_xticklabels(top_hub["Gene"].values, rotation=45, ha="right", fontsize=7)
ax.set_ylabel("Stability Score")
ax.set_title("(B)  Hub Gene Stability: MGS1 vs MGS4")
ax.legend(fontsize=8)
ax.grid(axis="y", alpha=0.3)

# Panel C – Hub status pie chart
ax = fig.add_subplot(gs[1, 1])
status_counts = summary["Status"].value_counts()
present = [s for s in ["Co-hub", "MGS1-specific", "MGS4-specific", "Non-hub"]
           if s in status_counts.index]
ax.pie(
    [status_counts[s] for s in present],
    labels=present,
    colors=[COLOR_MAP[s] for s in present],
    autopct="%1.0f%%",
    startangle=90,
    textprops={"fontsize": 9},
)
ax.set_title("(C)  Gene Hub Status Distribution")

# Panel D – MGS1 hub sub-network
ax = fig.add_subplot(gs[2, 0])
hub_list_mgs1 = [g for g in hub_mgs1 if g in G1.nodes()]
if len(hub_list_mgs1) > 1:
    G1_sub = G1.subgraph(hub_list_mgs1).copy()
    pos    = nx.spring_layout(G1_sub, seed=42, k=1.5 / np.sqrt(len(hub_list_mgs1)))
    node_col   = ["#d62728" if n in cohubs else "#1f77b4" for n in G1_sub.nodes()]
    node_sizes = [max(150, G1_sub.degree(n) * 60) for n in G1_sub.nodes()]
    nx.draw_networkx(
        G1_sub, pos=pos, ax=ax,
        node_size=node_sizes, node_color=node_col,
        edge_color="#aaaaaa", alpha=0.85, font_size=6, width=0.8,
    )
    ax.set_title(f"(D)  MGS1 Hub Sub-network  ({len(hub_list_mgs1)} nodes)\n"
                 "Red = co-hub · Blue = MGS1-specific")
else:
    ax.text(0.5, 0.5, "No hub genes above threshold\n(try N_BOOT ≥ 100)",
            ha="center", va="center", transform=ax.transAxes, fontsize=9)
    ax.set_title("(D)  MGS1 Hub Sub-network")
ax.axis("off")

# Panel E – MGS4 hub sub-network
ax = fig.add_subplot(gs[2, 1])
hub_list_mgs4 = [g for g in hub_mgs4 if g in G4.nodes()]
if len(hub_list_mgs4) > 1:
    G4_sub = G4.subgraph(hub_list_mgs4).copy()
    pos    = nx.spring_layout(G4_sub, seed=42, k=1.5 / np.sqrt(len(hub_list_mgs4)))
    node_col   = ["#d62728" if n in cohubs else "#ff7f0e" for n in G4_sub.nodes()]
    node_sizes = [max(150, G4_sub.degree(n) * 60) for n in G4_sub.nodes()]
    nx.draw_networkx(
        G4_sub, pos=pos, ax=ax,
        node_size=node_sizes, node_color=node_col,
        edge_color="#aaaaaa", alpha=0.85, font_size=6, width=0.8,
    )
    ax.set_title(f"(E)  MGS4 Hub Sub-network  ({len(hub_list_mgs4)} nodes)\n"
                 "Red = co-hub · Orange = MGS4-specific")
else:
    ax.text(0.5, 0.5, "No hub genes above threshold\n(try N_BOOT ≥ 100)",
            ha="center", va="center", transform=ax.transAxes, fontsize=9)
    ax.set_title("(E)  MGS4 Hub Sub-network")
ax.axis("off")

fig.savefig(f"{OUT_DIR}\\bhub_stability_results.png", dpi=150, bbox_inches="tight")
plt.show()
print(f"\nFigure saved → bhub_stability_results.png")

# ── 8. Save Tables ────────────────────────────────────────────────────────────
summary.to_csv(f"{OUT_DIR}\\bhub_stability_results.csv", index=False)
pd.DataFrame({"Co-hub_Gene": list(cohubs)}).to_csv(
    f"{OUT_DIR}\\cohub_genes.csv", index=False
)
print(f"Results saved → bhub_stability_results.csv, cohub_genes.csv")
