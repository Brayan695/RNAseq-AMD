"""
AMD Patient Similarity Network (PSN) construction.
Nodes = patients. Edges = fused similarity from gene expression + clinical metadata.
Inputs:  aak100_cpmdat.csv (166 patients x 82 genes, CPM)
         MetaSheet_1_4_17.csv (166 patients x clinical/genetic covariates)
Outputs: patient_similarity_adjacency.csv   (166x166 fused similarity matrix)
         patient_similarity_edgelist.csv    (sparse k-NN graph edge list)
         patient_graph.graphml              (for Cytoscape/Gephi/networkx)
         patient_graph_plot.png             (layout colored by mgs_level, for validation only)
         graph_summary.txt                  (basic stats)
"""
import pandas as pd
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import StandardScaler
import networkx as nx
import matplotlib.pyplot as plt

RANDOM_SEED = 42
K_NEIGHBORS = 10          # edges per patient in final sparse graph
EXPR_WEIGHT = 0.6         # fusion weight for expression-based similarity
META_WEIGHT = 0.4         # fusion weight for metadata-based similarity

# ---------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------
expr = pd.read_csv("/mnt/user-data/uploads/aak100_cpmdat.csv")
meta = pd.read_csv("/mnt/user-data/uploads/MetaSheet_1_4_17.csv")

expr = expr.sort_values("sample_id").reset_index(drop=True)
meta = meta.sort_values("sample_id").reset_index(drop=True)
assert (expr["sample_id"].values == meta["sample_id"].values).all(), "sample_id mismatch"

sample_ids = expr["sample_id"].astype(str).values
n_patients = len(sample_ids)
print(f"Loaded {n_patients} patients, {expr.shape[1]-1} genes, {meta.shape[1]-1} metadata fields")

# ---------------------------------------------------------------
# 2. Expression features: log2(CPM+1), z-score per gene
# ---------------------------------------------------------------
gene_cols = [c for c in expr.columns if c.startswith("ENSG")]
X_expr = np.log2(expr[gene_cols].values + 1)
X_expr = StandardScaler().fit_transform(X_expr)

# ---------------------------------------------------------------
# 3. Metadata features: covariates only (mgs_level held out as validation label)
# ---------------------------------------------------------------
label_col = "mgs_level"
covariate_cols = [c for c in meta.columns if c not in ("sample_id", label_col)]
X_meta = meta[covariate_cols].values.astype(float)
X_meta = StandardScaler().fit_transform(X_meta)

print(f"Expression feature matrix: {X_expr.shape}")
print(f"Metadata feature matrix:   {X_meta.shape} (columns: {covariate_cols})")

# ---------------------------------------------------------------
# 4. Similarity per view (cosine), then fuse
# ---------------------------------------------------------------
sim_expr = cosine_similarity(X_expr)
sim_meta = cosine_similarity(X_meta)

# rescale each view's similarity to [0,1] before fusing so one view can't dominate by scale
def rescale_01(S):
    S = S.copy()
    np.fill_diagonal(S, np.nan)
    lo, hi = np.nanmin(S), np.nanmax(S)
    S = (S - lo) / (hi - lo)
    np.fill_diagonal(S, 1.0)
    return S

sim_expr_r = rescale_01(sim_expr)
sim_meta_r = rescale_01(sim_meta)

sim_fused = EXPR_WEIGHT * sim_expr_r + META_WEIGHT * sim_meta_r
sim_fused_df = pd.DataFrame(sim_fused, index=sample_ids, columns=sample_ids)
sim_fused_df.to_csv("/home/claude/amd_psn/patient_similarity_adjacency.csv")

# ---------------------------------------------------------------
# 5. Sparsify via mutual k-NN thresholding
# ---------------------------------------------------------------
def knn_graph(sim_matrix, ids, k):
    n = sim_matrix.shape[0]
    G = nx.Graph()
    G.add_nodes_from(ids)
    S = sim_matrix.copy()
    np.fill_diagonal(S, -np.inf)
    for i in range(n):
        nn_idx = np.argsort(S[i])[::-1][:k]
        for j in nn_idx:
            w = sim_matrix[i, j]
            G.add_edge(ids[i], ids[j], weight=float(w))
    return G

G = knn_graph(sim_fused, sample_ids, K_NEIGHBORS)

# attach mgs_level + a couple other useful fields as node attributes (for coloring/validation, NOT used in similarity)
for i, sid in enumerate(sample_ids):
    G.nodes[sid]["mgs_level"] = int(meta.loc[i, "mgs_level"])
    G.nodes[sid]["sex_bin"] = int(meta.loc[i, "sex_bin"])

nx.write_graphml(G, "/home/claude/amd_psn/patient_graph.graphml")

edge_rows = [(u, v, d["weight"]) for u, v, d in G.edges(data=True)]
edge_df = pd.DataFrame(edge_rows, columns=["patient_1", "patient_2", "weight"]).sort_values("weight", ascending=False)
edge_df.to_csv("/home/claude/amd_psn/patient_similarity_edgelist.csv", index=False)

# ---------------------------------------------------------------
# 6. Validation plot: does the graph separate by mgs_level without using it as input?
# ---------------------------------------------------------------
pos = nx.spring_layout(G, seed=RANDOM_SEED, weight="weight", k=0.25)
colors = ["#d64545" if G.nodes[n]["mgs_level"] == 4 else "#3f7cd6" for n in G.nodes]

plt.figure(figsize=(9, 9))
nx.draw_networkx_edges(G, pos, alpha=0.15, width=0.6)
nx.draw_networkx_nodes(G, pos, node_color=colors, node_size=60, linewidths=0.5, edgecolors="white")
plt.title(f"AMD Patient Similarity Network (n={n_patients})\nred = mgs_level 4 (advanced), blue = mgs_level 1 (mild) — label NOT used to build graph", fontsize=11)
plt.axis("off")
plt.tight_layout()
plt.savefig("/home/claude/amd_psn/patient_graph_plot.png", dpi=200)
plt.close()

# ---------------------------------------------------------------
# 7. Summary stats
# ---------------------------------------------------------------
degrees = [d for _, d in G.degree()]
same_label_edges = sum(1 for u, v in G.edges() if G.nodes[u]["mgs_level"] == G.nodes[v]["mgs_level"])
total_edges = G.number_of_edges()

summary = f"""AMD Patient Similarity Network — Summary
==========================================
Patients (nodes):        {n_patients}
Genes used (expression):  {len(gene_cols)}
Metadata covariates used: {len(covariate_cols)}  -> {covariate_cols}
Held-out validation label: mgs_level (1=mild, n={sum(meta.mgs_level==1)}; 4=advanced, n={sum(meta.mgs_level==4)})

Fusion weights: expression={EXPR_WEIGHT}, metadata={META_WEIGHT}
k for k-NN sparsification: {K_NEIGHBORS}

Graph edges (directed k-NN unioned into undirected): {total_edges}
Mean node degree: {np.mean(degrees):.2f}  (min {min(degrees)}, max {max(degrees)})
Fraction of edges connecting patients with the SAME mgs_level: {same_label_edges/total_edges:.3f}
  (baseline if random: {(105/166)**2 + (61/166)**2:.3f})

Files written:
  patient_similarity_adjacency.csv  - full 166x166 fused similarity matrix
  patient_similarity_edgelist.csv   - sparse k-NN edge list (patient_1, patient_2, weight)
  patient_graph.graphml             - graph object for Cytoscape/Gephi/networkx
  patient_graph_plot.png            - spring-layout visualization colored by mgs_level
"""
with open("/home/claude/amd_psn/graph_summary.txt", "w") as f:
    f.write(summary)
print(summary)
