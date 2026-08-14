"""
AMD Patient Similarity Network -- Node2Vec embedding + spectral clustering.

Follows the same method as:
  - Patock, Ratnapriya & Barman, "A Graphical Method for Identifying Gene
    Clusters from RNA Sequencing Data" (arXiv:2511.09590)
  - https://github.com/arkobarman/graphical_gene_clustering_RNASeq
    (utils.py: calculate_coexpression_graph -> graph_to_node2vec_embeddings
     -> plot_dendrogram_with_embeddings)

The repo's pipeline builds a GENE co-expression network (nodes = genes,
edges = thresholded CS-CORE co-expression), embeds genes with node2vec, and
spectral-clusters the embeddings. Here the same three-stage pipeline is
transposed so that nodes = PATIENTS instead of genes:

  CS-CORE gene-gene co-expression  ->  Pearson correlation between patients'
                                        expression profiles across the AMD
                                        gene panel (CS-CORE's statistical
                                        model does not transpose to patients,
                                        so plain correlation is used instead,
                                        thresholded exactly as in the repo)
  calculate_coexpression_graph      ->  calculate_patient_similarity_graph
  graph_to_node2vec_embeddings      ->  graph_to_node2vec_embeddings (same)
  plot_dendrogram_with_embeddings   ->  spectral_cluster_and_plot (same
                                         SpectralClustering(affinity='rbf')
                                         + UMAP visualization, plus a second
                                         plot colored by mgs_level as a
                                         held-out validation check)

Inputs (relative to this script's parent folder):
  ../Dataset/aak100_cpmdat.csv      166 patients x 81 AMD genes, CPM
  ../Dataset/MetaSheet_1_4_17.csv   166 patients x clinical covariates
                                     (mgs_level used only for validation,
                                     never as a feature)

Outputs (written next to this script):
  patient_correlation_matrix.csv
  patient_similarity_graph.graphml
  patient_node2vec_embeddings.csv
  patient_cluster_labels.csv
  patient_umap_clusters.png
  patient_umap_mgs_level.png
"""
import pathlib

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from node2vec import Node2Vec
from sklearn.cluster import SpectralClustering
from umap import UMAP

HERE = pathlib.Path(__file__).resolve().parent
DATASET_DIR = HERE.parent / "Dataset"

RANDOM_STATE = 0

# --- hyperparameters (mirrors the repo's pipeline.ipynb hyperparameter cells) ---
THRESHOLD = 0.24     # tau: edge if |correlation| > tau (same role as repo's threshold=0.2)
DIMENSIONS = 32       # node2vec embedding dimensionality
WALK_LENGTH = 10
NUM_WALKS = 200
P = 1
Q = 1
WINDOW = 10
N_CLUSTERS = 4        # matches the 4 MGS severity levels for interpretability


# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
def load_data():
    expr = pd.read_csv(DATASET_DIR / "aak100_cpmdat.csv")
    meta = pd.read_csv(DATASET_DIR / "MetaSheet_1_4_17.csv")

    expr = expr.sort_values("sample_id").reset_index(drop=True)
    meta = meta.sort_values("sample_id").reset_index(drop=True)
    assert (expr["sample_id"].values == meta["sample_id"].values).all(), "sample_id mismatch"

    gene_cols = [c for c in expr.columns if c.startswith("ENSG")]
    patient_ids = expr["sample_id"].astype(str).values
    return expr, meta, gene_cols, patient_ids


# ---------------------------------------------------------------------------
# 2. Patient-patient correlation matrix
#    (log2(CPM+1), per-gene z-score, then Pearson correlation between
#    patients across the gene panel -- the patient-axis analogue of the
#    repo's calc_CSCORE step)
# ---------------------------------------------------------------------------
def compute_patient_correlation(expr, gene_cols):
    X = np.log2(expr[gene_cols].values + 1)
    X = (X - X.mean(axis=0)) / X.std(axis=0)
    corr = np.corrcoef(X)
    return corr


# ---------------------------------------------------------------------------
# 3. Threshold graph construction
#    (patient analogue of utils.calculate_coexpression_graph: undirected
#    edge added wherever |correlation| exceeds tau)
# ---------------------------------------------------------------------------
def calculate_patient_similarity_graph(corr_matrix, patient_ids, threshold=THRESHOLD):
    n = corr_matrix.shape[0]
    G = nx.Graph()
    G.add_nodes_from(patient_ids)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(corr_matrix[i, j]) > threshold:
                G.add_edge(patient_ids[i], patient_ids[j], weight=float(corr_matrix[i, j]))
    return G


# ---------------------------------------------------------------------------
# 4. Node2Vec embeddings (identical to utils.graph_to_node2vec_embeddings)
# ---------------------------------------------------------------------------
def graph_to_node2vec_embeddings(G, dimensions=DIMENSIONS, walk_length=WALK_LENGTH,
                                  num_walks=NUM_WALKS, p=P, q=Q, workers=4,
                                  window=WINDOW, min_count=1, batch_words=4):
    node2vec = Node2Vec(G, dimensions=dimensions, walk_length=walk_length,
                         num_walks=num_walks, p=p, q=q, workers=workers, seed=RANDOM_STATE)
    model = node2vec.fit(window=window, min_count=min_count, batch_words=batch_words)
    embeddings = np.array([model.wv[str(node)] for node in G.nodes()])
    return embeddings


# ---------------------------------------------------------------------------
# 5. Spectral clustering + UMAP visualization
#    (patient analogue of utils.plot_dendrogram_with_embeddings)
# ---------------------------------------------------------------------------
def spectral_cluster_and_plot(embeddings, patient_ids, mgs_levels, n_clusters=N_CLUSTERS,
                               random_state=RANDOM_STATE):
    spectral = SpectralClustering(n_clusters=n_clusters, random_state=random_state, affinity="rbf")
    labels = spectral.fit_predict(embeddings)

    umap_model = UMAP(n_components=2, random_state=random_state)
    x_reduced = umap_model.fit_transform(embeddings)

    # plot 1: colored by cluster assignment
    plt.figure(figsize=(9, 8))
    scatter = plt.scatter(x_reduced[:, 0], x_reduced[:, 1], c=labels, s=45, cmap="viridis")
    plt.title(f"Patient Similarity Network -- Spectral Clustering (k={n_clusters})\n"
              f"embedding dim={embeddings.shape[1]}, p={P}, q={Q}, tau={THRESHOLD}")
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.colorbar(scatter, label="cluster")
    plt.tight_layout()
    plt.savefig(HERE / "patient_umap_clusters.png", dpi=200)
    plt.close()

    # plot 2: colored by mgs_level (held-out validation label, NOT used to build the graph)
    plt.figure(figsize=(9, 8))
    colors = ["#3f7cd6" if m == 1 else "#d64545" for m in mgs_levels]
    plt.scatter(x_reduced[:, 0], x_reduced[:, 1], c=colors, s=45, edgecolors="white", linewidths=0.4)
    plt.title("Patient Similarity Network -- UMAP colored by mgs_level\n"
              "(blue=1 mild, red=4 advanced; label NOT used to build the graph)")
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.tight_layout()
    plt.savefig(HERE / "patient_umap_mgs_level.png", dpi=200)
    plt.close()

    return labels


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
def main():
    expr, meta, gene_cols, patient_ids = load_data()
    print(f"Loaded {len(patient_ids)} patients, {len(gene_cols)} genes")

    corr = compute_patient_correlation(expr, gene_cols)
    corr_df = pd.DataFrame(corr, index=patient_ids, columns=patient_ids)
    corr_df.to_csv(HERE / "patient_correlation_matrix.csv")

    G = calculate_patient_similarity_graph(corr, patient_ids, threshold=THRESHOLD)
    n_edges = G.number_of_edges()
    n_isolated = sum(1 for _, d in G.degree() if d == 0)
    n_components = nx.number_connected_components(G)
    print(f"Graph: {n_edges} edges, {n_components} connected component(s), "
          f"{n_isolated} isolated node(s), mean degree {2 * n_edges / len(patient_ids):.2f}")
    nx.write_graphml(G, HERE / "patient_similarity_graph.graphml")

    embeddings = graph_to_node2vec_embeddings(G)
    emb_df = pd.DataFrame(embeddings, index=list(G.nodes()))
    emb_df.to_csv(HERE / "patient_node2vec_embeddings.csv")

    mgs_by_id = dict(zip(meta["sample_id"].astype(str), meta["mgs_level"]))
    mgs_levels = [mgs_by_id[n] for n in G.nodes()]

    labels = spectral_cluster_and_plot(embeddings, list(G.nodes()), mgs_levels)

    result_df = pd.DataFrame({
        "sample_id": list(G.nodes()),
        "cluster": labels,
        "mgs_level": mgs_levels,
    })
    result_df.to_csv(HERE / "patient_cluster_labels.csv", index=False)

    print("\nCluster sizes:")
    print(result_df["cluster"].value_counts().sort_index())
    print("\nCluster x mgs_level cross-tab:")
    print(pd.crosstab(result_df["cluster"], result_df["mgs_level"]))
    print("\nWrote outputs to", HERE)


if __name__ == "__main__":
    main()
