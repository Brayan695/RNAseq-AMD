# ======================================================================
# Gene Co-Expression Network Analysis using WGCNA + MEGENA
# ======================================================================

# 1. Load Required Libraries
library(RColorBrewer)
library(MEGENA)
library(visNetwork)
library(readr)
library(dplyr)
library(WGCNA)
library(DESeq2)
library(doParallel)
library(igraph)

allowWGCNAThreads()

# 2. Load and Prepare Expression Data
cpm_data = read.csv(
  "C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/aak100_cpmdat.csv",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

pheno_data = cpm_data$mgs_level
cpm_data$mgs_level = NULL

datExpr = as.data.frame(t(cpm_data))
datExpr = mutate_all(datExpr, as.numeric)

# 3. Load and Prepare Gene Annotation
gene_annotation = read.delim(
  "C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/gene_info.tsv",
  header = TRUE,
  stringsAsFactors = FALSE
)

ensembl_ids_in_expr = rownames(datExpr)
gene_annotation_filtered = gene_annotation %>%
  filter(ensembl_gene_id %in% ensembl_ids_in_expr) %>%
  select(ensembl_gene_id, external_gene_name)

gene_annotation_filtered = gene_annotation_filtered %>%
  group_by(external_gene_name) %>%
  mutate(external_gene_name = if_else(
    n() > 1,
    paste0(external_gene_name, "_", ensembl_gene_id),
    external_gene_name
  )) %>%
  ungroup()

gene_map = setNames(
  gene_annotation_filtered$external_gene_name,
  gene_annotation_filtered$ensembl_gene_id
)

datExpr_final = datExpr[rownames(datExpr) %in% names(gene_map), ]
rownames(datExpr_final) = gene_map[rownames(datExpr_final)]

cat("\nFinal Expression Matrix (Genes x Samples):", dim(datExpr_final), "\n")

# 4. Quality Control
gsg = goodSamplesGenes(datExpr_final, verbose = 3)

if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0)
    printFlush(paste("Removing genes:", paste(colnames(datExpr_final)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples) > 0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr_final)[!gsg$goodSamples], collapse = ", ")))
  datExpr_final = datExpr_final[gsg$goodSamples, gsg$goodGenes]
}

# 5. Choose Soft-Thresholding Power (β)
powers = c(1:10, seq(12, 20, 2))
sft = pickSoftThreshold(
  data = t(datExpr_final),
  powerVector = powers,
  verbose = 5,
  networkType = "signed"
)

par(mfrow = c(1, 2))
plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (Power)",
  ylab = "Scale-Free Topology Fit (signed R²)",
  type = "n",
  main = "Scale Independence"
)
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.85, col = "red")

plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (Power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = "Mean Connectivity"
)
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = 0.9, col = "red")

softPower = 5

# 6. Compute Adjacency Matrix
Adj_Matrix = adjacency(
  datExpr = t(datExpr_final),
  power = softPower,
  type = "signed",
  corFnc = "cor"
)

Adj_Matrix = as.matrix(Adj_Matrix)
Adj_Matrix[is.na(Adj_Matrix)] = 0
Adj_Matrix = (Adj_Matrix + t(Adj_Matrix)) / 2

threshold = quantile(Adj_Matrix[Adj_Matrix > 0], 0.95)
Adj_Matrix[Adj_Matrix < threshold] = 0

# 7. Run MEGENA Pipeline
edge_list = which(Adj_Matrix > 0, arr.ind = TRUE)
edges = data.frame(
  from = rownames(Adj_Matrix)[edge_list[, 1]],
  to = colnames(Adj_Matrix)[edge_list[, 2]],
  weight = Adj_Matrix[edge_list]
)
edges = edges[edges$from != edges$to, ]
edges = edges[!duplicated(t(apply(edges[, 1:2], 1, sort))), ]
g = igraph::graph_from_data_frame(edges, directed = FALSE)

meg = do.MEGENA(
  g = g,
  mod.pval = 0.05,
  hub.pval = 0.05,
  min.size = 10
)

module_summary = MEGENA.ModuleSummary(meg)
modules = module_summary$modules

modules_df = data.frame(
  id = unlist(modules),
  group = rep(names(modules), sapply(modules, length))
)

nodes = data.frame(
  id = V(g)$name,
  value = degree(g)
)

nodes = merge(nodes, modules_df, by = "id", all.x = TRUE)
nodes$group[is.na(nodes$group)] = "No Module"
nodes$title = paste0("<b>Gene:</b> ", nodes$id,
                     "<br><b>Module:</b> ", nodes$group,
                     "<br><b>Degree:</b> ", nodes$value)
colnames(edges)[3] = "weight"

visNetwork(nodes, edges, main = "MEGENA Network from WGCNA Adjacency") %>%
  visIgraphLayout(layout = "layout_with_fr") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visLegend()

# Modularity
membership = rep(NA, vcount(g))
names(membership) = V(g)$name
for (m in names(modules)) {
  nodes_mod = modules[[m]]
  membership[nodes_mod] = m
}
membership[is.na(membership)] = "NoModule"
comm_fact = as.factor(membership)
modularity_score = modularity(g, comm_fact, weights = E(g)$weight)
print(modularity_score)

# Other Metrics
conductance = function(graph, nodes) {
  edgelist = igraph::as_data_frame(graph, what = "edges")
  external_edges = sum((edgelist$from %in% nodes & !(edgelist$to %in% nodes)) |
                         (edgelist$to %in% nodes & !(edgelist$from %in% nodes)))
  total_degree = sum(igraph::degree(graph, v = nodes))
  if (total_degree == 0) return(NA_real_)
  return(external_edges / total_degree)
}

conductance_values = sapply(unique(membership), function(m) {
  if (m == "NoModule") return(NA_real_)
  nodes_in_mod = names(membership[membership == m])
  conductance(g, nodes_in_mod)
})

density_values = sapply(unique(membership), function(m) {
  if (m == "NoModule") return(NA)
  subg = induced_subgraph(g, vids = names(membership[membership == m]))
  edge_density(subg)
})

transitivity_values = sapply(unique(membership), function(m) {
  if (m == "NoModule") return(NA)
  subg = induced_subgraph(g, vids = names(membership[membership == m]))
  transitivity(subg, type = "global")
})

avg_cor = sapply(unique(membership), function(m) {
  if (m == "NoModule") return(NA)
  genes_in_mod = names(membership[membership == m])
  if (length(genes_in_mod) < 2) return(NA)
  cor_mat = cor(t(datExpr_final[genes_in_mod, ]))
  mean(cor_mat[lower.tri(cor_mat)], na.rm = TRUE)
})

module_metrics = data.frame(
  Module = unique(membership),
  Conductance = conductance_values,
  Density = density_values,
  Transitivity = transitivity_values,
  AvgCorrelation = avg_cor
)

print(module_metrics)
