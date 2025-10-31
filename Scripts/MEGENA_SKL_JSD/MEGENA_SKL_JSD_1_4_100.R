# ===========================================================
# MEGENA with Information-Theoretic Similarity (SKL or √JSD)
# for Mixed Classes
# ===========================================================

library(RColorBrewer)
library(MEGENA)
library(visNetwork)
library(readr)
library(igraph)
library(dplyr)

# 1. Load and Prepare Data
genes = read.csv("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/aak100_cpmdat.csv",
                 check.names = FALSE, stringsAsFactors = FALSE)
info  = read.delim("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/gene_info.tsv",
                   check.names = FALSE, stringsAsFactors = FALSE)

expr = genes[, !(colnames(genes) %in% c("mgs_level"))]

expr_mat = as.matrix(expr)
expr_num = suppressWarnings(apply(expr_mat, 2, as.numeric))
expr_num = t(expr_num)

colnames(expr_num) = rownames(expr)
rownames(expr_num) = colnames(expr)

# 2. Divergence Functions (base-2 logs)
KL_div2 = function(p, q) {
  p = p / sum(p); q = q / sum(q)
  eps = 1e-12
  p[p <= 0] = eps; q[q <= 0] = eps
  sum(p * log2(p / q))
}

SKL_div2 = function(p, q) 0.5 * (KL_div2(p, q) + KL_div2(q, p))

JSD_div2 = function(p, q) {
  p = p / sum(p); q = q / sum(q)
  eps = 1e-12
  p[p <= 0] = eps; q[q <= 0] = eps
  m = 0.5 * (p + q)
  0.5 * KL_div2(p, m) + 0.5 * KL_div2(q, m)
}

# 3. Compute Pairwise Divergences
genes_ids = rownames(expr_num)
n = length(genes_ids)

SKL = matrix(0, n, n, dimnames = list(genes_ids, genes_ids))
JSD = matrix(0, n, n, dimnames = list(genes_ids, genes_ids))

for (i in 1:(n - 1)) {
  p = expr_num[i, ]
  for (j in (i + 1):n) {
    q = expr_num[j, ]
    skl = SKL_div2(p, q)
    jsd = JSD_div2(p, q)
    SKL[i, j] = SKL[j, i] = skl
    JSD[i, j] = JSD[j, i] = jsd
  }
}
diag(SKL) = 0; diag(JSD) = 0

use_metric = "SKL"  # or "JSD"
if (use_metric == "JSD") {
  D = sqrt(JSD)  # √JSD = metric
} else {
  D = SKL
}
D = D / max(D)
S = 1 - D  # similarity

# 4. Create MEGENA Input
edges_df = as.data.frame(as.table(S))
colnames(edges_df) = c("i", "j", "w")
edges_df = subset(edges_df, i != j)

# 5. PFN + MEGENA
pfn = MEGENA::calculate.PFN(edges_df)
pfn = pfn %>%
  dplyr::group_by(row, col) %>%
  dplyr::summarise(weight = mean(weight), .groups = "drop")

pfn_g = igraph::graph_from_data_frame(d = pfn, directed = FALSE)
pfn_g = igraph::simplify(pfn_g, remove.multiple = TRUE,
                         edge.attr.comb = list(weight = "mean", "first"))

meg = MEGENA::do.MEGENA(
  g = pfn_g,
  mod.pval = 0.05,
  hub.pval = 0.05,
  min.size = 10
)

# 6. Visualization
nodes = data.frame(id = V(pfn_g)$name, value = igraph::degree(pfn_g))
edges = pfn
colnames(edges) = c("from", "to", "weight")

module_summary = MEGENA::MEGENA.ModuleSummary(meg)
modules = module_summary$modules
modules_df = data.frame(id = unlist(modules),
                        group = rep(names(modules), sapply(modules, length)))

nodes = merge(nodes, modules_df, by = "id", all.x = TRUE)
nodes$group[is.na(nodes$group)] = "No Module"
nodes = merge(nodes, info, by.x = "id", by.y = "ensembl_gene_id", all.x = TRUE)
colnames(nodes)[colnames(nodes) == "external_gene_name"] = "label"
nodes$label[is.na(nodes$label)] = nodes$id[is.na(nodes$label)]

nodes$title = paste0("<b>Gene:</b> ", nodes$label,
                     "<br><b>Module:</b> ", nodes$group,
                     "<br><b>Degree:</b> ", nodes$value)
nodes_unique = nodes[!duplicated(nodes$id), ]

visNetwork(nodes_unique, edges,
           main = paste("MEGENA Network Visualization (MGS1 + MGS4,", use_metric, ")")) %>%
  visIgraphLayout(layout = "layout_with_fr") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visLegend()

# 7. Modularity and Other Metrics
membership = rep(NA, vcount(pfn_g))
names(membership) = V(pfn_g)$name
for (m in names(modules)) {
  nds = modules[[m]]
  membership[nds] = m
}
membership[is.na(membership)] = "NoModule"
comm_fact = as.factor(membership)

modularity_score = modularity(pfn_g, comm_fact, weights = E(pfn_g)$weight)
print(modularity_score)

conductance = function(graph, nodes) {
  edgelist = igraph::as_data_frame(graph, what = "edges")
  external_edges = sum((edgelist$from %in% nodes & !(edgelist$to %in% nodes)) |
                         (edgelist$to %in% nodes & !(edgelist$from %in% nodes)))
  total_degree = sum(igraph::degree(graph, v = nodes))
  if (total_degree == 0) return(NA_real_)
  external_edges / total_degree
}

conductance_values = sapply(unique(membership), function(m) {
  if (m == "NoModule") return(NA_real_)
  nodes_in_mod = names(membership[membership == m])
  conductance(pfn_g, nodes_in_mod)
})

density_values = sapply(unique(membership), function(m) {
  if (m == "NoModule") return(NA)
  subg = induced_subgraph(pfn_g, vids = names(membership[membership == m]))
  edge_density(subg)
})

transitivity_values = sapply(unique(membership), function(m) {
  if (m == "NoModule") return(NA)
  subg = induced_subgraph(pfn_g, vids = names(membership[membership == m]))
  transitivity(subg, type = "global")
})

avg_cor = sapply(unique(membership), function(m) {
  if (m == "NoModule") return(NA)
  genes_in_mod = names(membership[membership == m])
  if (length(genes_in_mod) < 2) return(NA)
  cor_mat = cor(t(expr_num[genes_in_mod, ]))
  mean(cor_mat[lower.tri(cor_mat)], na.rm = TRUE)
})

data.frame(
  Module = unique(membership),
  Conductance = conductance_values,
  Density = density_values,
  Transitivity = transitivity_values,
  AvgCorrelation = avg_cor
)
