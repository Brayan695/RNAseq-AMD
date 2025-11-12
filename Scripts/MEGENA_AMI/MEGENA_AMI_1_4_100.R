library(RColorBrewer)
library(MEGENA)
library(visNetwork)
library(readr)
library(igraph)

# 1. Load data
genes = read.csv("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/aak100_cpmdat.csv", row.names = 1)
info = read.delim("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/gene_info.tsv")
distance = read.csv("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/AMI/ami_edges_QNglobal.csv")

# 2. Construct Planar Filtered Network (PFN)
pfn = calculate.PFN(distance)

# 3. Create igraph object
pfn_g = igraph::graph_from_data_frame(d = pfn, directed = FALSE)

# 4. Run MEGENA
meg = do.MEGENA(
  g = pfn_g,
  mod.pval = 0.05,
  hub.pval = 0.05,
  min.size = 10
)

# 5. Prepare data for visNetwork
nodes = data.frame(id = V(pfn_g)$name, value = igraph::degree(pfn_g))
edges = pfn
colnames(edges) = c("from", "to", "weight")

module_summary = MEGENA.ModuleSummary(meg)
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

# 6. Interactive plot
visNetwork(nodes_unique, edges, main = "MEGENA AMI 100 Class (MGS4)") %>%
  visIgraphLayout(layout = "layout_with_fr") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visLegend()

# 7. Modularity
membership = rep(NA, vcount(pfn_g))
names(membership) = V(pfn_g)$name
for (m in names(modules)) {
  mod_nodes = modules[[m]]
  membership[mod_nodes] = m
}
membership[is.na(membership)] = "NoModule"

comm_fact = as.factor(membership)
modularity_score = modularity(pfn_g, comm_fact, weights = E(pfn_g)$weight)
print(modularity_score)

# 8. Other Metrics
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

# Replace expr_num with genes (your expression matrix)
avg_cor = sapply(unique(membership), function(m) {
  if (m == "NoModule") return(NA)
  genes_in_mod = names(membership[membership == m])
  genes_in_mod = intersect(genes_in_mod, colnames(genes))
  if (length(genes_in_mod) < 2) return(NA)
  cor_mat = cor(genes[, genes_in_mod], use = "pairwise.complete.obs")
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
