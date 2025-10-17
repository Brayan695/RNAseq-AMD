library(RColorBrewer)
library(MEGENA)
library(visNetwork)
library(readr)

# 1. Load data
genes = read.csv("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/aak500_cpmdat.csv")
info = read.delim("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/gene_info.tsv")

genes_control = subset(genes, mgs_level == 'MGS1')

# Remove non-expression column(s)
expr = genes_control[, !(colnames(genes_control) %in% c("mgs_level"))]

# Convert to numeric
expr_num = apply(expr, 2, as.numeric)
expr_num = matrix(expr_num, nrow = nrow(expr), ncol = ncol(expr))
rownames(expr_num) = expr$X
colnames(expr_num) = colnames(expr)

expr_num = expr_num[,-1]

# Transpose
expr_num = t(expr_num)

# 2. Calculate correlation
corr = calculate.correlation(expr_num)

# 3. Construct Planar Filtered Network (PFN)
# Note: If this step fails with an "empty network" error, you can add 'doPerm = FALSE'
pfn = calculate.PFN(corr)

# 4. Create the igraph object
pfn_g = igraph::graph_from_data_frame(d = pfn, directed = FALSE)

# 5. Run MEGENA on igraph
meg = do.MEGENA(
  g = pfn_g,
  mod.pval = 0.05,
  hub.pval = 0.05,
  min.size = 10
)

# 6. Prepare data for visNetwork
# Create nodes and edges data frames from the igraph object
nodes = data.frame(id = V(pfn_g)$name,
                   value = igraph::degree(pfn_g))

edges = pfn
colnames(edges) = c("from", "to")
colnames(edges)[3] = "weight"

# Merge gene symbols and module groups into the nodes data frame
module_summary = MEGENA.ModuleSummary(meg)
modules = module_summary$modules

modules_df = data.frame(id = unlist(modules),
                        group = rep(names(modules), sapply(modules, length)))

nodes = merge(nodes, modules_df, by = "id", all.x = TRUE)
nodes$group[is.na(nodes$group)] = "No Module"

nodes = merge(nodes, info, by.x = "id", by.y = "ensembl_gene_id", all.x = TRUE)
colnames(nodes)[colnames(nodes) == "external_gene_name"] = "label"
nodes$label[is.na(nodes$label)] = nodes$id[is.na(nodes$label)] # Use Ensembl ID as label if gene name is missing

# Add tooltips for more information on hover
nodes$title = paste0("<b>Gene:</b> ", nodes$label,
                     "<br><b>Module:</b> ", nodes$group,
                     "<br><b>Degree:</b> ", nodes$value)

# This keeps the FIRST occurrence of a duplicate ID.
nodes_unique = nodes[!duplicated(nodes$id), ]

# 7. Generate the interactive plot
visNetwork(nodes_unique, edges, main = "MEGENA Network Visualization 500 (MGS1)") %>%
  visIgraphLayout(layout = "layout_with_fr") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visLegend()





