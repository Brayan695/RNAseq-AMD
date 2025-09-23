library(RColorBrewer)
library(MEGENA)
library(visNetwork)

# 1. Load data
genes = read.csv("C:/Users/braya/OneDrive/Desktop/RNAseq-AMD/Data/aak100_cpmdat.csv")

# Remove non-expression column(s)
expr = genes[, !(colnames(genes) %in% c("mgs_level"))]

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

edges = to.data.frame(pfn_g)
colnames(edges) = c("from", "to")

# Merge gene symbols and module groups into the nodes data frame
module_summary = MEGENA.ModuleSummary(meg)
modules = module_summary$modules

modules_df = data.frame(id = unlist(modules),
                        group = rep(names(modules), sapply(modules, length)))

nodes = merge(nodes, modules_df, by = "id", all.x = TRUE)
nodes$group[is.na(nodes$group)] = "No Module"

nodes = merge(nodes, gene_info, by.x = "id", by.y = "ensembl_gene_id", all.x = TRUE)
colnames(nodes)[colnames(nodes) == "external_gene_name"] = "label"
nodes$label[is.na(nodes$label)] = nodes$id[is.na(nodes$label)] # Use Ensembl ID as label if gene name is missing

# Add tooltips for more information on hover
nodes$title = paste0("<b>Gene:</b> ", nodes$label,
                     "<br><b>Module:</b> ", nodes$group,
                     "<br><b>Degree:</b> ", nodes$value)

# 7. Generate the interactive plot
visNetwork(nodes, edges, main = "MEGENA Network Visualization") %>%
  visIgraphLayout(layout = "layout_with_fr") %>%
  visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
  visLegend()





