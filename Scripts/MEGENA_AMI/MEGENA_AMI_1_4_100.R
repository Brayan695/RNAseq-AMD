library(RColorBrewer)
library(MEGENA)
library(visNetwork)
library(readr)
library(igraph)

# 1. Load data
genes = read.csv("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/aak100_cpmdat.csv", row.names = 1)
info = read.delim("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/gene_info.tsv")
distance = read.csv("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/AMI/ami_edges.csv")

# 2. Construct Planar Filtered Network (PFN)
pfn = calculate.PFN(distance)

# 3. Create igraph object
pfn_g = igraph::graph_from_data_frame(d = pfn, directed = FALSE)

# 4. Run MEGENA
meg = do.MEGENA(
  g = pfn_g,
  mod.pval = 0.05,
  hub.pval = 0.05,
  min.size = 10,
  n.perm = 100
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
visNetwork(nodes_unique, edges, main = "MEGENA AMI 100 Class (MGS1 + MGS4)") %>%
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

# ============================================================
# Eigengene Extraction and Association with Phenotype (MEGENA)
# ============================================================

library(dplyr)
library(ggplot2)
library(pheatmap)

# ---- 1. Load your original data ----
genes = read.csv("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/aak100_cpmdat.csv",
                 check.names = FALSE, stringsAsFactors = FALSE)

# The samples are rows, genes are columns, and "mgs_level" is the phenotype
class_labels = genes$mgs_level
expr = genes[, !(colnames(genes) %in% c("mgs_level"))]

# Ensure numeric and transpose (samples × genes)
expr_mat = as.matrix(expr)
expr_num = suppressWarnings(apply(expr_mat, 2, as.numeric))
rownames(expr_num) = genes$X  # sample names if present
expr_num = as.data.frame(expr_num)

# ---- 2. Load MEGENA results ----
modules = module_summary$modules

# ---- 3. Compute eigengenes manually (PC1 per module) ----
eigengenes = list()
for (m in names(modules)) {
  genes_in_mod = intersect(modules[[m]], colnames(expr_num))
  if (length(genes_in_mod) < 2) next
  sub_expr = expr_num[, genes_in_mod, drop = FALSE]
  pca = prcomp(sub_expr, scale. = TRUE, center = TRUE)
  eigengenes[[m]] = pca$x[, 1]
}

ME_df = as.data.frame(eigengenes)
rownames(ME_df) = rownames(expr_num)

# ---- 4. Add phenotype labels ----
ME_df$Phenotype = class_labels
ME_df$Phenotype = factor(ME_df$Phenotype)

# ---- 5.1. Wilcoxon tests (paired if possible, else unpaired) ----

# Helper: rank-biserial correlation effect size for unpaired Wilcoxon
rank_biserial_unpaired <- function(x, g) {
  g <- droplevels(factor(g))
  stopifnot(nlevels(g) == 2)
  
  # keep complete cases
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]
  g <- g[ok]
  
  g1 <- levels(g)[1]
  g2 <- levels(g)[2]
  
  x1 <- x[g == g1]
  x2 <- x[g == g2]
  n1 <- length(x1); n2 <- length(x2)
  if (n1 == 0 || n2 == 0) return(NA_real_)
  
  # pooled ranks (average ranks for ties)
  r <- rank(x, ties.method = "average")
  
  # rank sum for group1
  W1 <- sum(r[g == g1])
  
  # Mann–Whitney U for group1
  U1 <- W1 - n1 * (n1 + 1) / 2
  
  # Common Language Effect Size (probability of superiority) for group2 over group1:
  # A = P(X2 > X1) + 0.5 P(X2 = X1) = 1 - U1/(n1*n2)
  A <- 1 - (U1 / (n1 * n2))
  
  # Rank-biserial correlation: r_rb = 2A - 1   (guaranteed in [-1, 1])
  r_rb <- 2 * A - 1
  
  # Numerical safety
  r_rb <- max(-1, min(1, r_rb))
  return(as.numeric(r_rb))
}

# Helper: matched-pairs rank-biserial correlation for paired signed-rank
rank_biserial_paired <- function(d) {
  d <- d[!is.na(d)]
  if (length(d) == 0) return(NA_real_)
  # r_rb = (W+ - W-) / (W+ + W-)
  rp <- rank(abs(d), ties.method = "average")
  Wpos <- sum(rp[d > 0])
  Wneg <- sum(rp[d < 0])
  if ((Wpos + Wneg) == 0) return(0)
  (Wpos - Wneg) / (Wpos + Wneg)
}

# ---- OPTIONAL: define pairing if you have it ----
# If your input data has a donor/pairing column, set it here.
# Examples:
# genes$donor_id  (common)
# genes$pair_id
# If none exists, we will treat groups as independent (rank-sum).
# pair_col <- NULL
# if ("donor_id" %in% colnames(genes)) pair_col <- "donor_id"
# if ("pair_id"  %in% colnames(genes)) pair_col <- "pair_id"
# 
# if (!is.null(pair_col)) {
#   ME_df$PairID <- genes[[pair_col]]
# } else {
#   ME_df$PairID <- NA
# }

# Determine if we can run a *valid* paired test:
# - exactly 2 phenotype levels
# - at least one PairID that appears once in each phenotype
can_do_paired <- FALSE
if (!all(is.na(ME_df$PairID))) {
  tab <- table(ME_df$PairID, ME_df$Phenotype)
  # keep PairIDs that have >=1 in BOTH phenotypes
  keep_pairs <- rownames(tab)[apply(tab > 0, 1, all)]
  if (length(keep_pairs) >= 3) {  # require at least a few pairs
    can_do_paired <- TRUE
  }
}

test_type_used <- if (can_do_paired) "paired_signed_rank" else "unpaired_rank_sum"
message("Wilcoxon test type: ", test_type_used)

module_stats <- lapply(names(eigengenes), function(mod) {
  x <- ME_df[[mod]]
  g <- droplevels(ME_df$Phenotype)
  
  # group summaries
  control_vals <- x[g == "MGS1"]
  late_vals    <- x[g == "MGS4"]
  
  Control_mean   <- mean(control_vals, na.rm = TRUE)
  Late_mean      <- mean(late_vals, na.rm = TRUE)
  Control_median <- median(control_vals, na.rm = TRUE)
  Late_median    <- median(late_vals, na.rm = TRUE)
  
  if (can_do_paired) {
    # build paired vectors by PairID intersection
    tab <- table(ME_df$PairID, ME_df$Phenotype)
    keep_pairs <- rownames(tab)[apply(tab > 0, 1, all)]
    sub <- ME_df[ME_df$PairID %in% keep_pairs, c("PairID", "Phenotype", mod)]
    # one value per PairID per phenotype (if duplicates exist, take mean)
    sub_agg <- sub %>%
      group_by(PairID, Phenotype) %>%
      summarise(val = mean(.data[[mod]], na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = Phenotype, values_from = val)
    
    d <- sub_agg$MGS4 - sub_agg$MGS1
    wt <- suppressWarnings(wilcox.test(sub_agg$MGS4, sub_agg$MGS1, paired = TRUE, exact = FALSE))
    
    data.frame(
      Module = mod,
      Test = "Wilcoxon signed-rank (paired)",
      N_pairs = nrow(sub_agg),
      Control_mean = Control_mean,
      Late_mean = Late_mean,
      Control_median = Control_median,
      Late_median = Late_median,
      Median_diff_Late_minus_Control = median(d, na.rm = TRUE),
      W_statistic = as.numeric(wt$statistic),
      p_value = wt$p.value,
      Effect_rank_biserial = rank_biserial_paired(d)
    )
  } else {
    wt <- suppressWarnings(wilcox.test(x ~ g, exact = FALSE))  # rank-sum
    data.frame(
      Module = mod,
      Test = "Wilcoxon rank-sum (unpaired)",
      N_MGS1 = sum(g == "MGS1" & !is.na(x)),
      N_MGS4 = sum(g == "MGS4" & !is.na(x)),
      Control_mean = Control_mean,
      Late_mean = Late_mean,
      Control_median = Control_median,
      Late_median = Late_median,
      Median_diff_Late_minus_Control = Late_median - Control_median,
      W_statistic = as.numeric(wt$statistic),
      p_value = wt$p.value,
      Effect_rank_biserial = rank_biserial_unpaired(x, g)
    )
  }
})

module_results <- do.call(rbind, module_stats)
module_results$FDR <- p.adjust(module_results$p_value, method = "BH")

module_results <- module_results[order(module_results$p_value), ]
print(module_results)

# ---- Optional: quick plot for the top hit ----
# top_mod <- module_results$Module[1]
# ggplot(ME_df, aes(x = Phenotype, y = .data[[top_mod]])) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(width = 0.15, height = 0) +
#   labs(title = paste0("Eigengene: ", top_mod, " (", module_results$Test[1], ")"),
#        y = "Module Eigengene (PC1)") +
#   theme_classic()


# ---- 5.2. Compare eigengenes between groups ----
module_stats = lapply(names(eigengenes), function(mod) {
  t_res = t.test(ME_df[[mod]] ~ ME_df$Phenotype)
  data.frame(
    Module = mod,
    Control_mean = mean(ME_df[[mod]][ME_df$Phenotype == "MGS1"], na.rm = TRUE),
    Late_mean = mean(ME_df[[mod]][ME_df$Phenotype == "MGS4"], na.rm = TRUE),
    t_statistic = t_res$statistic,
    p_value = t_res$p.value
  )
})
module_results = do.call(rbind, module_stats)
module_results$FDR = p.adjust(module_results$p_value, method = "BH")

print(module_results[order(module_results$p_value), ])

# ---- 6. Correlation-style heatmap ----
# Convert phenotype to numeric (0 = Control, 1 = Late)
# pheno_num = ifelse(ME_df$Phenotype == "MGS4", 1, 0)
# cor_mat = sapply(ME_df[, names(eigengenes)], function(x) cor(x, pheno_num, method = "pearson"))
# p_mat = sapply(ME_df[, names(eigengenes)], function(x) cor.test(x, pheno_num)$p.value)
# 
# pheatmap(matrix(cor_mat, nrow = 1),
#          color = colorRampPalette(c("blue", "white", "red"))(100),
#          display_numbers = matrix(round(p_mat, 3), nrow = 1),
#          main = "Module–Phenotype Correlation (MGS1 vs MGS4)",
#          cluster_rows = FALSE, cluster_cols = FALSE)
# 
# # ---- 7. Boxplot per module ----
# ME_long = reshape2::melt(ME_df, id.vars = "Phenotype")
# ggplot(ME_long, aes(x = Phenotype, y = value, fill = Phenotype)) +
#   geom_boxplot(outlier.shape = NA, alpha = 0.6) +
#   facet_wrap(~ variable, scales = "free_y") +
#   theme_bw(base_size = 12) +
#   labs(title = "Module Eigengene Expression by MGS Level",
#        x = "Group", y = "Module Eigengene (PC1)") +
#   scale_fill_brewer(palette = "Set2")
