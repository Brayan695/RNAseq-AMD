library(RColorBrewer)
library(MEGENA)
library(visNetwork)
library(readr)
library(igraph)
library(dplyr)


# Load and filter data

big_data = read.csv("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/gene_input.csv",
                    row.names = 1, check.names = FALSE)
cont_late = subset(big_data, mgs_level %in% c("MGS1", "MGS4"))


# Divergence + MEGENA function

JSD_megena_run <- function(expr_mat) {
  
  # Prepare numeric matrix
  expr_num = suppressWarnings(apply(expr_mat, 2, as.numeric))
  expr_num = t(expr_num)
  rownames(expr_num) = colnames(expr_mat)
  colnames(expr_num) = rownames(expr_mat)
  
  # ============================
  # Divergence Functions
  # ============================
  KL_div2 = function(p, q) {
    p = p / sum(p); q = q / sum(q)
    eps = 1e-12
    p[p <= 0] = eps; q[q <= 0] = eps
    sum(p * log2(p / q))
  }
  JSD_div2 = function(p, q) {
    p = p / sum(p); q = q / sum(q)
    eps = 1e-12
    p[p <= 0] = eps; q[q <= 0] = eps
    m = 0.5 * (p + q)
    0.5 * KL_div2(p, m) + 0.5 * KL_div2(q, m)
  }
  
  # ============================
  # Compute Pairwise Divergence
  # ============================
  genes_ids = rownames(expr_num)
  n = length(genes_ids)
  
  JSDmat = matrix(0, n, n, dimnames = list(genes_ids, genes_ids))
  
  for (i in 1:(n - 1)) {
    p = expr_num[i, ]
    for (j in (i + 1):n) {
      q = expr_num[j, ]
      jsd = JSD_div2(p, q)
      JSDmat[i, j] = JSDmat[j, i] = jsd
    }
  }
  
  # Convert to similarities
  D = sqrt(JSDmat)
  D = D / max(D)
  S = 1 - D
  
  # ============================
  # PFN network
  # ============================
  edges_df = as.data.frame(as.table(S))
  colnames(edges_df) = c("i", "j", "w")
  edges_df = subset(edges_df, i != j)
  
  pfn = MEGENA::calculate.PFN(edges_df)
  pfn = pfn %>%
    group_by(row, col) %>%
    summarise(weight = mean(weight), .groups = "drop")
  
  g = graph_from_data_frame(pfn, directed = FALSE)
  g = simplify(g, remove.multiple = TRUE,
               edge.attr.comb = list(weight = "mean", "first"))
  
  # ============================
  # MEGENA
  # ============================
  meg = MEGENA::do.MEGENA(
    g = g,
    mod.pval = 0.05,
    hub.pval = 0.05,
    min.size = 10,
    n.perm = 100
  )
  
  # Return hub gene names
  hubs = meg$hub.output$module.hubs$hub
  return(hubs)
}


# RUN 100 EXPERIMENTS

set.seed(123)

all_hubs = c()  # store every hub gene discovered

for (k in 1:100) {
  message("Running experiment ", k, " of 100...")
  
  # sample 100 random genes (excluding mgs_level)
  gene_pool = colnames(cont_late)[colnames(cont_late) != "mgs_level"]
  sampled_genes = sample(gene_pool, 100)
  
  # subset expression
  sampled_expr = cont_late[, sampled_genes, drop = FALSE]
  
  # run pipeline
  hubs = JSD_megena_run(sampled_expr)
  
  # accumulate hub names
  all_hubs = c(all_hubs, hubs)
}


# FINAL: Hub gene frequency

hub_frequency = sort(table(all_hubs), decreasing = TRUE)

hub_frequency

