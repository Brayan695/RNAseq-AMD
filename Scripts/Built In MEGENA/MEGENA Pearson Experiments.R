library(RColorBrewer)
library(MEGENA)
library(visNetwork)
library(readr)
library(igraph)
library(dplyr)
library(doParallel)
library(foreach)

# ============================================================
# 1. LOAD DATA
# ============================================================

cat("Loading gene expression matrix...\n")

genes_df = read.csv(
  "C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/gene_input.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Optional phenotype filter
genes_df = genes_df[genes_df$mgs_level == "MGS4", ]

expr_df = genes_df[, !(colnames(genes_df) %in% c("sample_id", "mgs_level")), drop = FALSE]

expr_df[] = lapply(expr_df, function(x) as.numeric(as.character(x)))

gene_names = colnames(expr_df)
expr_matrix = as.matrix(expr_df)

cat("Loaded expression matrix with", ncol(expr_matrix), "genes and", nrow(expr_matrix), "samples.\n")

# ============================================================
# 2. GENES OF INTEREST
# ============================================================

core_genes = c(
  "S100B", "SYNC", "PSMB8-AS1", "C1S", "BTN3A3", "HCP5", "HLA-DMB", "LCTL",
  "C3", "PTGS1", "HLA-DQA1", "MOXD1", "F5", "CFB", "C7", "CYP21A2", "PSMB8",
  "MR1", "FGF1", "ANKRD35", "B2M", "CXCL12", "C1R", "RARRES3", "HLA-B",
  "HLA-DMA", "PTX3", "CCDC3", "HLA-DRA", "HLA-DPA1", "SERPING1", "PSMB9",
  "ISLR", "PLCG2", "GYG2", "PODN", "RNASET2", "TIMP4", "PABPC1L", "GDPD2",
  "BTBD19", "CTSZ", "P2RY6", "GAB3", "TNFRSF11B", "NNMT", "VSTM4", "FBLN5",
  "LIX1", "CYP4F12", "TMEM98", "SERPINA5", "ENPP7P11", "IGFN1", "COPZ2",
  "LST1", "ACOX2", "HFE", "SLC9A9", "SLC47A2", "SCRG1", "LINC01411", "CUBN",
  "ANO3", "TMSB4X", "APOBEC3C", "C1RL", "IGFBP7", "RILPL2", "TRIM22", "DAAM2",
  "CYTH4", "C2", "MAOB", "CCND1", "POM121L9P", "DLGAP1-AS2", "GALNT15",
  "PLEKHS1", "HLA-DOA", "ARSD"
)

core_genes = core_genes[core_genes %in% gene_names]

cat("Core genes included every run:\n")
print(core_genes)

sampling_pool = setdiff(gene_names, core_genes)

subset_size = 500
runs = 100

pearson_output_folder = "C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/Pearson/pearson_runs_2/"
dir.create(pearson_output_folder, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 3. FUNCTION: COMPUTE PEARSON FOR ONE SUBSET
# ============================================================

compute_pearson_for_subset = function(selected_genes, run_id, expr_matrix, gene_names, output_folder) {
  
  cat("Running Pearson for experiment", run_id, "with", length(selected_genes), "genes...\n")
  
  idx = match(selected_genes, gene_names)
  sub_expr = expr_matrix[, idx, drop = FALSE]
  
  gene_vars = apply(sub_expr, 2, function(x) var(x, na.rm = TRUE))
  keep = !is.na(gene_vars) & gene_vars > 0
  
  filtered_genes = selected_genes[keep]
  filtered_expr = sub_expr[, keep, drop = FALSE]
  
  if (ncol(filtered_expr) < 2) {
    cat("Run", run_id, ": fewer than 2 variable genes after filtering.\n")
    return(NULL)
  }
  
  corr_mat = cor(filtered_expr, use = "pairwise.complete.obs", method = "pearson")
  
  upper_idx = which(upper.tri(corr_mat), arr.ind = TRUE)
  
  edges_df = data.frame(
    from = filtered_genes[upper_idx[, 1]],
    to = filtered_genes[upper_idx[, 2]],
    weight = corr_mat[upper.tri(corr_mat)],
    stringsAsFactors = FALSE
  )
  
  edges_df = edges_df[!is.na(edges_df$weight), ]
  
  outfile = file.path(output_folder, paste0("pearson_edges_run_", run_id, ".csv"))
  write.csv(edges_df, outfile, row.names = FALSE)
  
  cat("Saved", outfile, "\n")
  
  return(outfile)
}

# ============================================================
# 4. MAIN LOOP — RUN 100 PEARSON EXPERIMENTS
# ============================================================

saved_files = vector("list", runs)

for (run_id in seq_len(runs)) {
  
  remaining_needed = subset_size - length(core_genes)
  
  random_genes = sample(sampling_pool, remaining_needed, replace = FALSE)
  selected_genes = c(core_genes, random_genes)
  
  outfile = compute_pearson_for_subset(
    selected_genes = selected_genes,
    run_id = run_id,
    expr_matrix = expr_matrix,
    gene_names = gene_names,
    output_folder = pearson_output_folder
  )
  
  saved_files[[run_id]] = outfile
}

saved_files = unlist(saved_files)

cat("\nAll Pearson experiments completed!\n")
cat("List of saved files:\n")
print(saved_files)

# ============================================================
# 5. RUN MEGENA ON ALL PEARSON EDGE FILES
# ============================================================

pearson_dir = pearson_output_folder
pearson_files = list.files(pearson_dir, pattern = "\\.csv$", full.names = TRUE)

cat("Found", length(pearson_files), "Pearson edge files.\n")

safe_megena = function(edges_df) {
  
  pfn_sub = tryCatch({
    MEGENA::calculate.PFN(edges_df)
  }, error = function(e) NULL)
  
  if (is.null(pfn_sub) || nrow(pfn_sub) == 0) return(character(0))
  
  g = graph_from_data_frame(pfn_sub, directed = FALSE)
  g = simplify(
    g,
    remove.multiple = TRUE,
    remove.loops = TRUE,
    edge.attr.comb = list(weight = "mean", "first")
  )
  
  if (vcount(g) < 15) return(character(0))
  
  meg = tryCatch({
    MEGENA::do.MEGENA(
      g = g,
      mod.pval = 0.05,
      hub.pval = 0.05,
      min.size = 10,
      n.perm = 100
    )
  }, error = function(e) NULL)
  
  if (is.null(meg)) return(character(0))
  if (is.null(meg$hub.output)) return(character(0))
  if (is.null(meg$hub.output$hub.list)) return(character(0))
  
  hubs = meg$hub.output$hub.list
  
  if (is.atomic(hubs)) return(character(0))
  
  unique(unlist(hubs))
}

num_cores = max(1, parallel::detectCores() - 1)
cl = makeCluster(num_cores)
registerDoParallel(cl)

cat("Running MEGENA on", length(pearson_files), "Pearson files using", num_cores, "cores...\n")

all_hubs = foreach(
  f = pearson_files,
  .combine = "c",
  .packages = c("MEGENA", "igraph", "dplyr")
) %dopar% {
  edges_df = read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
  safe_megena(edges_df)
}

stopCluster(cl)

hub_frequency = sort(table(all_hubs), decreasing = TRUE)
print(hub_frequency)

hub_freq_df = as.data.frame(hub_frequency, stringsAsFactors = FALSE)
colnames(hub_freq_df) = c("gene", "frequency")

output_file = "C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/Pearson/Pearson_hub_frequency_4_81.csv"
write.csv(hub_freq_df, output_file, row.names = FALSE)

cat("Hub frequency table saved to:\n", output_file, "\n")