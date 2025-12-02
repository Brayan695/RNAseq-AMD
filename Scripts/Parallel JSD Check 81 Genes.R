library(RColorBrewer)
library(MEGENA)
library(visNetwork)
library(readr)
library(igraph)
library(dplyr)
library(doParallel)
library(foreach)

# ------------------------------------------------------------
# 1. LOAD DATA
# ------------------------------------------------------------
big_data = read.csv(
  "C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/gene_input.csv",
  row.names = 1,
  check.names = FALSE
)

cont_late = subset(big_data, mgs_level %in% c("MGS1", "MGS4"))
cont_late$sample_id = NULL

small_data = read.csv(
  "C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/aak100_cpmdat.csv",
  row.names = 1,
  check.names = FALSE
)

# baseline ENSG IDs from small_data
baseline_ensg = colnames(small_data)[-82]


# ------------------------------------------------------------
# 2. LOAD GENE SYMBOL MAPPING
# ------------------------------------------------------------
gene_info = read.delim(
  "C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/gene_info.tsv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

lookup = gene_info[, c("ensembl_gene_id", "external_gene_name")]

# get their symbols
baseline_symbols = lookup[match(baseline_ensg, lookup$ensembl_gene_id),
                           "external_gene_name"]

# remove NA symbols
baseline_symbols = baseline_symbols[!is.na(baseline_symbols)]

# SYMBOL → ENSG map (needed to convert for sampling)
symbol_to_ensg = lookup$ensembl_gene_id
names(symbol_to_ensg) = lookup$external_gene_name

# ENSG → SYMBOL map (needed to convert MEGENA hubs)
ensg_to_symbol = lookup$external_gene_name
names(ensg_to_symbol) = lookup$ensembl_gene_id

# core genes *in SYMBOLS*
core_symbols = baseline_symbols



# ------------------------------------------------------------
# 3. JSD + MEGENA WRAPPER
# ------------------------------------------------------------
JSD_megena_run = function(expr_mat) {
  
  expr_num = suppressWarnings(apply(expr_mat, 2, as.numeric))
  expr_num = t(expr_num)
  
  rownames(expr_num) = colnames(expr_mat) # genes
  colnames(expr_num) = rownames(expr_mat) # samples
  
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
  
  gene_ids = rownames(expr_num)
  n = length(gene_ids)
  JSDmat = matrix(0, n, n, dimnames = list(gene_ids, gene_ids))
  
  for (i in 1:(n - 1)) {
    p = expr_num[i, ]
    for (j in (i + 1):n) {
      q = expr_num[j, ]
      JSDmat[i, j] = JSDmat[j, i] = JSD_div2(p, q)
    }
  }
  
  D = sqrt(JSDmat)
  D = D / max(D)
  S = 1 - D
  
  edges_df = as.data.frame(as.table(S))
  colnames(edges_df) = c("i", "j", "w")
  edges_df = subset(edges_df, i != j)
  
  pfn = MEGENA::calculate.PFN(edges_df)
  pfn = pfn %>% group_by(row, col) %>%
    summarise(weight = mean(weight), .groups = "drop")
  
  g = graph_from_data_frame(pfn, directed = FALSE)
  g = simplify(g, remove.multiple = TRUE,
                edge.attr.comb = list(weight = "mean", "first"))
  
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
  if (!"hub.output" %in% names(meg)) return(character(0))
  
  hubs = meg$hub.output$hub.list
  
  if (is.null(hubs) || is.atomic(hubs)) return(character(0))
  
  # IMPORTANT: remove repeated hub names within the same run
  unique(as.character(unlist(hubs)))
}



# ------------------------------------------------------------
# 4. PARALLEL SETUP
# ------------------------------------------------------------
num_cores = parallel::detectCores() - 1
cl = makeCluster(num_cores)
registerDoParallel(cl)
cat("Running on", num_cores, "cores.\n")


# ------------------------------------------------------------
# 5. SAMPLING — baseline SYMBOLS always included
# ------------------------------------------------------------

gene_pool = colnames(cont_late)[colnames(cont_late) != "mgs_level"]

# Convert core SYMBOLS → ENSG (only keep those actually present)
core_ensg = symbol_to_ensg[core_symbols]          # convert symbols to ENSG
core_ensg = core_ensg[core_ensg %in% gene_pool]   # keep valid ENSG IDs

target_size = 500

all_hubs = foreach(
  k = 1:100,
  .combine = "c",
  .packages = c("MEGENA", "igraph", "dplyr")
) %dopar% {
  
  # forced baseline genes (now ENSG)
  forced_genes = core_ensg
  
  # remaining ENSG IDs not in the core set
  remaining_pool = setdiff(gene_pool, forced_genes)
  
  # number of additional genes needed
  n_needed = max(0, target_size - length(forced_genes))
  
  # random fill
  remaining_genes = if (n_needed > 0) {
    sample(remaining_pool, n_needed)
  } else character(0)
  
  sampled_genes = c(forced_genes, remaining_genes)
  
  # extract expression matrix columns
  sampled_expr = cont_late[, sampled_genes, drop = FALSE]
  
  # run MEGENA
  JSD_megena_run(sampled_expr)
}

stopCluster(cl)

all_hubs = as.character(all_hubs)


# ------------------------------------------------------------
# 6. FREQUENCY TABLE — FIXED (<=100 always)
# ------------------------------------------------------------
hub_frequency = sort(table(all_hubs), decreasing = TRUE)

hub_freq_df = as.data.frame(hub_frequency)
colnames(hub_freq_df) = c("gene_ENSG", "frequency")


# ------------------------------------------------------------
# 7. SAVE
# ------------------------------------------------------------
output_file = "C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/JSD_hub_frequency2.csv"
write.csv(hub_freq_df, output_file, row.names = FALSE)

cat("Hub frequency table saved to:\n", output_file, "\n")
