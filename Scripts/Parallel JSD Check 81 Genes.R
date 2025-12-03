library(RColorBrewer)
library(MEGENA)
library(visNetwork)
library(readr)
library(igraph)
library(dplyr)
library(doParallel)
library(foreach)

big_data = read.csv("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/gene_input.csv",
                    row.names = 1, check.names = FALSE)
cont_late = subset(big_data, mgs_level %in% c("MGS1", 'MGS4'))
cont_late$sample_id = NULL

JSD_megena_run = function(expr_mat) {
  expr_num = suppressWarnings(apply(expr_mat, 2, as.numeric))
  expr_num = t(expr_num)
  rownames(expr_num) = colnames(expr_mat)
  colnames(expr_num) = rownames(expr_mat)
  
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
  
  unique(unlist(hubs))
}

num_cores = parallel::detectCores() - 1
cl = makeCluster(num_cores)
registerDoParallel(cl)

cat("Running on", num_cores, "cores\n")

gene_pool = colnames(cont_late)[colnames(cont_late) != "mgs_level"]

all_hubs = foreach(k = 1:100, .combine = "c",
                   .packages = c("MEGENA", "igraph", "dplyr")) %dopar% {
                     core_genes = c("S100B", "SYNC", "PSMB8-AS1", "C1S", "BTN3A3", "HCP5", "HLA-DMB", "LCTL",
                                    "C3", "PTGS1", "HLA-DQA1", "MOXD1", "F5", "CFB", "C7", "CYP21A2", "PSMB8",
                                    "MR1", "FGF1", "ANKRD35", "B2M", "CXCL12", "C1R", "RARRES3", "HLA-B",
                                    "HLA-DMA", "PTX3", "CCDC3", "HLA-DRA", "HLA-DPA1", "SERPING1", "PSMB9",
                                    "ISLR", "PLCG2", "GYG2", "PODN", "RNASET2", "TIMP4", "PABPC1L", "GDPD2",
                                    "BTBD19", "CTSZ", "P2RY6", "GAB3", "TNFRSF11B", "NNMT", "VSTM4", "FBLN5",
                                    "LIX1", "CYP4F12", "TMEM98", "SERPINA5", "ENPP7P11", "IGFN1", "COPZ2",
                                    "LST1", "ACOX2", "HFE", "SLC9A9", "SLC47A2", "SCRG1", "LINC01411", "CUBN",
                                    "ANO3", "TMSB4X", "APOBEC3C", "C1RL", "IGFBP7", "RILPL2", "TRIM22", "DAAM2",
                                    "CYTH4", "C2", "MAOB", "CCND1", "POM121L9P", "DLGAP1-AS2", "GALNT15",
                                    "PLEKHS1", "HLA-DOA", "ARSD")
                     core_genes = core_genes[core_genes %in% gene_pool]
                     remaining_pool = setdiff(gene_pool, core_genes)
                     remaining_genes = sample(remaining_pool, 500 - length(core_genes))
                     sampled_genes = c(core_genes, remaining_genes)
                     sampled_expr = cont_late[, sampled_genes, drop = FALSE]
                     JSD_megena_run(sampled_expr)
                   }

stopCluster(cl)

hub_frequency = sort(table(all_hubs), decreasing = TRUE)
hub_frequency

hub_freq_df = as.data.frame(hub_frequency)
colnames(hub_freq_df) = c("gene", "frequency")

output_file = "C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/JSD_hub_frequency_1_4_81.csv"
write.csv(hub_freq_df, output_file, row.names = FALSE)

cat("Hub frequency table saved to:\n", output_file, "\n")
