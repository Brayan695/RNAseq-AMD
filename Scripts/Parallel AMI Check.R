library(RColorBrewer)
library(MEGENA)
library(visNetwork)
library(readr)
library(igraph)
library(dplyr)
library(doParallel)
library(foreach)

ami_dir = "C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/AMI/ami_runs_HO_4/"
ami_files = list.files(ami_dir, pattern = "\\.csv$", full.names = TRUE)

cat("Found", length(ami_files), "AMI edge files.\n")

safe_megena = function(edges_df) {
  pfn_sub = tryCatch({
    MEGENA::calculate.PFN(edges_df)
  }, error = function(e) NULL)
  
  if (is.null(pfn_sub) || nrow(pfn_sub) == 0) return(character(0))
  
  g = graph_from_data_frame(pfn_sub, directed = FALSE)
  g = simplify(g, remove.multiple = TRUE,
               edge.attr.comb = list(weight = "mean", "first"))
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
  if (is.null(meg$hub.output$hub.list)) return(character(0))
  
  hubs = meg$hub.output$hub.list
  if (is.atomic(hubs)) return(character(0))
  
  unique(unlist(hubs))
}

num_cores = parallel::detectCores() - 1
cl = makeCluster(num_cores)
registerDoParallel(cl)

cat("Running MEGENA on", length(ami_files), "AMI files using", num_cores, "cores...\n")

all_hubs = foreach(f = ami_files, .combine = "c",
                   .packages = c("MEGENA", "igraph", "dplyr")) %dopar% {
                     edges_df = read.csv(f, check.names = FALSE)
                     safe_megena(edges_df)
                   }

stopCluster(cl)

hub_frequency = sort(table(all_hubs), decreasing = TRUE)
hub_frequency

hub_freq_df = as.data.frame(hub_frequency)
colnames(hub_freq_df) = c("gene", "frequency")

output_file = "C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/AMI/AMI_hub_frequency_4.csv"
write.csv(hub_freq_df, output_file, row.names = FALSE)

cat("Hub frequency table saved to:\n", output_file, "\n")
