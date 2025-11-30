library(RColorBrewer)
library(MEGENA)
library(visNetwork)
library(readr)
library(igraph)
library(dplyr)

big_data = read.csv("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/gene_input.csv", row.names = 1)
cont_late = subset(big_data, mgs_level %in% c('MGS1', "MGS4"))

