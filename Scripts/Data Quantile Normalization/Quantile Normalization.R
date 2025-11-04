library(RColorBrewer)
library(MEGENA)
library(visNetwork)
library(readr)
library(igraph)
library(dplyr)
library(preprocessCore)
library(qsmooth)

# Control + Late ----

# 1. Load and Prepare Data
genes = read.csv("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/aak100_cpmdat.csv",
                 check.names = FALSE, stringsAsFactors = FALSE, row.names = 1)
info  = read.delim("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/gene_info.tsv",
                   check.names = FALSE, stringsAsFactors = FALSE)

expr = genes[, !(colnames(genes) %in% c("mgs_level"))]

expr_mat = as.matrix(expr)
expr_num = suppressWarnings(apply(expr_mat, 2, as.numeric))
expr_num = t(expr_num)

colnames(expr_num) = rownames(expr)
rownames(expr_num) = colnames(expr)

class_labels = factor(genes$mgs_level)

# QN
quantile_norm = function(expr_num) {
  expr_num = as.matrix(expr_num)
  normalize.quantiles(expr_num)
}

expr_QN_global = quantile_norm(genes[,-82])
rownames(expr_QN_global) = rownames(genes[,-82])
colnames(expr_QN_global) = colnames(genes[,-82])

expr_QN_global = t(expr_QN_global)


# Class-Specific

quantile_norm_by_class = function(expr_num, classes) {
  classes = factor(classes)
  norm_matrix = matrix(NA, nrow(expr_num), ncol(expr_num),
                        dimnames = dimnames(expr_num))
  
  for (cl in levels(classes)) {
    idx = which(classes == cl)
    sub_mat = expr_num[, idx, drop = FALSE]
    norm_matrix[, idx] = normalize.quantiles(as.matrix(sub_mat))
  }
  norm_matrix
}

expr_QN_class = quantile_norm_by_class(expr_num, class_labels)

# QN (QSmooth)

expr_qsmooth = qsmooth(expr_num, class_labels)
expr_QS = qsmoothData(expr_qsmooth)

# Control ----

genes = read.csv("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/aak100_cpmdat.csv",
                 check.names = FALSE, stringsAsFactors = FALSE, row.names = 1)
info  = read.delim("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/gene_info.tsv",
                   check.names = FALSE, stringsAsFactors = FALSE)

genes_control = subset(genes, mgs_level == "MGS1")
expr = genes_control[, !(colnames(genes_control) %in% c("mgs_level"))]

expr_mat = as.matrix(expr)
expr_num = suppressWarnings(apply(expr_mat, 2, as.numeric))
expr_num = t(expr_num)

colnames(expr_num) = rownames(expr)
rownames(expr_num) = colnames(expr)

class_labels = factor(genes_control$mgs_level)

# QN
quantile_norm = function(expr_num) {
  expr_num = as.matrix(expr_num)
  normalize.quantiles(expr_num)
}

expr_QN_global = quantile_norm(genes[,-82])
rownames(expr_QN_global) = rownames(genes[,-82])
colnames(expr_QN_global) = colnames(genes[,-82])

expr_QN_global = t(expr_QN_global)


# Class-Specific

quantile_norm_by_class = function(expr_num, classes) {
  classes = factor(classes)
  norm_matrix = matrix(NA, nrow(expr_num), ncol(expr_num),
                       dimnames = dimnames(expr_num))
  
  for (cl in levels(classes)) {
    idx = which(classes == cl)
    sub_mat = expr_num[, idx, drop = FALSE]
    norm_matrix[, idx] = normalize.quantiles(as.matrix(sub_mat))
  }
  norm_matrix
}

expr_QN_class = quantile_norm_by_class(expr_num, class_labels)


# Late ----

genes = read.csv("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/aak100_cpmdat.csv",
                 check.names = FALSE, stringsAsFactors = FALSE, row.names = 1)
info  = read.delim("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/gene_info.tsv",
                   check.names = FALSE, stringsAsFactors = FALSE)

genes_control = subset(genes, mgs_level == "MGS4")
expr = genes_control[, !(colnames(genes_control) %in% c("mgs_level"))]

expr_mat = as.matrix(expr)
expr_num = suppressWarnings(apply(expr_mat, 2, as.numeric))
expr_num = t(expr_num)

colnames(expr_num) = rownames(expr)
rownames(expr_num) = colnames(expr)

class_labels = factor(genes_control$mgs_level)

# QN
quantile_norm = function(expr_num) {
  expr_num = as.matrix(expr_num)
  normalize.quantiles(expr_num)
}

expr_QN_global = quantile_norm(genes[,-82])
rownames(expr_QN_global) = rownames(genes[,-82])
colnames(expr_QN_global) = colnames(genes[,-82])

expr_QN_global = t(expr_QN_global)


# Class-Specific

quantile_norm_by_class = function(expr_num, classes) {
  classes = factor(classes)
  norm_matrix = matrix(NA, nrow(expr_num), ncol(expr_num),
                       dimnames = dimnames(expr_num))
  
  for (cl in levels(classes)) {
    idx = which(classes == cl)
    sub_mat = expr_num[, idx, drop = FALSE]
    norm_matrix[, idx] = normalize.quantiles(as.matrix(sub_mat))
  }
  norm_matrix
}

expr_QN_class = quantile_norm_by_class(expr_num, class_labels)








































































