# ---------------------------------------------------------------------------
# Cross-check gene expression between gene_input.csv and aak100_cpmdat.csv
#
# What this does:
#   1. Filter gene_input rows to MGS levels 1 and 4 (MGS labels come from
#      MetaSheet.csv, which is row-aligned with gene_input).
#   2. Extract the gene columns from the aak100 file (dropping sample_id and
#      the trailing mgs_level label column).
#   3. Confirm the expression values agree across both files.
#
# Note on gene naming: gene_input uses gene SYMBOLS (TSPAN6, ...) while aak100
# uses ENSG IDs. They are not directly comparable by name, so each aak100 gene
# is matched to its gene_input column by value. The match itself is the check.
# ---------------------------------------------------------------------------

# ---- paths -----------------------------------------------------------------
gene_input_path = "C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/gene_input.csv"
aak100_path     = "C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/aak100_cpmdat.csv"
meta_path       = "C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/MetaSheet.csv"
tol             = 1e-3   # tolerance for "identical" values

# ---- read ------------------------------------------------------------------
# check.names = FALSE keeps original column names (symbols / ENSG IDs) intact.
gene_input = read.csv(gene_input_path, header = TRUE,
                       check.names = FALSE, stringsAsFactors = FALSE)
aak100     = read.csv(aak100_path, header = TRUE,
                       check.names = FALSE, stringsAsFactors = FALSE)
# MetaSheet has a non-UTF-8 byte in a column we don't use; read as latin1.
meta       = read.csv(meta_path, header = TRUE, check.names = FALSE,
                       stringsAsFactors = FALSE, fileEncoding = "latin1")

# gene_input has no sample/MGS column; it is assumed row-aligned with MetaSheet.
stopifnot(nrow(gene_input) == nrow(meta))

# ---- subset rows: MGS 1 and 4 (via MetaSheet) ------------------------------
keep         = meta$mgs_level %in% c(1, 4)
gene_input_r = gene_input[keep, , drop = FALSE]          # rows subset
gi_sample_id = as.character(meta$sample_id[keep])

# ---- the 100 file: gene matrix + its sample ids ----------------------------
# layout: sample_id | 81 ENSG gene columns | mgs_level
gene_cols      = setdiff(names(aak100), c("sample_id", "mgs_level"))
aak_sample_id  = as.character(aak100$sample_id)

# ---- align rows by sample_id (put both in the 100 file's order) ------------
common = aak_sample_id[aak_sample_id %in% gi_sample_id]
stopifnot(length(common) == length(aak_sample_id))       # every 100-file sample present

A = data.matrix(aak100[match(common, aak_sample_id), gene_cols, drop = FALSE])  # 166 x 81
G = data.matrix(gene_input_r[match(common, gi_sample_id), , drop = FALSE])      # 166 x ~18057

# ---- map "columns of the 100 file" onto gene_input columns (by value) ------
map_idx     = integer(length(gene_cols))
n_candidate = integer(length(gene_cols))
for (j in seq_along(gene_cols)) {
  col_max_diff   = apply(abs(G - A[, j]), 2, max)
  cand           = which(col_max_diff <= tol)
  n_candidate[j] = length(cand)
  map_idx[j]     = if (length(cand) >= 1) cand[which.min(col_max_diff[cand])] else NA_integer_
}

# sanity on the mapping itself
if (any(n_candidate == 0))
  warning("Genes with NO matching gene_input column: ",
          paste(gene_cols[n_candidate == 0], collapse = ", "))
if (any(n_candidate > 1))
  warning("Genes matching MORE THAN ONE gene_input column (ambiguous): ",
          paste(gene_cols[n_candidate > 1], collapse = ", "))
if (anyDuplicated(map_idx[!is.na(map_idx)]))
  warning("Two 100-file genes mapped to the SAME gene_input column.")

# ---- build the subset of gene_input (rows x the 100-file columns) ----------
# columns reordered + renamed to the 100 file's gene order/IDs
gene_input_subset = G[, map_idx, drop = FALSE]
colnames(gene_input_subset) = gene_cols                  # 166 x 81, matches A's layout

# ---- element-wise comparison against the 100 file --------------------------
stopifnot(dim(gene_input_subset) == dim(A))
diff_mat   = abs(A - gene_input_subset)
max_diff   = max(diff_mat)
n_mismatch = sum(diff_mat > tol)

cat("subset dimensions (samples x genes):", paste(dim(gene_input_subset), collapse = " x "), "\n")
cat("rows aligned by sample_id:", length(common), "/", nrow(A), "\n")
cat("largest absolute difference across all cells:", max_diff, "\n")
cat("cells differing by more than tol:", n_mismatch, "of", length(diff_mat), "\n")

if (n_mismatch == 0) {
  cat("\nRESULT: fully aligned. Every expression value in the 100 file equals",
      "the corresponding cell of gene_input subset to the MGS 1/4 rows and the",
      "100-file genes.\n")
} else {
  cat("\nRESULT: mismatches found. Locations (sample_id, gene):\n")
  bad = which(diff_mat > tol, arr.ind = TRUE)
  print(data.frame(sample_id = common[bad[, "row"]],
                   gene       = gene_cols[bad[, "col"]],
                   aak100     = A[bad],
                   gene_input = gene_input_subset[bad],
                   abs_diff   = diff_mat[bad]))
}