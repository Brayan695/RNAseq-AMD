library(dplyr)

raw_meta = read.csv("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/MetaSheet.csv")
genes = read.csv("C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/aak100_cpmdat.csv")

meta = raw_meta %>% select(sample_id, A69S_rs10490924, 
                           Y402H_rs1061170, age, sex,
                           medical_history, ocular_history,
                           mgs_level)

binarize = function(x, empties = character(0), prefix = "") {
  x[is.na(x)] = ""
  x = trimws(x)
  x = ifelse(tolower(x) %in% tolower(empties), "", x)
  
  # protect commas inside parentheses: "smoker (1pk/day, 20yrs)" stays one item
  x = gsub(",(?=[^()]*\\))", ";", x, perl = TRUE)
  
  item_lists = lapply(strsplit(x, ","), function(v) {
    v = trimws(v)
    v[v != ""]
  })
  
  universe = unique(unlist(item_lists))
  
  # nothing to encode -> 0-column frame with the right number of rows
  if (length(universe) == 0) {
    return(data.frame(matrix(nrow = length(x), ncol = 0)))
  }
  
  mat = matrix(0L, nrow = length(item_lists), ncol = length(universe),
                dimnames = list(NULL, paste0(prefix, universe)))
  for (j in seq_along(universe)) {
    mat[, j] = as.integer(vapply(item_lists,
                                  function(v) universe[j] %in% v,
                                  logical(1)))
  }
  as.data.frame(mat, check.names = FALSE)
}

# ocular_history: blank and "-" -> all zeros
ocular_bin  = binarize(meta$ocular_history,  empties = c("", "-"), prefix = "oc_")
# medical_history: only blank means empty
medical_bin = binarize(meta$medical_history, empties = c(""),      prefix = "mh_")

meta = cbind(meta, ocular_bin, medical_bin)

meta = meta %>% select(-c(medical_history, ocular_history))

meta_1_4 = meta %>% filter(mgs_level %in% c(1,4))

# Checking if samples align
setequal(genes$sample_id, meta_1_4$sample_id)

write.csv(meta_1_4, "C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/MetaSheet_1_4.csv")

write.csv(meta, "C:/Users/Brayan Gutierrez/Desktop/RNAseq-AMD/Dataset/MetaSheet_Processed.csv")






