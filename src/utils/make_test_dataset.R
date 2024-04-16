#!/usr/bin/env -S Rscript --vanilla

# For reproducibility
set.seed(7)

args = commandArgs(trailingOnly = TRUE)

target_path <- args[1]
gois_path <- args[2]
proj <- args[3]

read.delim(file.path(target_path, paste0(proj, "_CountMatrix_genes_TPM.tsv"))) -> df
read.csv(gois_path, header = FALSE) |> unlist() -> GOIs

logic <- df$SYMBOL %in% GOIs
logic[sample(1:length(logic),5e3)] <- TRUE
df <- df[logic,]

file_name <- paste0("test_", proj, "_CountMatrix_genes_TPM.tsv")

cat("----------------------------------------------------\n")
cat("Saving: ", file_name,
    "\n         Test count matrix with", dim(df)[1]," rows.\n")
cat("----------------------------------------------------\n\n")

write.table(df, file.path(target_path, file_name),
            quote = FALSE, sep = "\t", row.names = FALSE)
