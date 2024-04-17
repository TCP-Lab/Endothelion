#!/usr/bin/env -S Rscript --vanilla

# For reproducibility
set.seed(7)

# Parse arguments
args = commandArgs(trailingOnly = TRUE)
target_path <- args[1]
gois_path <- args[2]
proj <- args[3]

# Load files
read.delim(file.path(target_path,
                     paste0(proj, "_CountMatrix_genes_TPM.tsv"))) -> df
read.csv(gois_path, header = FALSE) |> unlist() -> GOIs

# Find GOIs within the count matrix
logic <- df$SYMBOL %in% GOIs
# Add 5000 more random genes
# (statistically, some them will be GOIs already included in the previous step)
logic[sample(1:length(logic), 5e3)] <- TRUE
# Subset the count matrix to get a slim version of it
df <- df[logic,]

file_name <- paste0("test_", proj, "_CountMatrix_genes_TPM.tsv")
cat("----------------------------------------------------\n")
cat("Saving:", file_name,
    "\n        Test count matrix with", dim(df)[1], "rows.\n")
cat("----------------------------------------------------\n\n")

# Save the slim matrix within the target path
write.table(df, file.path(target_path, file_name),
            quote = FALSE, sep = "\t", row.names = FALSE)
