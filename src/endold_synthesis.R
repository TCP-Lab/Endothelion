#!/usr/bin/env -S Rscript --vanilla

# Endothelion Project
#
# Assumptions and format of expected arguments.
#
# Subsets of GOIS are in the same directory of the full list of GOIs

# Packages ---------------------------------------------------------------------

# Input Loading-----------------------------------------------------------------

# General error message
error_msg <- "\nERROR by endold_synthesis.R\n"

# Check if the correct number of arguments is provided from command-line
if (length(commandArgs(trailingOnly = TRUE)) != 5) {
  cat(error_msg,
      "One or more arguments are missing. Usage:\n\n",
      "Rscript endold_synthesis.R <endo_model> <central_tendency> \\\n",
      "                           <GOIs> <subGOIs_prefix> <out_dir>\n\n")
  quit(status = 1)
}

# Extract command-line arguments.
endo_model <- commandArgs(trailingOnly = TRUE)[1]
central_tendency <- commandArgs(trailingOnly = TRUE)[2]
gois_file <- commandArgs(trailingOnly = TRUE)[3]
subGOIs_prefix <- commandArgs(trailingOnly = TRUE)[4]
out_dir <- commandArgs(trailingOnly = TRUE)[5]

# Check if the target directory exists
if (! dir.exists(endo_model)) {
 cat(error_msg,
     " File \'", endo_model, "\' does not exist.\n", sep = "")
 quit(status = 2)
}

# Check the central_tendency statistic
# mean = arithmetic mean
# nwam = sample size (n)-weighted arithmetic mean
# sdwam = standard deviation (sd)-weighted arithmetic mean
# sewam = standard error (se)-weighted arithmetic mean
# median = well... median
if (! central_tendency %in% c("mean", "nwam", "sdwam", "sewam", "median")) {
  cat(error_msg,
      " Unknown central tendency statistic \'", central_tendency,
      "\'\n", sep = "")
  quit(status = 3)
}

# Check if the list of the Genes of Interest (GOIs) exists.
if (! file.exists(gois_file)) {
  cat(error_msg,
      " File \'", gois_file, "\' does not exist.\n", sep = "")
  quit(status = 4)
}

# Set the output folder
if (! dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Data Loading -----------------------------------------------------------------

# Get the list of files ending with '_profileReport.csv'
gois_stats_files <- list.files(path = endo_model,
                               pattern = "_profileReport\\.csv$",
                               full.names = TRUE, recursive = TRUE)

# Read the files into separate data frames with distinctive column names and
# store them in a list
gois_stats_list <- lapply(gois_stats_files, function(x){
  df <- read.csv(x)
  key <- grep("Symbol", colnames(df))
  colnames(df)[-key] <- paste0(colnames(df)[-key], "_", basename(dirname(x)))
  return(df)
})

# Merge the data frames based on the 'Symbol' column
gois_stats <- Reduce(\(x, y) merge(x, y, by = "Symbol", all = TRUE),
                     gois_stats_list)

# Cross-study Synthesis --------------------------------------------------------

# Compute cross-study stats (only 'mean' is implemented in this draft)
index <- grep("^Mean", colnames(gois_stats))
average_expr <- rowMeans(gois_stats[,index], na.rm = TRUE)
sd_expr <- apply(gois_stats[,index], 1, sd, na.rm = TRUE)
sem_expr <- sd_expr/sqrt(length(index))

# Final average expression
average_expression <- data.frame(Symbol = gois_stats$Symbol,
                                 Mean = average_expr,
                                 Std_Dev = sd_expr,
                                 SEM = sem_expr)

# Save full GOI set as CSV
write.csv(average_expression,
          file.path(out_dir, "AverExpress_ALL_GOIs.csv"),
          row.names = FALSE)

# Functional Subsetting --------------------------------------------------------

# Get the list of functional subsets to extract
subGOIs_files <- list.files(path = dirname(gois_file),
                            pattern = paste0("^", subGOIs_prefix, ".*\\.csv$"),
                            full.names = TRUE, recursive = TRUE)

# Read the files into separate data frames and store them in a named list (since
# 'subGOIs_files' is a character vector, its elements will be used to assign
# names to output list)
subGOIs_list <- sapply(subGOIs_files, \(x) read.csv(x, header = FALSE),
                       simplify = TRUE, USE.NAMES = TRUE)

# Extract GOI subsets and save them with distinctive file names
for (file_name in names(subGOIs_list)) {
  
  # Modify input file names to get output names
  file_name |> basename() |>
    sub(subGOIs_prefix, "AverExpress_", x=_) |>
    sub("\\.csv\\..*$", ".csv", x=_) -> file_label
  
  # Subset and save as CSV
  average_expression |> subset(Symbol %in% subGOIs_list[[file_name]]) |>
    write.csv(file = file.path(out_dir, file_label), row.names = FALSE)
}

# END --------------------------------------------------------------------------

cat("\n", basename(endo_model), " is done!\n", sep = "")
