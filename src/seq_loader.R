


# A constructor function for the bioModel S3 class
# 
# Assumptions:
#
# 1 - Both data and metadata have file names starting with the same series
#     (aka experiment or project) ID (both GEO and ENA are fine), followed by an
#     underscore ("_") and either the distinctive string "CountMatrix" or "meta"
#     (ignoring case) if they are counts or metadata, respectively. Additional
#     characters are allowed after this pattern, including either CSV or TSV
#     file extensions. Examples of good file name patterns are:
#
#         GSE205739_CountMatrix_genes_TPM.tsv
#         PRJNA847413_CountMatrix.csv
#         GSE76528_meta.csv
#
# 2 - Metadata files have proper column names as table header; among them only
#     `ena_run` is mandatory, as the only ID that is certainly unique to each
#     file.
#
# 3 - Count tables have proper column names as header, among which a mandatory
#     ID column matching the regex "gene.*id" or "transcript.*id" (usually being
#     "gene_id" or "transcript_id", respectively).
#
# 4 - In count table, the names of the columns containing the counts of each
#     different run feature the corresponding `ena_run` ID (and possibly other
#     text strings).
#

library(dplyr) # arrange, select
library(rlang)  # Injection operator `!!`

# Automatically adapt to CSV or TSV format
read.xsv <- function(file, header = TRUE) {
  if (grepl(".csv$", file, ignore.case = TRUE)) {
    read.csv(file, header = header)
  } else if (grepl(".tsv$", file, ignore.case = TRUE)) {
    read.delim(file, header = header)
  } else {
    stop("Non-compliant file extension.")
  }
}

# Define a new binary operator (with 2 aliases) to concatenate strings in pipe
`%|+>%` <- `%+>%` <- `%+%` <- \(x,y)paste0(x,y)










new_bioModel <- function(target_dir) {

  # Get all files inside target directory
  file_list <- list.files(path = target_dir,
                          pattern = "\\.[ct]sv$",
                          full.names = FALSE)




  # Clean series IDs
  file_list |> sub("_.*$", "", x=_) |> sort() |> unique() -> series_IDs
  
  # Build up the bioModel object
  lapply(series_IDs, function(series_ID, files = file_list) {
    
    # Load data-metadata pairs (also sort metadata by `ena_run`)
    series_ID %+% "_CountMatrix.*" |> grep(files, ignore.case=T, value=T) |> read.xsv() -> counts_df
    series_ID %+% "_meta.*" |> grep(files, ignore.case=T, value=T) |> read.xsv() |> arrange(ena_run) -> meta_df
    
    # Convert rows to list
    meta_df |> split(seq(nrow(meta_df))) |> setNames(meta_df$ena_run) -> meta_list
    
    # Find ID column in `counts_df`
    "gene.*id|transcript.*id" |> grep(colnames(counts_df), ignore.case=T) -> ids_index
    
    # Build up a `series` object
    lapply(meta_list, function(run, coun_df = counts_df, id_ndx = ids_index) {
      # Look for run's count data...
      run$ena_run |> grep(colnames(coun_df)) -> coun_ndx
      # ...and add both counts (if present) and IDs to each run-list as a data frame
      coun_df |> select(IDs = !!id_ndx, counts = !!coun_ndx) |> list(genes=_) |> append(run, values=_)
    }) -> series
    
    # Add annotation to each series
    meta_df$ena_run |> paste(collapse = "|") |> grep(colnames(counts_df), invert=T) -> annot_index
    counts_df |> select(!!annot_index) |> list(annotation=_) |> append(series, values=_) -> series
    
    # Make the series list an S3 object and return it
    structure(series, class = "bioSeries")
    
  }) |> setNames(series_IDs) -> base_list
  
  # Make the base list an S3 object and return it
  structure(base_list, class = "bioModel")
}





countMatrix.bioSeries <- function(series, annot = FALSE) {
  # Find run elements
  grep("SRR", names(series), ignore.case=T) -> run_index
  
  # Extract counts, restore run ID names, then merge into one data frame
  series[run_index] |> lapply(\(run) {
    counts_df <- run$genes
    colnames(counts_df)[colnames(counts_df) == "counts"] <- run$ena_run
    counts_df
  }) |> Reduce(\(x, y) merge(x, y, by = "IDs", all = TRUE), x=_) -> count_matrix
  
  if (annot) {
    # Get annotation
    annot <- series$annotation
    "gene.*id|transcript.*id" |> grep(colnames(annot), ignore.case=T) -> ids_index
    count_matrix <- merge(annot, count_matrix,
                          by.x = ids_index, by.y = "IDs", all.y = TRUE)
  }
  return(count_matrix)
}





rowStats.bioSeries <- function(series, annot = FALSE) {
  
  # Get a log-transformed count matrix
  count_matrix <- countMatrix.bioSeries(series, annot = annot)
  count_index <- sapply(count_matrix, is.numeric)
  count_matrix[,count_index] <- log2(count_matrix[,count_index] + 1)
  
  # Compute descriptive statistics (mean, SD, SEM)
  average_counts <- rowMeans(count_matrix[,count_index], na.rm = TRUE)
  sd_counts <- apply(count_matrix[,count_index], 1, sd, na.rm = TRUE)
  sem_counts <- sd_counts/sqrt(sum(count_index))
  
  # Assemble results
  row_stats <- count_matrix[,!count_index]
  row_stats$Mean <- average_counts
  row_stats$Std_Dev <- sd_counts
  row_stats$SEM <- sem_counts
  
  return(row_stats)
}




a <- countMatrix.bioSeries(s$GSE195781, T)

b <- rowStats.bioSeries(s$GSE195781, T)




rowStats.bioModel <- function(model) {
  
}



