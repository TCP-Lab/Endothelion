
# Assumptions:
#
# 1 - Both data and metadata have file names starting with the same series
#     (experiment, project) ID (both GEO and ENA are fine), followed by an
#     underscore ("_") and either the distinctive string "CountMatrix" or "meta"
#     if they are counts or metadata, respectively. Additional characters are
#     allowed after this pattern, including either CSV or TSV file extensions.
#     Examples of good file name patterns are:
#
#         GSE205739_CountMatrix_genes_TPM.tsv
#         GSE76528_meta.csv
#
# 2 - Metadata files have proper column names as table header, among which
#     `ena_run` is mandatory, as the only ID that is certainly unique to each
#     file.
#
# 3 - Count tables have proper column names as header, among which a mandatory
#     ID column mtching the regex "gene.*id" or "transcript.*id" (usually being
#     "gene_id" or "transcript_id", respectively).
#
# 4 - Columns in cout table are named according to the `ena_run` ID




library(dplyr)


setwd("../data/in/Lines/hCMEC_D3/")


target_dir <- "."

file_list <- list.files(path = target_dir,
                        pattern = "\\.[ct]sv$",
                        full.names = FALSE)


file_list |> sub("_.*$", "", x=_) |> unique() -> series_IDs


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

# Define a new binary operator (with an alias) to concatenate strings in pipe
`%|+>%` <- `%+>%` <- \(x,y)paste0(x,y)





lapply(series_IDs, function(series_ID, files = file_list) {
  
  # Load count data and associated metadata
  series_ID %|+>% "_CountMatrix.*" |> grep(files, value=T) |> read.xsv() -> counts_df
  series_ID %|+>% "_meta.*" |> grep(files, value=T) |> read.xsv() -> meta_df
  
  # Sort by `ena_run` and convert rows to list
  meta_df |> arrange(ena_run) -> meta_df
  meta_df |> split(seq(nrow(meta_df))) |> setNames(meta_df$ena_run) -> meta_list
  
  # Find the ID column in count table
  grep("gene.*id|transcript.*id", colnames(counts_df)) -> ids_index
  
  # Add read counts and gene IDs to each run in `meta_list` to build up the
  # `series` object
  lapply(meta_list, function(run, counts_df = counts_df, ids_index = ids_index){
    # Find the correct SRR run column in count table
    run$ena_run |> grep(colnames(counts_df)) -> counts_index
    # Add gene IDs and read count data (if present) as a data frame to each run 
    run |> append(list(gene = data.frame(IDs = counts_df[ids_index],
                                         counts = counts_df[counts_index])))
    }) -> series
  
  
  meta_df$ena_run |> paste(collapse = "|") |>
    grep(colnames(counts_df), invert = TRUE) -> annot_index
  
  series |> append(list(annotation = counts_df[,annot_index]))

}) |> setNames(series_IDs) -> endoModel



