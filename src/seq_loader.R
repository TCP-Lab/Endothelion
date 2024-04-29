setwd("../data/in/Lines/hCMEC_D3/")


target_dir <- "."

file_list <- list.files(path = target_dir,
                        pattern = "\\.[ct]sv$",
                        full.names = FALSE)

# Both data and metadata file names have to start with the same series
# (experiment) ID, followed by an underscore and either "CountMatrix" or "meta"
# if counts or metadata, respectively. More characters are allowed after them,
# including ANY file extension.
# meta need to have proper column names as heading, with ena_run as mandatory
#
file_list |> sub("_.*$", "", x=_) |> unique() -> series_IDs



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
  
  files |> grep(paste0(series_ID, "_CountMatrix.*"), x=_, value = TRUE) |>
    read.xsv() -> counts_df
  
  files |> grep(paste0(series_ID, "_meta.*"), x=_, value = TRUE) |>
    read.xsv() -> meta_df
  
  # Convert rows to list
  meta_df |> split(seq(nrow(meta_df))) |> setNames(meta_df$ena_run) -> meta_list
  
  grep("gene.*id|transcript.*id", colnames(counts_df)) -> ids_index
  
  lapply(meta_list, function(run){
    
    run$ena_run |> grep(colnames(counts_df)) -> counts_index
    
    run |> append(list(
      gene = data.frame(IDs = counts_df[ids_index],
                        counts = counts_df[counts_index])))
    
    }) -> series
  
  meta_df$ena_run |> paste(collapse = "|") |>
    grep(colnames(counts_df), invert = TRUE) -> annot_index
  
  series |> append(list(annotation = counts_df[,annot_index]))

}) |> setNames(series_IDs) -> endoModel



