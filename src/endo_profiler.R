#!/usr/bin/env Rscript

# Endothelion Project
#
# Assumptions and formats of expected arguments.
#
# This second version of 'endo_profiler' leverages xSeries and xModel S3 classes
# as provided by the SeqLoader package (https://github.com/TCP-Lab/SeqLoader).
# For this reason, the formal requirements on data are the same of SeqLoader,
# the only additions for Endothelion being that counts are supposed to be
#  - already normalized in TPM;
#  - provided with gene symbol annotation ('SYMBOL' column);
#
# GOIs are assumed to be listed as plain gene symbols, arranged in a single
# column (CSV or TSV) with no header.
#
# When `threshold_adapt == "true"`, a GMM is used to find the expression
# threshold as the decision boundary separating the two sub-populations of
# expressed and unexpressed genes. In this case, the integer entered as the
# `threshold_value` indicates the number of Gaussian components to be used in
# the mixture for the adaptive calculation of the expression threshold (thr). If
# the algorithm returns a thr value less than 1, thr is coerced to 1 by default.
# In contrast, when `threshold_adapt == "false"`, the `threshold_value` is the
# real value to be used as constant threshold in all the experiments.

# --- Packages -----------------------------------------------------------------

library(ggplot2)
library(r4tcpl)
library(dplyr, warn.conflicts = FALSE)
library(tidyr)
library(DBI)
library(RSQLite)
# library(httr)

# Ugly temporary solution while waiting to package the SeqLoader
dest_file <- "./src/seq_loader.R"
if (! file.exists(dest_file)) {
  "https://raw.githubusercontent.com/TCP-Lab/SeqLoader/main/seq_loader.R" |>
    httr::GET(httr::write_disk(dest_file, overwrite = TRUE)) -> download_status
}
source(dest_file) |> suppressMessages()

# Function loading
source("./src/endo_functions.R")

# --- Input Parsing ------------------------------------------------------------

# General error message
error_msg <- "\nERROR by endo_profiler.R\n"

# Check if the correct number of arguments is provided from command-line
if (length(commandArgs(trailingOnly = TRUE)) != 7) {
  cat(error_msg,
      "One or more arguments are missing. Usage:\n\n",
      "Rscript endo_profiler.R <in_path> <central_tendency> \\\n",
      "                        <threshold_adapt> <threshold_value> \\\n",
      "                        <GOIs> <db_path> <out_dir>\n\n")
  quit(status = 1)
}

# Extract command-line arguments.
target_dir <- commandArgs(trailingOnly = TRUE)[1]
descriptive <- toupper(commandArgs(trailingOnly = TRUE)[2])
threshold_adapt <- commandArgs(trailingOnly = TRUE)[3]
threshold_value <- as.numeric(commandArgs(trailingOnly = TRUE)[4])
gois_file <- commandArgs(trailingOnly = TRUE)[5]
db_path <- commandArgs(trailingOnly = TRUE)[6]
out_dir <- commandArgs(trailingOnly = TRUE)[7]

# # Interactive debug (from the project root directory)
# target_dir <- "./data/in/Lines/test_hCMEC_D3/"
# target_dir <- "./data/in/Lines/hCMEC_D3/"
# descriptive <- "MEAN"
# threshold_adapt <- "false"
# threshold_value <- 1
# gois_file <- "./data/in/ICT_set_v2.csv"
# db_path <- "./data/MTPDB.sqlite"
# out_dir <- "./data/out/Lines/test_hCMEC_D3/"

# Check if the target directory exists
if (! dir.exists(target_dir)) {
  cat(error_msg,
      " Directory \'", target_dir, "\' does not exist.\n", sep = "")
  quit(status = 2)
}

# Check central tendency descriptor
if (! descriptive %in% c("MEAN", "MEDIAN", "NWMEAN")) {
  cat(error_msg,
      " Invalid \'descriptive\' parameter \'", descriptive, "\'.\n",
      " It must be one of the central tendency metrics currently implemented",
      " by SeqLoader.\n", sep = "")
  quit(status = 3)
}

# Check adaptive threshold logical flag
if (! threshold_adapt %in% c("true", "false")) {
  cat(error_msg,
      " Invalid \'threshold_adapt\' parameter \'", threshold_adapt, "\'.\n",
      " It must be one of the two Bash logical values true or false.\n", sep="")
  quit(status = 4)
}

# Check threshold_value input type
if (is.na(threshold_value)) {
  cat(error_msg,
      " Invalid data type for \'threshold_value\' parameter.\n",
      " A single numeric value is expected here.\n", sep = "")
  quit(status = 5)
}

# Check if the list of the Genes of Interest (GOIs) exists.
if (! file.exists(gois_file)) {
  cat(error_msg,
      " File \'", gois_file, "\' does not exist.\n", sep = "")
  quit(status = 6)
}

# Check if the list of the Genes of Interest (GOIs) exists.
if (! file.exists(db_path)) {
    cat(error_msg,
        " File \'", db_path, "\' does not exist.\n", sep = "")
    quit(status = 7)
}

# --- xModel Loading -----------------------------------------------------------

# Construct xModel from data
echo("\nSTEP 01 :: xModel Loading", "green")
model <- new_xModel(target_dir)

echo("\nFull-Model facts", "yellow")
cat("Source   :", target_dir, "\n")
factTable(model)

echo("\nReduced-Model facts", "yellow")
model |> pruneRuns() |> keepRuns("extra == 1") -> model
factTable(model)

# Normalization check: sum(TPMs) == 10^6
model |> sapply(\(series) {
  if (!all((abs(totalCounts(series) - 1e6) < 5) |
           (totalCounts(series) == 0))) {
    cat("\nWARNING:\n Bad TPM normalization in series ",
        attr(series, "own_name"), "...\n Check counts in the matrix!\n\n",
        sep = "")
    totalCounts(series) |> print()
  }
}) |> invisible()

# Annotation check: Endothelion needs SYMBOLs
model |> sapply(\(series) {
  if(!hasName(series$annotation, "SYMBOL")) {
    stop("\nMissing gene symbol annotation... Stop execution.")
  }
}) |> invisible()

# Per-Series detailed report
echo("\nDetailed report", "yellow")
model |> sapply(factTable) |> invisible()

# --- Threshold Computation ----------------------------------------------------

# Compute threshold
echo("\nSTEP 02 :: Threshold Computation", "green")
model |> lapply(threshold, names(model)) -> thr

# --- GOI Extraction -----------------------------------------------------------

# Intersect with GOIs
echo("\nSTEP 03 :: GOI Extraction", "green")

# Load the list of GOIs
# NOTE: use 'r4tcpl::TGS' to access the full transportome, or a subset of it!
gois_file |> read.delim(header = FALSE) -> gois

echo("\nGenes of Interest", "yellow")
cat("Source:", gois_file, "\nLoaded a list of", nrow(gois), "GOIs\n")

# Select Runs and subset genes
model |> subsetGenes("SYMBOL", gois) -> slim_model

echo("\nSlim model facts", "yellow")
factTable(slim_model)

# Endothelion works at SYMBOL level. Remove possible duplicated gene symbols
# (i.e., multiple IDs mapping to the same SYMBOL) by keeping the most expressed
# one. To do this, sort duplicated SYMBOLS by descending (grand) 'Mean' and use
# 'distinct' to keep the first entry only.
#
# NOTE: here slim_model[[1]]$annotation is used for annotation, while waiting
#       for 'Model-level annotation synthesis'.
dnues(slim_model[[1]]$annotation$SYMBOL)[1] -> duplicated_SYMBOLS
if (duplicated_SYMBOLS > 0) {
    cat("\n", duplicated_SYMBOLS, "duplicated gene symbol(s) found! Collapsing...\n")
    
    merge(slim_model[[1]]$annotation, geneStats(slim_model),
          by = "IDs", all.y = TRUE) |>
        arrange(SYMBOL, desc(Mean)) |> distinct(SYMBOL, .keep_all = TRUE) |>
        select(IDs) -> IDs_of_unique_SYMBOLs
    
    subsetGenes(slim_model, "IDs", IDs_of_unique_SYMBOLs) -> slim_model
    echo("\nSlim model facts", "yellow")
    factTable(slim_model)
}

# --- Per-Series Stats ---------------------------------------------------------

# Get series-specific stats for all GOIs and save as CSV
echo("\nSTEP 04 :: Per-Series Stats", "green")
slim_model |> lapply(\(series) {
  series |> geneStats(annot = TRUE) -> all_gois_stats
  series_ID <- attr(series, "own_name")
  report_name <- paste0(series_ID, "_profileReport.csv")
  # Series sub-dirs should be already there from STEP 00, but just in case...
  dir.create(file.path(out_dir, series_ID),
             recursive = TRUE, showWarnings = FALSE)
  cat("\nSaving:", report_name)
  write.csv(all_gois_stats,
            file.path(out_dir, series_ID, report_name),
            row.names = FALSE)
  return(all_gois_stats)
}) |> set_own_names() -> all_gois_stats
cat("\n")

# --- Bar Charting -------------------------------------------------------------

# Make subgroups by retrieving data from the MTP-DB

if (Sys.info()["sysname"] == "Windows") {
    # Can't access the DB directly on WSL when running from Windows...
    db_path <- file.path(Sys.getenv("USERPROFILE"), "Desktop", "MTPDB.sqlite")
    file.copy(from = "./data/MTPDB.sqlite", to = db_path)
}

# Connect to the MTP-DB
connection <- dbConnect(SQLite(), dbname = db_path)

# Pores (Ion Channels + AQPs)
query_pores <-
    "SELECT DISTINCT
    	gene_names.hugo_gene_symbol
    FROM
    	channels JOIN gene_names ON channels.ensg = gene_names.ensg
    UNION SELECT DISTINCT
    	gene_names.hugo_gene_symbol
    FROM
    	aquaporins JOIN gene_names ON aquaporins.ensg = gene_names.ensg"

# Transporters (SLCs + pumps)
query_trans <-
    "SELECT DISTINCT
    	gene_names.hugo_gene_symbol
    FROM
    	solute_carriers JOIN gene_names ON solute_carriers.ensg = gene_names.ensg
    UNION SELECT DISTINCT
    	gene_names.hugo_gene_symbol
    FROM
    	pumps JOIN gene_names ON pumps.ensg = gene_names.ensg"

# Make the calls
pores <- dbGetQuery(connection, query_pores)
trans <- dbGetQuery(connection, query_trans)

# Check sets here!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Disconnect from the MTP-DB
dbDisconnect(connection)

# Apply patches (for MTP-DB ver. 1.25.24)
pores |> unlist() |> c("CATSPERG",
                       "CATSPERB",
                       "CATSPERD",
                       "CATSPERE",
                       "CATSPERZ") |> unname() -> pores
trans |> unlist() |> c("SLC4A6") |> unname() -> trans

# Only expressed (aka "high") GOIs
all_gois_stats |> names() |> sapply(\(series_name) {
    all_gois_stats[[series_name]] |> filter(Mean > thr[[series_name]])
}, simplify = FALSE, USE.NAMES = TRUE) -> high_gois_stats

# Only expressed ICTs, pores, transporters, and receptors (GPCRs and RTKs)
# NOTE: use dplyr::filter instead of base::subset to preserve attributes!
high_gois_stats |> lapply(filter, SYMBOL %in% c(pores, trans)) -> high_ICT_stats
high_gois_stats |> lapply(filter, SYMBOL %in% pores) -> high_pores_stats
high_gois_stats |> lapply(filter, SYMBOL %in% trans) -> high_trans_stats
high_gois_stats |> lapply(filter, !SYMBOL %in% c(pores, trans)) -> high_rex_stats

# Make a comprehensive named structure (list of lists)
gois_stats <- list(all_GOIs = all_gois_stats,
                   high_GOIs = high_gois_stats,
                   high_ICTs = high_ICT_stats,
                   high_pores = high_pores_stats,
                   high_Transporters = high_trans_stats,
                   high_Receptors = high_rex_stats)

# Draw a bar chart for each Series
echo("\nSTEP 05 :: Series Bar Charting", "green")

# Set Limits
all_gois_stats |> sapply(\(series)series$Mean |> max()) |> which.max() -> s_indx
all_gois_stats[[s_indx]]$Mean |> which.max() -> g_indx
all_gois_stats[[s_indx]]$Mean[g_indx] -> y_max
all_gois_stats[[s_indx]]$Std_Dev[g_indx] -> y_max_sd
y_limit <- y_max + y_max_sd

# Draw a bar chart for each Series
gois_stats |> names() |> lapply(\(family_name) {
  echo(paste("\nCharting", family_name), "yellow")
  gois_stats[[family_name]] |>
    lapply(plot_barChart, family_name, y_limit, border = FALSE, thr)
}) |> invisible()

# --- Meta-analysis ------------------------------------------------------------

echo("\nSTEP 06 :: Model Synthesis (Meta-analysis)", "green")

slim_model |> geneStats(descriptive = eval(parse(text = descriptive)),
                        maic = "inclusive",
                        annot = FALSE) -> average_expression

cat("\nAdding gene annotation...\n")
# NOTE: again slim_model[[1]]$annotation is used for annotation, while waiting
#       for 'Model-level annotation synthesis'...
merge(slim_model[[1]]$annotation, average_expression,
      by = "IDs", all.y = TRUE) -> average_expression
    
# # Regenerating annotation from scratch may introduce some duplication issue...
# cat("\nRegenerating gene annotation from scratch...")
# # Change to `OrgDb_key="ENSEMBLTRANS"` when working at isoform level
# # NOTE: Possible 1:many mapping in IDs:SYMBOL will results in duplicated IDs and
# #       (likely) non-GOI genes. Thus, `filter` is used to remove genes not in
# #       set of GOIs (and non-unique IDs).
# average_expression |> add_annotation(OrgDb_key="ENSEMBL") |>
#     filter(SYMBOL %in% unlist(gois)) -> average_expression

# Save full GOI set as CSV
write.csv(average_expression,
          file.path(out_dir,
                    paste0(attr(model,"own_name"),"_AverExpress_ALL_GOIs.csv")),
          row.names = FALSE)

# --- Between Series Correlation -----------------------------------------------

echo("\nSTEP 07 :: Between Series Correlation", "green")

# Get the mean expression values for all the series
gois_stats$all_GOIs |> lapply(select, c("SYMBOL", "Mean")) |>
  lapply(\(series) {
    colnames(series)[colnames(series) == "Mean"] <- attr(series, "own_name")
    return(series)
  }) |> Reduce(\(x,y) merge(x,y, by="SYMBOL", all=TRUE), x=_) -> matrix_of_means

# Plot and save the pairwise correlation Chart
plot_label <- "Correlation_chart"
r4tcpl::savePlots(
  \(){custom_pairs(matrix_of_means[-1], color = "steelblue4")},
  width_px = 2000,
  figure_Name = plot_label,
  figure_Folder = out_dir)

n <- dim(matrix_of_means)[2] - 1
cat(paste0("\nSaving: ", plot_label, " (", n, "-by-", n, ")"))

# Set the expression cutoff (set filter to 0 to plot all the GOIs), merge the
# `average_expression` with the `matrix_of_means`, and sort by decreasing
# average ICT expression.
average_expression |>
    select(c("SYMBOL", "Mean")) |>
    filter(Mean >= 1) |>
    merge(matrix_of_means, by = "SYMBOL", all.y = FALSE) |>
    arrange(desc(Mean)) -> matrix_of_means

# Convert SYMBOL to a factor to ensure that the x-axis values follow the order
# of the sorted dataframe instead of being alphabetically ordered
matrix_of_means$SYMBOL <- factor(matrix_of_means$SYMBOL,
                                 levels = matrix_of_means$SYMBOL)

# Melt the data frame (convert to long format, for ggplot)
matrix_of_means |> pivot_longer(cols = !matches("SYMBOL"),
                                names_to = "Source",
                                values_to = "Mean") -> matrix_of_means

# Make the 'expression fading plots'...
ct <- "95CI" # chart_type argument
profile_plots <- list(
  Profile_global = profilePlot(matrix_of_means, ct),
  Profile_pores  = profilePlot(matrix_of_means |> filter(SYMBOL %in% pores), ct),
  Profile_trans  = profilePlot(matrix_of_means |> filter(SYMBOL %in% trans), ct),
  Profile_rex    = profilePlot(matrix_of_means |> filter(!SYMBOL %in% c(pores, trans)), ct))
# ...and save them
for (name in names(profile_plots)) {
  r4tcpl::savePlots(
    \(){print(profile_plots[[name]])},
    width_px = 1024,
    ratio = 1/sqrt(2), # A4 ratio
    #ratio = 1/(nrow(profile_plots[[name]]$data)/500 + 1), # empirical relation
    figure_Name = paste0(name, "_", ct),
    figure_Folder = out_dir)
  cat(paste("\nSaving:", name, "-", ct))
}
cat("\n")

# Functional Subsetting --------------------------------------------------------

echo("\nSTEP 08 :: Functional Subsetting", "green")

# Get the list of functional subsets to extract
subGOIs_prefix <- "ICT_subset_"
subGOIs_files <- list.files(path = dirname(gois_file),
                            pattern = paste0("^", subGOIs_prefix, ".*\\.csv$"),
                            full.names = TRUE, recursive = TRUE)

# Read the files into separate data frames and store them in a named list
subGOIs_files |> sapply(read.csv, header = FALSE) -> subGOIs_list

# Extract GOI subsets and save them with distinctive file names
for (file_name in names(subGOIs_list)) {
  # Modify input file names to get output names
  file_name |> basename() |>
    sub(subGOIs_prefix, "AverExpress_", x=_) |>
    sub("\\.csv\\..*$", ".csv", x=_) -> file_label
  # Subset and save as CSV
  average_expression |> filter(SYMBOL %in% subGOIs_list[[file_name]]) |>
    write.csv(file = file.path(out_dir, file_label), row.names = FALSE)
  cat("\nSaving:", file_label)
  # Check for completeness
  if (length(setdiff(subGOIs_list[[file_name]],average_expression$SYMBOL)) > 0) {
    cat("\nWARNING by ", file_label, ":",
        "\n Missing gene symbols in GOI main list:\n  ", sep = "")
    cat(setdiff(subGOIs_list[[file_name]],average_expression$SYMBOL), sep="\n  ")
  }
}
cat("\n")

# --- END ----------------------------------------------------------------------
echo(paste0("\n", attr(model, "own_name"), " is done!\n"), "green")







