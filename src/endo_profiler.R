#!/usr/bin/env -S Rscript --vanilla

# Endothelion Project
#
# Assumptions and format the of expected arguments.
#
# This new version of 'endo_profiler' leverages S3 xSeries and xModel classes as
# provided by the SeqLoader package (https://github.com/TCP-Lab/SeqLoader). For
# this reason the formal requirements on data are the same of SeqLoader, the
# only additions being that counts are supposed to be
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
# library(tidyr)
# library(httr)
# library(AnnotationDbi)
# library(org.Hs.eg.db)

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
error_msg <- "\nERROR by endo_profiler_x.R\n"

# Check if the correct number of arguments is provided from command-line
if (length(commandArgs(trailingOnly = TRUE)) != 6) {
  cat(error_msg,
      "One or more arguments are missing. Usage:\n\n",
      "Rscript endo_profiler_x.R <in_path> <central_tendency> \\\n",
      "                        <threshold_adapt> <threshold_value> \\\n",
      "                        <GOIs> <out_dir>\n\n")
  quit(status = 1)
}

# Extract command-line arguments.
target_dir <- commandArgs(trailingOnly = TRUE)[1]
descriptive <- toupper(commandArgs(trailingOnly = TRUE)[2])
threshold_adapt <- commandArgs(trailingOnly = TRUE)[3]
threshold_value <- as.numeric(commandArgs(trailingOnly = TRUE)[4])
gois_file <- commandArgs(trailingOnly = TRUE)[5]
out_dir <- commandArgs(trailingOnly = TRUE)[6]

# Check if the target directory exists
if (! dir.exists(target_dir)) {
  cat(error_msg,
      " Directory \'", target_dir, "\' does not exist.\n", sep = "")
  quit(status = 2)
}

# Check adaptive threshold logical flag
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

# --- Per-Series Stats ---------------------------------------------------------

# Get series-specific stats for all GOIs and save as CSV
echo("\nSTEP 04 :: Per-Series Stats", "green")
slim_model |> lapply(\(series) {
  series_ID <- attr(series, "own_name")
  series |> geneStats(annot = TRUE) -> all_gois_stats
  report_name <- paste0(series_ID, "_profileReport.csv")
  cat("\nSaving:", report_name)
  write.csv(all_gois_stats,
            file.path(out_dir, series_ID, report_name),
            row.names = FALSE)
  return(all_gois_stats)
}) |> set_own_names() -> all_gois_stats
cat("\n")

# --- Bar Charting -------------------------------------------------------------

# Make subgroups
# NOTE: use dplyr::filter instead of base::subset to preserve attributes!

# Only expressed (aka "high") GOIs
all_gois_stats |> names() |> sapply(\(series_name) {
    all_gois_stats[[series_name]] |> filter(Mean > thr[[series_name]])
  }, simplify = FALSE, USE.NAMES = TRUE) -> high_gois_stats

# A (temporary) patch for the TGS dataset from r4tcpl v.1.5.1
patch4TGS <- c("MCUB",
               "MCUR1",
               "KCNRG",
               "KCNIP1",
               "KCNIP2",
               "KCNIP3",
               "KCNIP4",
               "CATSPERE",
               "CATSPERZ")
ICs <- c(TGS$ICs, patch4TGS)
trans <- TGS$trans

# Only expressed ICTs, ICs, transporters, and receptors (GPCRs and RTKs)
high_gois_stats |> lapply(filter, SYMBOL %in% c(ICs, trans)) -> high_ICT_stats
high_gois_stats |> lapply(filter, SYMBOL %in% ICs) -> high_IC_stats
high_gois_stats |> lapply(filter, SYMBOL %in% trans) -> high_trans_stats
high_gois_stats |> lapply(filter, !SYMBOL %in% c(ICs, trans)) -> high_rex_stats

# Make a comprehensive named structure (list of lists)
gois_stats <- list(all_GOIs = all_gois_stats,
                   high_GOIs = high_gois_stats,
                   high_ICTs = high_ICT_stats,
                   high_ICs = high_IC_stats,
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

echo("\nSTEP 06 :: Model Synthesis", "green")

slim_model |> geneStats(descriptive = eval(parse(text = descriptive)),
                        maic = "inclusive",
                        annot = FALSE) -> average_expression

cat("\nRegenerating gene annotation from scratch...")
# Change to `OrgDb_key="ENSEMBLTRANS"` when working at isoform level
average_expression |> add_annotation(OrgDb_key="ENSEMBL") -> average_expression

# Save full GOI set as CSV
write.csv(average_expression,
          file.path(out_dir,
                    paste0(attr(model,"own_name"),"_AverExpress_ALL_GOIs.csv")),
          row.names = FALSE)





# --- correlation ----

# Correlation ScatterPlot Matrix (with fixed text size)
custom_pairs <- function(data_set, color = "gray15") {
  
  # Customize lower panel (correlation values)
  panel_cor <- function(x, y) {
    default_usr <- par("usr")
    on.exit(par(usr = default_usr))
    par(usr = c(0, 1, 0, 1))
    r <- round(cor(x, y), digits = 3)
    text(0.5, 0.5, r, cex = 2)
  }
  
  # Customize upper panel (scatter plots)
  panel_points <- function(x, y) {
    points(x, y, pch = 19, cex = 1, col = color)
  }
  
  # Create the plots
  pairs(data_set,
        cex.labels = 3,
        font.labels = 4,
        lower.panel = panel_cor,
        upper.panel = panel_points)
}

# Correlation ScatterPlot Matrix (with fixed text size)
custom_pairs(data_set = data.frame(GSE138309 = gois_stats$all_GOIs$GSE138309$Mean,
                                   GSE139133 = gois_stats$all_GOIs$GSE139133$Mean,
                                   GSE195781 = gois_stats$all_GOIs$GSE195781$Mean,
                                   GSE205739 = gois_stats$all_GOIs$GSE205739$Mean,
                                   GSE76528 = gois_stats$all_GOIs$GSE76528$Mean),
             color = "steelblue4")






# Sort the average reference by decreasing ICT expression
average_expression |>
  select(c("SYMBOL", "Mean")) |> arrange(desc(Mean)) -> reference

# Convert SYMBOL to a factor to ensure that the x-axis values follow the order
# of the sorted dataframe instead of being alphabetically ordered
reference$SYMBOL <- factor(reference$SYMBOL, levels = reference$SYMBOL)


gois_stats$all_GOIs$GSE138309 |> select(c("SYMBOL", "Mean")) -> df1
gois_stats$all_GOIs$GSE139133 |> select(c("SYMBOL", "Mean")) -> df2
gois_stats$all_GOIs$GSE195781 |> select(c("SYMBOL", "Mean")) -> df3

# By using `levels(reference$SYMBOL)` as levels for all other factors ensures
# that all of them will be plotted following the gene ranking of the reference
# (i.e., the average) data frame.
df1$SYMBOL <- factor(df1$SYMBOL, levels = levels(reference$SYMBOL))
df2$SYMBOL <- factor(df2$SYMBOL, levels = levels(reference$SYMBOL))
df3$SYMBOL <- factor(df3$SYMBOL, levels = levels(reference$SYMBOL))

# Combine all data into one dataframe for easier plotting
combined <- bind_rows(
  mutate(reference, Source = "Reference"),
  mutate(df1, Source = "DF1"),
  mutate(df2, Source = "DF2"),
  mutate(df3, Source = "DF3")
)


# Plot the reference dataframe
# Prepare the Frame
gg_frame <-
  ggplot(combined, aes(x = SYMBOL, y = Mean, group = Source, color = Source)) + 
  theme_bw(base_size = 15, base_rect_size = 1.5) +
  theme(axis.text.x = element_text(size = 10, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "none") +
  xlab("Genes of Interest") +
  ylab(substitute(log[2]*(x+1), list(x = "TPM"))) +
  ggtitle(label = "Fading Trend")

gg_points <- gg_frame + geom_point() # + geom_line()



print(gg_points)







# Functional Subsetting --------------------------------------------------------

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
  # Check for completeness
  if (length(setdiff(subGOIs_list[[file_name]],average_expression$SYMBOL)) > 0) {
    cat("\nWARNING by ", file_label, ":",
        "\n Missing these gene symbols in GOI main list:\n  ", sep = "")
    cat(setdiff(subGOIs_list[[file_name]],average_expression$SYMBOL), sep="\n  ")
  }
}



# --- END ----------------------------------------------------------------------
echo(paste0("\n", attr(model, "own_name"), " is done!\n"), "green")
















