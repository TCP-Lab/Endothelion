#!/usr/bin/env -S Rscript --vanilla

# Endothelion Project
#
# Assumptions and format the of expected arguments.
#
# Count matrices are assumed to be the end product of the sequential application
# of many x.FASTQ modules. In particular, counts are supposed to be
#  - already normalized (ncounts), the normalization metric being the trailing
#    part of each sample name used for column heading;
#  - expressed in linear scale (not yet log-transformed);
#  - provided with ENSG IDs ('gene_id' column);
#  - provided with gene symbol annotation ('SYMBOL' column);
#
# In addition, the files containing the count matrices are supposed to be
#  - saved in TSV format;
#  - named according to the following pattern:
#           <SeriesID>_CountMatrix_<level>_<metric>.tsv 
#    E.g.,  GSE138309_CountMatrix_genes_TPM.tsv
#
# GOIs are assumed to be listed as plain gene symbols, arranged in a single
# column (CSV or TSV) with no header.
#
# When `threshold_adapt == "true"`, a GMM is used to find the expression
# threshold as the  decision boundary separating the two sub-populations of
# expressed and unexpressed genes. In this case, the integer entered as the
# `threshold_value` indicates the number of Gaussian components to be used in
# the mixture for the adaptive calculation of the expression threshold (thr). If
# the algorithm returns a thr value less than 1, thr is coerced to 1 by default.
# In contrast, when `threshold_adapt == "false"`, the `threshold_value` is the
# real value to be used as constant threshold in all the experiments.
#

# --- Packages -----------------------------------------------------------------

library(ggplot2)
library(r4tcpl)
# library(tidyr)
# library(httr)

# Ugly temporary solution while waiting to package the SeqLoader
dest_file <- "./src/seq_loader.R"
if (! file.exists(dest_file)) {
  "https://raw.githubusercontent.com/TCP-Lab/SeqLoader/main/seq_loader.R?token=GHSAT0AAAAAACMPIXKIYKZQIHPQUSZLEVAKZS7CMHA" |>
    httr::GET(httr::write_disk(dest_file, overwrite = TRUE)) -> download_status
}
source(dest_file)

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

# --- Gene Set -----------------------------------------------------------------

# Load the list of GOIs
gois_file |> read.delim(header = FALSE) -> gois

# NOTE
# Use 'r4tcpl::TGS' dataset to access the full transportome, or a subset of it,
# e.g.:
# gois <- r4tcpl::TGS$ICs

# --- Load xModel --------------------------------------------------------------

# Construct xModel from data
model <- new_xModel(target_dir)

# Select Runs and subset genes
model |> pruneRuns() |> keepRuns("extra == 1") |>
  subsetGenes("SYMBOL", gois) -> slim_model







# Threshold --------------------------------------------------------------------

# Adaptive expression threshold
if (threshold_adapt == "true") {
  # Subset the numeric columns and take their log2
  only_counts <- log2(ncounts[,-1] + 1)
  
  # Make box-plots of count distributions
  savePlots(
    \(){boxplot(only_counts)},
    figure_Name = paste0(GEO_id, "_boxplot"),
    figure_Folder = out_subdir,
    pdf_out = FALSE)
  
  # Find the expression threshold adaptively
  # Filter the dataset keeping only those genes that are detected in the majority
  # of the samples, compute their average expression, then use that distribution
  # of mean log-counts to fit the GMM.
  gmm <- GMM_divide(
    rowMeans(only_counts)[rowSums(only_counts > 0) > sample_size/2],
    G = threshold_value)
  
  # Set the new expression threshold as the right-most decision boundary
  thr <- gmm$boundary[threshold_value*(threshold_value-1)/2]
  
  # Make density plots with GMM overlaid
  savePlots(
    \(){
      # Density curves
      count_density(only_counts,
                    remove_zeros = TRUE,
                    xlim = c(-1,10),
                    titles = c(paste0("Kernel Density Plot\n", GEO_id), ""),
                    col = "gray20")
      # Plot the GMM
      for (i in 1:threshold_value) {
        lines(gmm$x, gmm$components[,i], col = "dodgerblue")
      }
      lines(gmm$x, rowSums(gmm$components), col = "firebrick2")
      # Plot the expression threshold
      y_lim <- par("yaxp")[2]
      lines(c(thr, thr), c(0, 1.5*y_lim), col = "darkslategray", lty = "dashed")
      original_adj <- par("adj") # Store the original value of 'adj'
      par(adj = 0) # Set text justification to left
      text(x = thr + 0.3, y = 0.8*y_lim,
           labels = paste("Decision Boundary =", round(thr, digits = 2)),
           cex = 1.1)
      par(adj = original_adj) # Restore the original 'adj' value
    },
    figure_Name = paste0(GEO_id, "_threshold"),
    figure_Folder = out_subdir)
  
  if (thr < 1) {
    cat("\nWARNING:\n Adaptive threshold from GMM returned",
        round(thr, digits = 2), "...been coerced to 1.")
    thr <- 1
  }
} else {
  # Fixed expression threshold (non-adaptive mode)
  thr <- threshold_value
}



# Intersection -----------------------------------------------------------------



ncounts <- geneStats(model$...)
gois_ncounts <- geneStats(slim_model$...)






# Statistics -------------------------------------------------------------------


gois_all_stats <- ncounts
  
# Saving as CSV
write.csv(gois_all_stats,
          file.path(out_subdir,
                    paste0(GEO_id, "_log2", count_type, "_profileReport.csv")),
          row.names = FALSE)

# Only expressed GOIs
gois_all_stats |> subset(Mean > thr) -> gois_high_stats

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

# Only expressed ICTs
gois_high_stats |> subset(Symbol %in% c(ICs, trans)) -> ICT_high_stats

# Only expressed ICs
gois_high_stats |> subset(Symbol %in% ICs) -> IC_high_stats

# Only expressed transporters
gois_high_stats |> subset(Symbol %in% trans) -> trans_high_stats

# Only expressed receptors (GPCRs and RTKs)
gois_high_stats |> subset(!Symbol %in% c(ICs, trans)) -> GPCRTK_high_stats

# Bar Chart --------------------------------------------------------------------

# Colors
line_color <- "gray17"
err_color <- "gray17"

# Limits
y_max <- max(gois_all_stats$Mean)
y_max_sd <- gois_all_stats$Std_Dev[gois_all_stats$Mean == y_max]
y_limit <- ceiling(y_max + y_max_sd)

gois_stats <- list(allGOIs = gois_all_stats,
                   highGOIs = gois_high_stats,
                   highICTs = ICT_high_stats,
                   highICs = IC_high_stats,
                   highTrans = trans_high_stats,
                   highGPCRTKs = GPCRTK_high_stats)

for (name in names(gois_stats)) {
  
  # Prepare the Frame
  gg_frame <-
    ggplot(data = gois_stats[[name]],
           aes(x = Symbol, y = Mean, fill = Symbol)) +
    theme_bw(base_size = 15, base_rect_size = 1.5) +
    theme(axis.text.x = element_text(size = 10, angle = 90,
                                     vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 14),
          legend.position = "none") +
    scale_y_continuous(expand = c(0.01, 0.02),
                       breaks = seq(0, y_limit, 0.5)) +
    xlab("Genes of Interest") +
    ylab(substitute(log[2]*(x+1), list(x = count_type))) +
    ggtitle(label = paste0(GEO_id, " (n = ", sample_size,")"))
  
  # Draw the Bars
  gg_bars <- gg_frame +
    geom_bar(stat = "identity", width = 0.75) +
    geom_errorbar(aes(ymin = Mean - Std_Dev,
                      ymax = Mean + Std_Dev),
                  linewidth = 1.0, width = 0.5, color = err_color)
  
  # # Alternative with borders
  # gg_bars <- gg_frame +
  #   geom_bar(stat = "identity", width = 0.7,
  #            color = line_color, linewidth = 0.1) +
  #   geom_errorbar(aes(ymin = Mean - Std_Dev,
  #                     ymax = Mean + Std_Dev),
  #                 linewidth = 1.0, width = 0.5, color = err_color)
  
  # Add the Expression Threshold
  gg_thr <- gg_bars +
    geom_hline(yintercept = thr,
               linetype = "dashed",
               color = line_color,
               linewidth = 1)
  
  # Save the Chart
  savePlots(
    \(){print(gg_thr)},
    width_px = 2000,
    figure_Name = paste0(GEO_id, "_log2", count_type, "_", name, "_chart"),
    figure_Folder = out_subdir)
}

# END --------------------------------------------------------------------------

cat("\n", GEO_id, " is done!\n", sep = "")







out_subdir <- file.path(out_dir, GEO_id)
# Set the output folder
if (! dir.exists(out_subdir)) {
  dir.create(out_subdir, recursive = TRUE)
}


# Load count matrix and check header
if (!("gene_id" %in% colnames(ncounts) && "SYMBOL" %in% colnames(ncounts))) {
  stop("\n Bad formatted count table... Stop executing.")
}


# Normalization check: sum(TPMs) == 10^6
if (any(abs(totalCounts(series) - 1e6) > 5)) {
  cat("\nWARNING:\n Bad TPM normalization... Check counts in matrix!\n")
}





if (dnues2(gois_ncounts$SYMBOL)[1] > 0) {
  stop("\n Cannot handle duplicated gene symbols...")
}

# Possibly collapse duplicated Gene Symbols (keeping the most expressed)
#DEGs <- DEGs[order(DEGs$adj_pval), ] ## from GOZER; to be adapted
#DEGs <- DEGs[!duplicated(DEGs$GeneSymbol), ]
