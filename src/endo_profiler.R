#
# Endothelion Project
#
# Assumptions and format of expected arguments.
#
# Count matrices are assumed to be the end product of the sequential application
# of the different modules from the x.FASTQ framework. In particular, counts are
# supposed to be
#  - already normalized (ncounts), the normalization metric being the trailing
#    part of each sample name used for column heading;
#  - expressed in linear scale (not yet log-transformed);
#  - provided with ENSG IDs ('gene_id' column);
#  - provided with gene symbol annotation ('SYMBOL' column);
#
# In addition, the files containing the count matrices are supposed to be
#  - saved in TSV format;
#  - named according to the following pattern:
#
#       <GEOid>_<metric>_<depth>_<PE/SE>.tsv
#
#    E.g., a study with an average depth of 40 megareads per sample, paired-end,
#
#       GSE76528_TPM_40_PE.tsv
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

# Packages ---------------------------------------------------------------------

library(ggplot2)
library(r4tcpl)
# library(tidyr)

# Input Loading-----------------------------------------------------------------

# General error message
error_msg <- "\nERROR by endo_profiler.R\n"

# Check if the correct number of arguments is provided from command-line
if (length(commandArgs(trailingOnly = TRUE)) != 6) {
  cat(error_msg,
      "One or more arguments are missing. Usage:\n\n",
      "Rscript endo_profiler.R <count_matrix> <count_type> \\\n",
      "                        <threshold_adapt> <threshold_value> \\\n",
      "                        <GOIs> <out_dir>\n\n")
  quit(status = 1)
}

# Extract command-line arguments.
count_file <- commandArgs(trailingOnly = TRUE)[1]
count_type <- commandArgs(trailingOnly = TRUE)[2]
threshold_adapt <- commandArgs(trailingOnly = TRUE)[3]
threshold_value <- as.numeric(commandArgs(trailingOnly = TRUE)[4])
gois_file <- commandArgs(trailingOnly = TRUE)[5]
out_dir <- commandArgs(trailingOnly = TRUE)[6]

# Check if the target file exists
if (! file.exists(count_file)) {
 cat(error_msg,
     " File \'", count_file, "\' does not exist.\n", sep = "")
 quit(status = 2)
}

# Check the metric
if (! count_type %in% c("expected_count", "TPM", "FPKM")) {
  cat(error_msg,
      " Unknown metric \'", count_type, "\'\n", sep = "")
  quit(status = 3)
}

# Check adaptive threshold logical flag
if (! threshold_adapt %in% c("true", "false")) {
  cat(error_msg,
      " Invalid \'threshold_adapt\' parameter \'", threshold_adapt, "\'.\n",
      " It must be one of the two Bash logical values true or false.\n", sep = "")
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

# Dirs & Bases -----------------------------------------------------------------

# Set the output folder
if (! dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Retrieve the GEO ID from input file name
count_file |> basename2() |> {\(x)strsplit(x,"_")[[1]][1]}() -> GEO_id

# Count Data -------------------------------------------------------------------

# Load and heading check
count_file |> read.delim() -> ncounts
if (!("gene_id" %in% colnames(ncounts) && "SYMBOL" %in% colnames(ncounts))) {
  stop("ERROR: Bad formatted count table... Stop executing.")
}

# Subset
regex <- paste0("^SYMBOL$|_", count_type, "$")
ncounts <- ncounts[,grep(regex, colnames(ncounts))]
sample_size <- dim(ncounts)[2] - 1

# Normalization check: sum(TPMs) == 10^6
if (count_type == "TPM" && any(abs(colSums(ncounts[,-1]) - 1e6) > 5)) {
  stop("ERROR: Bad TPM normalization...Stop executing.")
}

# Explode the dataframe (collapsed by cc_assembler.R from x.FASTQ)
# tidyr::separate_rows(ncounts, SYMBOL, sep = ",") |> as.data.frame() -> ncounts ###########################

# Threshold --------------------------------------------------------------------

# Adaptive expression threshold
if (threshold_adapt == "true") {
  # Subset the numeric columns and take their log2
  only_counts <- log2(ncounts[,-1] + 1)
  
  # Make box-plots of count distributions
  savePlots(
    \(){boxplot(only_counts)},
    figure_Name = paste0(GEO_id, "_boxplot"),
    figure_Folder = out_dir,
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
    figure_Folder = out_dir)
  
  if (thr < 1) {
    cat("\nWARNING:\n Adaptive threshold from GMM returned",
        round(thr, digits = 2), "...been coerced to 1.")
    thr <- 1
  }
} else {
  # Fixed expression threshold (non-adaptive mode)
  thr <- threshold_value
}

# Gene Set ---------------------------------------------------------------------

# Load
gois_file |> read.delim(header = FALSE) |> unlist() -> gois

# NOTE
# Use 'r4tcpl::TGS' dataset to access the full transportome, or a subset of it
# E.g., gois <- r4tcpl::TGS$ICs

# Intersection -----------------------------------------------------------------

# Row subsetting
gois_ncounts <- subset(ncounts, SYMBOL %in% gois)

if (setdiff(gois, ncounts$SYMBOL) |> length() > 0) {
  cat("\nWARNING:\n Can't find these Genes of Interest in the Count Matrix:",
      setdiff(gois, ncounts$SYMBOL), sep = "\n  ")
}

if (dnues2(gois_ncounts$SYMBOL)[1] > 0) {
  stop("ERROR: Cannot handle duplicated gene symbols...")
}

# Possibly collapse duplicated Gene Symbols (keeping the most expressed)
#DEGs <- DEGs[order(DEGs$adj_pval), ] ## from GOZER; to be adapted
#DEGs <- DEGs[!duplicated(DEGs$GeneSymbol), ]

# Statistics -------------------------------------------------------------------

# Log-transform
gois_ncounts[,-1] <- log2(gois_ncounts[,-1] + 1)

# Compute descriptive statistics (mean, SD, SEM)
average_ncounts <- rowMeans(gois_ncounts[,-1], na.rm = TRUE)
sd_ncounts <- apply(gois_ncounts[,-1], 1, sd, na.rm = TRUE)
sem_ncounts <- sd_ncounts/sqrt(sample_size)

# Final expression matrix
gois_expression <- data.frame(Symbol = gois_ncounts$SYMBOL,
                              Mean = average_ncounts,
                              Std_Dev = sd_ncounts,
                              SEM = sem_ncounts)

# Saving as CSV
write.csv(gois_expression,
          file.path(out_dir,
                    paste0(GEO_id, "_log2", count_type, "_profileReport.csv")))

# Bar Chart --------------------------------------------------------------------

# Colors
line_color <- "gray17"
err_color <- "gray17"

# Limits
y_max <- max(gois_expression$Mean)
y_max_sd <- gois_expression$Std_Dev[gois_expression$Mean == y_max]
y_limit <- ceiling(y_max + y_max_sd)

# Prepare the Frame
gg_frame <-
  ggplot(data = gois_expression,
         aes(x = Symbol, y = Mean, fill = Symbol)) +
  theme_bw(base_size = 15, base_rect_size = 1.5) +
  theme(axis.text.x = element_text(size = 7, angle = 90,
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
  figure_Name = paste0(GEO_id, "_log2", count_type, "_chart"),
  figure_Folder = out_dir)

# END --------------------------------------------------------------------------

cat("\n", GEO_id, " is done!\n", sep = "")
