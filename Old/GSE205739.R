

# Package Loading
library(ggplot2)
library(r4tcpl)

# Other required packages: c(stringr, clusterProfiler, org.Hs.eg.db)

# Dirnames and basenames -------------------------------------------------------
local_path <- "C:/Users/aleph/Desktop/"
local_path <- "C:/Users/FeAR/Desktop/"
count_file <- "GSE205739_FPKM.csv"
count_type <- "FPKM"
GEO_id <- stringr::str_extract(count_file, "GSE\\d{3,10}")
gois_file <- "GOIs.csv"
# ------------------------------------------------------------------------------

# Load data
ncounts <- read.csv(file = paste0(local_path, count_file), sep = ",")
gois <- read.csv(file = paste0(local_path, gois_file))[,1]
# NOTE: For the full channelome use cmatools' TGS: gois <- TGS$ICs

# Check Point
lms(ncounts)
lms(gois)

# Column subsetting ------------------------------------------------------------
colnames(ncounts)
gene_ids <- ncounts[,1]
ncounts <- ncounts[c(2:5)]
# ------------------------------------------------------------------------------
genome_size <- dim(ncounts)[1]
sample_size <- dim(ncounts)[2]
lms(gene_ids)
lms(ncounts)

# Possibly convert to TPM ------------------------------------------------------
colSums(ncounts)
ncounts <- fpkm2tpm(ncounts)
count_type <- "TPM"
lms(ncounts)
# Normalization check (== 10^6)
colSums(ncounts)
# ------------------------------------------------------------------------------



# log2 transformation

are_Log_Counts <- function(count_data)
{
  qx <- as.numeric(quantile(count_data, c(0.99, 1.0), na.rm = TRUE))
  print(qx)
  print(qx[2]-qx[1])
  
  if (qx[2]-qx[1] < 1e3) {
    warning("Data seems already log-transformed!")
  } else {
    print("Linear Data")
  }
}


are_Log_Counts(ncounts)
are_Log_Counts(ncounts[,1])
are_Log_Counts(ncounts[ncounts[,1] > 0, 1])


#' issues:
#'   add a RNA-Seq example data set
#' use it to add @examples in fpkm2tpm and count_density functions
#' 
#' 
#' test quick_chart
#' #' @examples
#' #' # Get a graphical representation



# Check distribution and log-scale transformation
boxplot(ncounts)
count_density(ncounts)
ncounts <- log2(ncounts + 1)
boxplot(ncounts)
count_density(ncounts, remove_zeros = F, ylim = c(0,0.1))
count_density(ncounts)

gmm <- GMM_divide(ncounts[ncounts[,1] != 0, 1], G = 3)
for (i in 1:gmm$fit$G) {
  lines(gmm$x, gmm$components[,i], col = "blue")
}
lines(gmm$x, rowSums(gmm$components), col = "red")
# rug(raw_data)

y_lim <- par("yaxp")[2]
lines(c(gmm$boundary[1], gmm$boundary[1]), c(0, 1.5*y_lim))





# Biological ID Translator to get Gene Symbols
more_ids <- clusterProfiler::bitr(gene_ids,
                                  fromType = "ENSEMBL",
                                  toType = "SYMBOL",
                                  OrgDb = "org.Hs.eg.db",
                                  drop = FALSE)
lms(more_ids)
dnues2(more_ids$ENSEMBL)
dnues2(more_ids$SYMBOL)

# Row subsetting
gois_ids <- more_ids[more_ids$SYMBOL %in% gois, ]

# Check Point
lms(gois_ids)
dnues2(gois_ids$ENSEMBL)
dnues2(gois_ids$SYMBOL)

# Possibly collapse duplicated Gene Symbols (keeping the most expressed)
#dnues2(more_ids[,1])
#DEGs <- DEGs[order(DEGs$adj_pval), ] ## from GOZER; to be adapted
#DEGs <- DEGs[!duplicated(DEGs$GeneSymbol), ]

# Are there any GOIs missing in the Seq?
venny(gois_ids$SYMBOL, gois, lab = c("GOIs in Seq", "GOIs"))$diff_BA

# Merge Symbols with ncounts
gois_ncounts <- merge(gois_ids,
                      data.frame(gene_ids = gene_ids, ncounts),
                      by.x = "ENSEMBL", by.y = "gene_ids",
                      all.x = TRUE)
lms(gois_ncounts)

# Compute descriptive statistics (mean and SD)
average_ncounts <- rowMeans(gois_ncounts[,-c(1,2)], na.rm = TRUE)
sd_ncounts <- apply(gois_ncounts[,-c(1,2)], 1, sd, na.rm = TRUE)
lms(average_ncounts)
lms(sd_ncounts)

# Final expression matrix
gois_expression <- data.frame(gois_ncounts[,c(1,2)], average_ncounts, sd_ncounts)
lms(gois_expression)

# Now plot ---------------------------------------------------------------------

# Colors
line_color <- "gray17"
point_color <- "steelblue4"
fill_col <- "slategray4"
err_color <- "gray17"

y_max <- max(gois_expression$average_ncounts)
y_max_sd <- gois_expression$sd_ncounts[gois_expression$average_ncounts == y_max]
y_limit <- ceiling(y_max + y_max_sd)
y_limit

gg_frame <-
  ggplot(data = gois_expression,
         aes(x = SYMBOL, y = average_ncounts, fill = SYMBOL)) +
  theme_bw(base_size = 15, base_rect_size = 1.5) +
  theme(axis.text.x = element_text(size = 7, angle = 90,
                                   vjust = 0.5, hjust = 1),
        axis.text.y = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.position = "none") +
  scale_y_continuous(expand = c(0.01, 0.02),
                     breaks = seq(0, y_limit, 0.5)) +
  xlab("Genes of Interest") + ylab(substitute(log[2]*(x+1),list(x=count_type))) +
  ggtitle(label = paste0(GEO_id, " (n_ctrl = ", sample_size,")"))

gg_bars <- gg_frame +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = average_ncounts - sd_ncounts,
                    ymax = average_ncounts + sd_ncounts),
                linewidth = 1.0, width = 0.5, color = err_color)

# # Alternative with borders
# gg_bars <- gg_frame +
#   geom_bar(stat = "identity", width = 0.7,
#            color = line_color, linewidth = 0.1) +
#   geom_errorbar(aes(ymin = average_ncounts - sd_ncounts,
#                     ymax = average_ncounts + sd_ncounts),
#                 linewidth = 1.0, width = 0.5, color = err_color)

gg_bars +
  geom_hline(yintercept = 1,
             linetype = "dashed",
             color = line_color,
             linewidth = 1)

# Print Plot
dev.print(device = png, filename = paste0(GEO_id, "_", count_type, ".png"),
          width = 2000, height = 984)



