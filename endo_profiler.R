

# Package Loading
library(ggplot2)
library(r4tcpl)

# --- Auxiliary Functions ------------------------------------------------------

# Custom log and log^(-1) functions
unlog <- function(x){round((2^x)-1)}
relog <- function(x){log2(x+1)}


# Other required packages: c(stringr, clusterProfiler, org.Hs.eg.db)

# Dirnames and basenames -------------------------------------------------------

local_path <- "C:/Users/FeAR/Desktop/New folder"

count_file <- "Lines/hCMEC D3/GSE76528_TPM.tsv"
count_file <- "Lines/hCMEC D3/GSE138309_TPM.tsv"
count_file |> basename2() |> {\(x)strsplit(x,"_")[[1]][1]}() -> GEO_id

gois_file <- "Analysis/ICT_set.csv"

# ------------------------------------------------------------------------------

# Load data
file.path(local_path, count_file) |> read.delim() -> ncounts
file.path(local_path, gois_file) |> read.delim(header = FALSE) -> gois

file.path(local_path, "RSEM_out") |> read.delim() -> ncounts

# NOTE: For the full channelome use cmatools' TGS: gois <- TGS$ICs

# Check Point
lms(ncounts)
lms(gois)

# Column subsetting ------------------------------------------------------------



gene_ids <- ncounts$gene_id
ncounts <- ncounts[,grep("_TPM$", colnames(ncounts))]


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



