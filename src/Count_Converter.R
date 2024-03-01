

# Is it possible to convert read counts to expression values via TPM and return these values?
# https://support.bioconductor.org/p/91218/

# Workshop on RNA-Seq
# https://nbisweden.github.io/workshop-RNAseq/2011/lab_preprocessing.html

# Raw counts to TPM in R
# https://www.biostars.org/p/335187/

# https://search.r-project.org/CRAN/refmans/DGEobj.utils/html/convertCounts.html

# https://www.biostars.org/p/160989/

# https://support.bioconductor.org/p/9146327/

# ---------------------------------------
# https://www.biostars.org/p/307603/
#   https://www.biostars.org/p/83901/
#   https://gist.github.com/slowkow/c6ab0348747f86e2748b
#   https://gist.github.com/slowkow/8101509
# ---------------------------------------

library(GenomicFeatures)
library(clusterProfiler)
library(cmatools)



local_path <- "C:/Users/aleph/Desktop/"

count_file <- paste0(local_path, "GSE76528_mRNA_Counts.txt") 
gtf_file <- paste0(local_path, "GRCh38_GTF.txt")

counts <- read.csv(file = count_file, sep = "\t")
lms(counts)



# 1. Import the GTF-file
txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")


# 2. Collect the exons per gene id
exons_per_gene <- exonsBy(txdb, by = "gene")
# 3. For each gene, reduce all the exons to a set of non overlapping exons,
#    calculate their lengths (widths) and sum then
exonic_gene_sizes <- as.data.frame(sum(width(GenomicRanges::reduce(exons_per_gene))))

lms(exonic_gene_sizes)


genes(exons_per_gene)



bitr(row.names(exonic_gene_sizes), fromType = "E") # no...


