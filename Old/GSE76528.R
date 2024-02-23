


library(clusterProfiler)
library(cmatools)

local_path <- "C:/Users/aleph/Desktop/"
count_file <- paste0(local_path, "GSE76528_mRNA_Counts.txt")
gois_file <- paste0(local_path, "lista_geni_di_mio_interesse.csv")

counts <- read.csv(file = count_file, sep = "\t")
gois <- read.csv(file = cois_file)[,1]
# Full Channelome
#gois <- TGS$ICs

lms(counts)
colnames(counts)

# Use GeneSymbols a rownames (if they are unique)
if (length(dup_report(counts[,2])) == 0) {
  row.names(counts) <- counts[,2]
}
# Keep just the columns of interest
counts <- counts[c(13:20)]

lms(counts)

boxplot(counts)

average_counts <- rowMeans(counts, na.rm = TRUE)

lms(average_counts)


# Subset
IC_expression <- average_counts[names(average_counts) %in% gois]


barplot(IC_expression, xlab = "Ion Channels", ylab = "Normalized Expression")

  
  
  
  
  