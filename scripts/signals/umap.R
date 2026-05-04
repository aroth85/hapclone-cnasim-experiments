library(signals)
library(ggplot2)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

print('reading data')
hmmcopy <- fread(file = args[1], sep = "\t")#, data.table=FALSE)
hmmcopy$chr <- as.character(hmmcopy$chr)
hmmcopy$cell_id <- as.character(hmmcopy$cell_id)

signals <- fread(file = args[2], sep = ",") # data.table=FALSE)
signals$chr <- as.character(signals$chr)
signals$cell_id <- as.character(signals$cell_id)

print('clustering Signals')
cluster <- umap_clustering(signals, hscn = TRUE)
write.table(cluster$clustering, file=args[3], sep='\t')

print('clustering HmmCopy')
cluster <- umap_clustering(hmmcopy)
write.table(cluster$clustering, file=args[4], sep='\t')

