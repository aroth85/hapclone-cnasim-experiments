library(signals)
library(ggplot2)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)

print('reading data')
hmmcopy <- fread(file = args[1], sep = "\t", header = TRUE)#, data.table=FALSE)
hmmcopy$chr <- as.character(hmmcopy$chr)
hmmcopy$cell_id <- as.character(hmmcopy$cell_id)

signals <- fread(file = args[2], sep = "\t") # data.table=FALSE)
signals$chr <- as.character(signals$chr)
signals$cell_id <- as.character(signals$cell_id)

print('clustering Signals')
cluster <- umap_clustering(signals, hscn = TRUE)
write.table(cluster$clustering, file=args[3], sep='\t')
#plot <- cluster$clustering %>% filter(clone_id != 0)

ggplot(cluster$clustering, aes(x=umap1, y=umap2, colour=clone_id)) +
  # to create a scatterplot
  geom_point(size=0.5) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave(args[4])

cluster <- umap_clustering(hmmcopy)
write.table(cluster$clustering, file=args[5], sep='\t')
ggplot(cluster$clustering, aes(x=umap1, y=umap2, colour=clone_id)) +
  # to create a scatterplot
  geom_point(size=0.5) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank())

ggsave(args[6])
