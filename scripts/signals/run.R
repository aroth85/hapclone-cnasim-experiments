library(data.table)
library(ggplot2)
library(signals)

args <- commandArgs(trailingOnly = TRUE)

print('reading data')
print(args[1])
cnbins <- fread(file = args[1], sep = "\t", header = TRUE, data.table = FALSE)
cnbins$chr <- as.character(cnbins$chr)
cnbins$cell_id <- as.character(cnbins$cell_id)


print(args[2])
haplotypes <- fread(file = args[2], sep = "\t", data.table = FALSE)
colnames(haplotypes) <- c(
  "chr",
  "start",
  "end",
  "cell_id",
  "hap_label",
  "allele1",
  "allele0",
  "totalcounts"
)
haplotypes$chr <- as.character(haplotypes$chr)
haplotypes$cell_id <- as.character(haplotypes$cell_id)

print('starting hscn')
hscn <- callHaplotypeSpecificCN(
  cnbins,
  haplotypes,
  firstpassfiltering = FALSE,
  mincells = 1,
  progressbar = FALSE
)
fwrite(hscn$data, file = args[3], sep = '\t')
