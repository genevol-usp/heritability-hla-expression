library(SeqArray)
library(SNPRelate)

args <- commandArgs(TRUE)
gds_file <- args[1]
sample_ids <- args[2]
out <- args[3]

gds <- seqOpen(gds_file)
ids <- readLines(sample_ids)

pca <- snpgdsPCA(gds, sample.id = ids, num.thread = 16L)

saveRDS(pca, sprintf("../data/pca_%s.rds", out))
