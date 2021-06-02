library(SNPRelate)
library(SeqArray)

gds_file <- commandArgs(TRUE)[1]
gds_out <- sub("\\.gds$", "_pruned.gds", gds_file)

gds <- seqOpen(gds_file)

set.seed(1)
pruned <- snpgdsLDpruning(gds, method = "corr", ld.threshold = sqrt(0.1))

pruned_snps <- unlist(pruned, use.names = FALSE)

seqSetFilter(gds, variant.id = pruned_snps)
seqExport(gds, gds_out)

unlink(gds_file)
