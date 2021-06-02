library(SeqArray)
library(SNPRelate)

args <- commandArgs(TRUE)
gds_file <- args[1]
out <- args[2] 
    
gds <- seqOpen(gds_file)

grm_gcta_obj <- snpgdsGRM(gds, method = "GCTA", num.thread = 16L)
grm_wg_obj <- snpgdsGRM(gds, method = "IndivBeta", num.thread = 16L)

grm_gcta <- grm_gcta_obj$grm
rownames(grm_gcta) <- grm_gcta_obj$sample.id
colnames(grm_gcta) <- grm_gcta_obj$sample.id

grm_wg <- grm_wg_obj$grm
rownames(grm_wg) <- grm_wg_obj$sample.id
colnames(grm_wg) <- grm_wg_obj$sample.id

saveRDS(grm_gcta, sprintf("../data/grm_gcta_%s.rds", out))
saveRDS(grm_wg, sprintf("../data/grm_wg_%s.rds", out))
