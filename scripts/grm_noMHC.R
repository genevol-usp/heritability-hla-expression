library(SeqArray)
library(SNPRelate)

gds_file <- "../data/allchrs_geuv_eur.gds"
out <- "dataset3_noMHC"
    
gds <- seqOpen(gds_file)
seqSetFilterChrom(gds, include = c(1:5, 7:22), verbose=TRUE)
variant_id <- seqGetData(gds, "variant.id")

grm_gcta_obj <- snpgdsGRM(gds, snp.id = variant_id,
                          method = "GCTA", num.thread = 16L)

grm_gcta <- grm_gcta_obj$grm
rownames(grm_gcta) <- grm_gcta_obj$sample.id
colnames(grm_gcta) <- grm_gcta_obj$sample.id

saveRDS(grm_gcta, sprintf("../data/grm_gcta_%s.rds", out))
