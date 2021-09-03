library(SeqArray)

gds_list_1 <- sprintf("/scratch/vitor/chr%d_pruned.gds", 1:22)
gds_list_2 <- sprintf("/scratch/vitor/chr%d_geuv_pruned.gds", 1:22)
gds_list_3 <- sprintf("/scratch/vitor/chr%d_geuv_eur_pruned.gds", 1:22)

gds_out_1 <- "/home/vitor/heritability-hla-expression/data/allchrs.gds"
gds_out_2 <- "/home/vitor/heritability-hla-expression/data/allchrs_geuv.gds"
gds_out_3 <- "/home/vitor/heritability-hla-expression/data/allchrs_geuv_eur.gds"

seqMerge(gds_list_1, gds_out_1)
seqMerge(gds_list_2, gds_out_2)
seqMerge(gds_list_3, gds_out_3)

unlink(gds_list_1)
unlink(gds_list_2)
unlink(gds_list_3)