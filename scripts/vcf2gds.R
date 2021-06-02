library(SeqArray)

vcf_file <- commandArgs(TRUE)[1]
gds_file <- commandArgs(TRUE)[2]

#convert
seqVCF2GDS(vcf_file, gds_file, fmt.import = "GT", verbose = FALSE)
