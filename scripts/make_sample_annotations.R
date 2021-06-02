library(tidyverse)
library(readxl)

annot_kgp <- 
    "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" %>%
    read_tsv() %>%
    select(sampleid = sample, pop, continental_pop = super_pop, sex = gender)

annot_geuvadis <- 
    "http://www.ebi.ac.uk/arrayexpress/files/E-GEUV-1/E-GEUV-1.sdrf.txt" %>%
    read_tsv() %>%
    select(sampleid = `Source Name`, lab = Performer) %>%
    distinct() %>%
    inner_join(annot_kgp, by = "sampleid")

annot_geuvadis %>%
    pull(sampleid) %>%
    write_lines("../data/ids_geuvadis.txt")

annot_geuvadis %>%
    filter(pop != "YRI") %>%
    pull(sampleid) %>%
    write_lines("../data/ids_geuvadis_EUR.txt")

write_tsv(annot_kgp, "../data/sample_annotation_1000g_tidy.tsv")
write_tsv(annot_geuvadis, "../data/sample_annotation_geuvadis_tidy.tsv")
