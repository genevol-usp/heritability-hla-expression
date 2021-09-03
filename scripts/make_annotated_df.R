library(tidyverse)
library(GENESIS)
library(Biobase)

hla_expression_std <- 
    "/raid/genevol/heritability/data/geuvadis_expression_std.bed" %>%
    read_tsv() %>%
    filter(id == "HLA-A") %>%
    pivot_longer(-(1:6), names_to = "sampleid") %>%
    select(sampleid, value)

sample_info <- "../data/sample_annotation_geuvadis_tidy.tsv" %>%
    read_tsv() %>%
    select(sampleid, population = pop, sex, lab)

expression_df <- 
    left_join(hla_expression_std, sample_info, by = "sampleid") %>%
    select(sampleid, population, sex, lab, value)

pca <- readRDS("../data/pca_dataset3.rds")

pca_df <- as.data.frame(pca$eigenvect) %>%
    as_tibble() %>%
    add_column(sampleid = pca$sample.id, .before = 1) %>%
    select(1:3)


metadata <- 
    c("sample identifier",
      "subject identifier",
      "laboratory of RNA sequencing",
      "PC 1",
      "PC 2",
      "standardized expression levels") %>%
    data.frame(labelDescription = .) 

annotphen <- expression_df %>%
    filter(population != "YRI") %>%
    mutate(sample.id = sampleid) %>%
    left_join(select(pca_df, sampleid, V1:V2), by = "sampleid") %>%
    select(sample.id, sampleid, lab, V1, V2, value) %>%
    as.data.frame() %>%
    AnnotatedDataFrame(metadata)

varMetadata(annotphen)

saveRDS(annotphen, "../data/annotdf.rds")
