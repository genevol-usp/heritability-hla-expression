library(tidyverse)
library(tidytext)
library(cowplot)
library(ggsci)
library(GGally)

# PCA
plot_parcoord <- function(dataset) {
    ggparcoord(dataset, columns = 3:10,
               groupColumn = "pop", scale = "uniminmax", alpha = .5) +
        scale_color_nejm() +
        theme(panel.grid = element_blank(),
              legend.position = "none") +
        labs(x = NULL, y = "")
}

plot_pca <- function(dataset, pcx, pcy) {
    ggplot(dataset, aes_string(pcx, pcy, color = "pop")) +
        geom_point(alpha = .5) +
        scale_color_nejm() +
        scale_x_continuous(breaks = scales::pretty_breaks(3)) +
        theme(panel.grid.minor = element_blank(),
              legend.position = "none")
}

get_pca_legend <- function(dataset) {
    pl <- ggplot(dataset, aes(V1, V2, color = pop)) +
        geom_point() +
        scale_color_nejm() +
        theme(legend.position = "top") +
        labs(color = "Population: ")
    
    get_legend(pl)
}


annot_geuv <- read_tsv("../data/sample_annotation_geuvadis_tidy.tsv") %>%
    select(sampleid, pop)

pca_1 <- readRDS("../data/pca_dataset1.rds")
pca_2 <- readRDS("../data/pca_dataset2.rds")
pca_3 <- readRDS("../data/pca_dataset3.rds")

pca_df_1 <- as.data.frame(pca_1$eigenvect) %>%
    as_tibble() %>%
    add_column(sampleid = pca_1$sample.id, .before = 1) %>%
    inner_join(annot_geuv, ., by = "sampleid")

pca_df_2 <- as.data.frame(pca_2$eigenvect) %>%
    as_tibble() %>%
    add_column(sampleid = pca_2$sample.id, .before = 1) %>%
    inner_join(annot_geuv, ., by = "sampleid")

pca_df_3 <- as.data.frame(pca_3$eigenvect) %>%
    as_tibble() %>%
    add_column(sampleid = pca_3$sample.id, .before = 1) %>%
    inner_join(annot_geuv, ., by = "sampleid")
          

# dataset 1
leg1 <- get_pca_legend(pca_df_1)
p1_1 <- plot_pca(pca_df_1, "V1", "V2")
p1_2 <- plot_pca(pca_df_1, "V2", "V3")
p1_3 <- plot_pca(pca_df_1, "V3", "V4")
p1_4 <- plot_pca(pca_df_1, "V4", "V5")

pcoord1 <- plot_parcoord(pca_df_1)

grid1 <- 
    plot_grid(NULL,
              plot_grid(p1_1, p1_2, p1_3, p1_4, nrow = 2),
              NULL,
              pcoord1,
              ncol = 1, rel_heights = c(.1, 1, .05, .6)) +
    draw_label("Dataset 1", x = 0.05, y = 0.98, hjust = 0)

# dataset 2
p2_1 <- plot_pca(pca_df_2, "V1", "V2")
p2_2 <- plot_pca(pca_df_2, "V2", "V3")
p2_3 <- plot_pca(pca_df_2, "V3", "V4")
p2_4 <- plot_pca(pca_df_2, "V4", "V5")

pcoord2 <- plot_parcoord(pca_df_2)

grid2 <- 
    plot_grid(NULL,
              plot_grid(p2_1, p2_2, p2_3, p2_4, nrow = 2),
              NULL,
              pcoord2,
              ncol = 1, rel_heights = c(.1, 1, .05, .6)) +
    draw_label("Dataset 2", x = 0.05, y = 0.98, hjust = 0)

# dataset 3
p3_1 <- plot_pca(pca_df_3, "V1", "V2")
p3_2 <- plot_pca(pca_df_3, "V2", "V3")
p3_3 <- plot_pca(pca_df_3, "V3", "V4")
p3_4 <- plot_pca(pca_df_3, "V4", "V5")

pcoord3 <- plot_parcoord(pca_df_3)

grid3 <- 
    plot_grid(NULL,
              plot_grid(p3_1, p3_2, p3_3, p3_4, nrow = 2),
              NULL,
              pcoord3,
              ncol = 1, rel_heights = c(.1, 1, .05, .6)) +
    draw_label("Dataset 3", x = 0.05, y = 0.98, hjust = 0)

plot_grid(leg1, 
          plot_grid(grid1, NULL, grid2, NULL, grid3, 
                    nrow = 1, rel_widths = c(1, .05, 1, .05, 1)),
          rel_heights = c(.1, 1), ncol = 1)

ggsave("pca.png", width = 12, height = 8)


# GRM
grm_gcta_1 <- readRDS("../data/grm_gcta_dataset1.rds")
grm_gcta_2 <- readRDS("../data/grm_gcta_dataset2.rds")
grm_gcta_3 <- readRDS("../data/grm_gcta_dataset3.rds")

grm_wg_1 <- readRDS("../data/grm_wg_dataset1.rds")
grm_wg_2 <- readRDS("../data/grm_wg_dataset2.rds")
grm_wg_3 <- readRDS("../data/grm_wg_dataset3.rds")

# Diagonal values
# Dataset 1
grm_gcta_1_df <- grm_gcta_1 %>%
    as.data.frame() %>%
    rownames_to_column("id1") %>%
    pivot_longer(-id1, names_to = "id2") %>%
    filter(id1 %in% annot_geuv$sampleid & id2 %in% annot_geuv$sampleid) %>%
    left_join(annot_geuv, by = c("id1" = "sampleid")) %>%
    left_join(annot_geuv, by = c("id2" = "sampleid")) %>%
    select(id1, id2, pop1 = pop.x, pop2 = pop.y, value)

grm_wg_1_df <- grm_wg_1 %>%
    as.data.frame() %>%
    rownames_to_column("id1") %>%
    pivot_longer(-id1, names_to = "id2") %>%
    filter(id1 %in% annot_geuv$sampleid & id2 %in% annot_geuv$sampleid) %>%
    left_join(annot_geuv, by = c("id1" = "sampleid")) %>%
    left_join(annot_geuv, by = c("id2" = "sampleid")) %>%
    select(id1, id2, pop1 = pop.x, pop2 = pop.y, value)

grm_1_df <- bind_rows("GCTA" = grm_gcta_1_df,
                    "Weir & Goudet" = grm_wg_1_df,
                    .id = "method")

# Dataset 2
grm_gcta_2_df <- grm_gcta_2 %>%
    as.data.frame() %>%
    rownames_to_column("id1") %>%
    pivot_longer(-id1, names_to = "id2") %>%
    left_join(annot_geuv, by = c("id1" = "sampleid")) %>%
    left_join(annot_geuv, by = c("id2" = "sampleid")) %>%
    select(id1, id2, pop1 = pop.x, pop2 = pop.y, value)

grm_wg_2_df <- grm_wg_2 %>%
    as.data.frame() %>%
    rownames_to_column("id1") %>%
    pivot_longer(-id1, names_to = "id2") %>%
    left_join(annot_geuv, by = c("id1" = "sampleid")) %>%
    left_join(annot_geuv, by = c("id2" = "sampleid")) %>%
    select(id1, id2, pop1 = pop.x, pop2 = pop.y, value)

grm_2_df <- bind_rows("GCTA" = grm_gcta_2_df,
                      "Weir & Goudet" = grm_wg_2_df,
                      .id = "method")

# Dataset 3
grm_gcta_3_df <- grm_gcta_3 %>%
    as.data.frame() %>%
    rownames_to_column("id1") %>%
    pivot_longer(-id1, names_to = "id2") %>%
    left_join(annot_geuv, by = c("id1" = "sampleid")) %>%
    left_join(annot_geuv, by = c("id2" = "sampleid")) %>%
    select(id1, id2, pop1 = pop.x, pop2 = pop.y, value)

grm_wg_3_df <- grm_wg_3 %>%
    as.data.frame() %>%
    rownames_to_column("id1") %>%
    pivot_longer(-id1, names_to = "id2") %>%
    left_join(annot_geuv, by = c("id1" = "sampleid")) %>%
    left_join(annot_geuv, by = c("id2" = "sampleid")) %>%
    select(id1, id2, pop1 = pop.x, pop2 = pop.y, value)

grm_3_df <- bind_rows("GCTA" = grm_gcta_3_df,
                      "Weir & Goudet" = grm_wg_3_df,
                      .id = "method")

grm_plot1 <- grm_1_df %>%
    filter(id1 == id2) %>%
    ggplot(aes(reorder_within(id1, value, method), value)) +
    geom_bar(aes(fill = pop1), stat = "identity", width = 1) +
    facet_wrap(~method, scales = "free_x", ncol = 2) +
    scale_fill_nejm() +
    geom_hline(yintercept = 1L, linetype = 2, size = 1, color = "grey") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_blank(),
          legend.position = "top") +
    labs(fill = "Population:", title = "Dataset 1")

grm_plot2 <- grm_2_df %>%
    filter(id1 == id2) %>%
    ggplot(aes(reorder_within(id1, value, method), value)) +
    geom_bar(aes(fill = pop1), stat = "identity", width = 1) +
    facet_wrap(~method, scales = "free_x", ncol = 2) +
    scale_fill_nejm() +
    geom_hline(yintercept = 1L, linetype = 2, size = 1, color = "grey") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_blank()) +
    labs(fill = "Population:", title = "Dataset 2")

grm_plot3 <- grm_3_df %>%
    filter(id1 == id2) %>%
    ggplot(aes(reorder_within(id1, value, method), value)) +
    geom_bar(aes(fill = pop1), stat = "identity", width = 1) +
    facet_wrap(~method, scales = "free_x", ncol = 2) +
    scale_fill_nejm() +
    geom_hline(yintercept = 1L, linetype = 2, size = 1, color = "grey") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_blank()) +
    labs(fill = "Population:", title = "Dataset 3")

plot_grid(get_legend(grm_plot1),
          grm_plot1 + theme(legend.position = "none"), 
          grm_plot2 + theme(legend.position = "none"), 
          grm_plot3 + theme(legend.position = "none"), 
          ncol = 1, rel_heights = c(.1, 1, 1, 1))

ggsave("./grm_diag.png", height = 6)

# Off-diagonal

plot_off_diag <- function(df, m) {
 
    grm_offdiag <- df %>%
        filter(method == m) %>%
        mutate(value = ifelse(id1 == id2, NA, value)) %>%
        arrange(dataset, method, pop1, pop2, id1, id2) %>%
        mutate(id1 = fct_inorder(id1),
               id2 = fct_inorder(id2))
    
    lines_v <- grm_offdiag %>%
        group_by(dataset, method, pop1) %>%
        filter(1:n() == last(1:n())) %>%
        ungroup() %>%
        select(dataset, method, id1, pop1)
    
    lines_h <- grm_offdiag %>%
        group_by(dataset, method, pop2) %>%
        filter(1:n() == last(1:n())) %>%
        ungroup() %>%
        select(dataset, method, id2, pop2)
    
    labels_v <- grm_offdiag %>%
        group_by(dataset, method, pop1) %>%
        filter(1:n() == floor(median(1:n()))) %>%
        ungroup() %>%
        select(dataset, method, id1, pop1)
    
    labels_h <- grm_offdiag %>%
        group_by(dataset, method, pop2) %>%
        filter(1:n() == floor(median(1:n()))) %>%
        ungroup() %>%
        select(dataset, method, id2, pop2)


    ggplot(grm_offdiag) +
    geom_tile(aes(id1, id2, fill = value)) +
    scale_fill_gradient2(name = "GRM values", 
                         breaks = scales::pretty_breaks(5)) +
    facet_wrap(method~dataset, scales = "free") +
    geom_hline(data = lines_h, aes(yintercept = id2)) +
    geom_vline(data = lines_v, aes(xintercept = id1)) +
    coord_cartesian(clip = "off") +
    geom_text(data = labels_v, aes(x = id1, y = 0, label = pop1), 
              vjust = 1.5) +
    geom_text(data = labels_h, aes(y = id2, x = 0, label = pop2), 
              vjust = -0.5, angle = 90) +
    theme_bw() +
    theme(panel.grid = element_blank(),
          panel.spacing = unit(1, "lines"),
          axis.text = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.margin = unit(c(1, 0, 1, 1), "lines"))
}

grm_all_df <- bind_rows(dataset1 = grm_1_df, 
                        dataset2 = grm_2_df, 
                        dataset3 = grm_3_df,
                        .id = "dataset")

plot_off_1 <- plot_off_diag(grm_all_df, "GCTA")
plot_off_2 <- plot_off_diag(grm_all_df, "Weir & Goudet")

plot_grid(plot_off_1, plot_off_2, ncol = 1)
ggsave("./grm_off.png", width = 10, height = 5)

