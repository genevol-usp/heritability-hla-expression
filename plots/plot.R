library(tidyverse)
library(tidytext)
library(cowplot)
library(ggsci)
library(GGally)
library(ggbeeswarm)

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


##### plots for presentations
# p_pres_1 <- grm_2_df %>%
#     filter(id1 == id2) %>%
#     ggplot(aes(reorder_within(id1, value, method), value)) +
#     geom_bar(aes(fill = pop1), stat = "identity", width = 1) +
#     facet_wrap(~method, scales = "free_x", ncol = 2) +
#     scale_fill_nejm() +
#     geom_hline(yintercept = 1L, linetype = 2, size = 1, color = "grey") +
#     theme_bw() +
#     theme(panel.grid = element_blank(),
#           axis.text.x = element_blank(),
#           axis.line.x = element_blank(),
#           axis.ticks.x = element_blank(),
#           axis.title = element_blank()) +
#     labs(fill = "Population:", title = "Diagonal values")
# 
# plot_off_diag <- function(df, m) {
#     
#     grm_offdiag <- df %>%
#         mutate(value = ifelse(id1 == id2, NA, value)) %>%
#         arrange(method, pop1, pop2, id1, id2) %>%
#         mutate(id1 = fct_inorder(id1),
#                id2 = fct_inorder(id2))
#     
#     lines_v <- grm_offdiag %>%
#         group_by(method, pop1) %>%
#         filter(1:n() == last(1:n())) %>%
#         ungroup() %>%
#         select(method, id1, pop1)
#     
#     lines_h <- grm_offdiag %>%
#         group_by(method, pop2) %>%
#         filter(1:n() == last(1:n())) %>%
#         ungroup() %>%
#         select(method, id2, pop2)
#     
#     labels_v <- grm_offdiag %>%
#         group_by(method, pop1) %>%
#         filter(1:n() == floor(median(1:n()))) %>%
#         ungroup() %>%
#         select(method, id1, pop1)
#     
#     labels_h <- grm_offdiag %>%
#         group_by(method, pop2) %>%
#         filter(1:n() == floor(median(1:n()))) %>%
#         ungroup() %>%
#         select(method, id2, pop2)
#     
#     ggplot(grm_offdiag) +
#         geom_tile(aes(id1, id2, fill = value)) +
#         scale_fill_gradient2(name = "GRM values", 
#                              breaks = scales::pretty_breaks(5)) +
#         facet_wrap(~method, scales = "free") +
#         geom_hline(data = lines_h, aes(yintercept = id2)) +
#         geom_vline(data = lines_v, aes(xintercept = id1)) +
#         coord_cartesian(clip = "off") +
#         geom_text(data = labels_v, aes(x = id1, y = 0, label = pop1), 
#                   vjust = 1.5, size = 2.5) +
#         geom_text(data = labels_h, aes(y = id2, x = 0, label = pop2), 
#                   vjust = -0.5, angle = 90, size = 2.5) +
#         theme_bw() +
#         theme(panel.grid = element_blank(),
#               panel.spacing = unit(1, "lines"),
#               axis.text = element_blank(),
#               axis.line = element_blank(),
#               axis.ticks = element_blank(),
#               axis.title = element_blank(),
#               plot.margin = unit(c(1, 0, 1, 1), "lines")) +
#         labs(title = "Off-diagonal values")
# }
# 
# p_pres_2 <- plot_off_diag(grm_2_df)
# 
# 
# plot_grid(p_pres_1, p_pres_2, ncol = 1)
# 
# ggsave("./pres_plot.png", width = 6, height = 4)
# 


# HWE test

hwe_df <- sprintf("/scratch/vitor/chr%s_geuv_eur.hwe", 1:22) %>%
  map_df(. %>% read_tsv() %>%
           select(CHR, POS, P_HWE))

hwe_cumul_pos <- hwe_df %>%
  group_by(CHR) %>%
  summarise(max_bp = max(POS)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(CHR, bp_add) %>%
  ungroup()

hwe_data <- hwe_df %>% 
  inner_join(hwe_cumul_pos, by = "CHR") %>% 
  mutate(pos_cumul = POS + bp_add)

axis_set <- hwe_data %>% 
  group_by(CHR) %>% 
  summarize(center = median(pos_cumul))

hwe_data_frac <- hwe_data %>%
  group_by(CHR) %>%
  sample_frac(., .25) %>%
  ungroup() %>%
  arrange(CHR, POS) %>%
  mutate(CHR = factor(CHR, levels = 1:22))


manhattan <- 
  ggplot(hwe_data_frac, aes(x = pos_cumul, y = -log10(P_HWE), 
                       color = as_factor(CHR))) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_hline(yintercept = -log10(0.001), linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05/nrow(hwe_data)), linetype = "dashed",
             color = "tomato3") + 
  scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) +
  scale_color_manual(values = rep(c("grey", "cornflowerblue"), nrow(axis_set))) +
  labs(x = NULL, 
       y = expression(paste("-log"[10], italic(P)))) + 
  theme_minimal() +
  theme( 
    legend.position = "none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5)
  )

ggsave("./hwe_geuv_eur.png", manhattan, width = 10, height = 5)


# Expression

hla_expression <- 
  "/raid/genevol/heritability/data/geuvadis_expression.bed.gz" %>%
  read_tsv() %>%
  filter(id == "HLA-A") %>%
  pivot_longer(-(1:6), names_to = "sampleid") %>%
  select(sampleid, value)

hla_expression_std <- 
  "/raid/genevol/heritability/data/geuvadis_expression_std.bed" %>%
  read_tsv() %>%
  filter(id == "HLA-A") %>%
  pivot_longer(-(1:6), names_to = "sampleid") %>%
  select(sampleid, value)

exp1 <- ggplot(hla_expression, aes(value)) +
  geom_histogram() +
  labs(x = NULL, title = "raw")

exp2 <- ggplot(hla_expression_std, aes(value)) +
  geom_histogram() +
  labs(x = NULL, title = "standardized")

exp3 <- left_join(hla_expression, hla_expression_std, 
                by = "sampleid", suffix = c("_raw", "_std")) %>%
  ggplot(aes(value_raw, value_std)) +
  geom_point() +
  labs(x = "Raw", y = "Standardized")

exp_grid_1 <- plot_grid(exp1, exp2, exp3, nrow = 1)

sample_info <- "../data/sample_annotation_geuvadis_tidy.tsv" %>%
  read_tsv() %>%
  select(sampleid, population = pop, sex, lab)

expression_df <- 
  left_join(hla_expression_std, sample_info, by = "sampleid") %>%
  select(sampleid, population, sex, lab, value)

exp_r2_1 <- ggplot(expression_df, aes(reorder(population, value), value)) +
  geom_quasirandom(aes(color = population), method = "smiley") +
  geom_boxplot(fill = NA, outlier.color = NA) +
  scale_color_nejm() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Population", y = "Std expression")

exp_r2_2 <- ggplot(expression_df, aes(reorder(lab, value), value)) +
  geom_quasirandom(color = "grey", method = "smiley") +
  geom_boxplot(fill = NA, outlier.color = NA) +
  scale_x_discrete(labels = function(x) sub("^([^_]+).*$", "\\1", x)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 8)) +
  labs(x = "Lab", y = " ")

exp_r2_3 <- ggplot(expression_df, aes(reorder(sex, value), value)) +
  geom_quasirandom(aes(color = sex), method = "smiley") +
  geom_boxplot(fill = NA, outlier.color = NA) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Sex", y = " ")

exp_grid_2 <- plot_grid(exp_r2_1, exp_r2_2, exp_r2_3, nrow = 1)

exp_grid <- plot_grid(exp_grid_1, NULL, exp_grid_2, 
                      nrow = 3,
                      rel_heights = c(.9, .05, 1), 
                      labels = c("A)", "", "B)"))

ggsave("expression_eda.png", exp_grid, width = 10, height = 5)

