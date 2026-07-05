library("tidyverse")
library("adegenet")
library("hierfstat")
library("ade4")
library("ggpubr")
library("naniar")

# constants

COLOR_LIMITS <- c("CEU", "GBR", "FIN", "TSI", "IBS", "MR",
                  "early_neolithic_farmers", "late_neolithic_farmers", "sakhtysh")
COLOR_LABELS <- c("CEU", "GBR", "FIN", "TSI", "IBS", "MR", "EF", "LF", "Sakhtysh")
COLOR_VALUES <- c("#C8A47A", "#A67C52", "#8A5E35", "#6F4424", "#4E2D14", "#2C1608",
                  "#0072B2", "#CC79A7", "#E69F00")

# functions

plot_pcoa <- function(fst_mat, subtitle,
                      flip_ax1 = FALSE, flip_ax2 = FALSE,
                      color_limits = COLOR_LIMITS,
                      color_labels = COLOR_LABELS,
                      color_values = COLOR_VALUES) {
  pco <- dudi.pco(quasieuclid(as.dist(fst_mat)), scannf = FALSE, full = TRUE)
  var_exp <- round(100 * pco$eig / sum(pco$eig), 1)
  df <- pco$tab
  df$Group <- rownames(df)
  if (flip_ax1) df$A1 <- -df$A1
  if (flip_ax2) df$A2 <- -df$A2
  ggplot(df, aes(x = A1, y = A2)) +
    geom_point(aes(color = Group), size = 3, alpha = 0.7) +
    labs(subtitle = subtitle, color = "") +
    xlab(paste0("Ax1 (", var_exp[1], "%)")) +
    ylab(paste0("Ax2 (", var_exp[2], "%)")) +
    scale_color_manual(limits = color_limits,
                       labels = color_labels,
                       values = color_values) +
    theme_light(base_size = 11) +
    theme(
      legend.text = element_text(size = 9, margin = margin(l = 0.5, t = -3)),
      legend.position = "bottom",
      legend.key.width = unit(5, "mm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    guides(colour = guide_legend(nrow = 1))
}

run_fisher_allele_test <- function(data, group_col, col1, col2,
                                   groups = c("Early Volosovo", "Late Volosovo"),
                                   alpha = 0.05) {
  df <- data %>%
    select(all_of(c(group_col, col1, col2))) %>%
    filter(.data[[group_col]] %in% groups) %>%
    drop_na()
  allele_list <- unique(c(df[[col1]], df[[col2]]))
  results <- list()
  for (allele in allele_list) {
    carrier <- as.integer(df[[col1]] == allele | df[[col2]] == allele)
    test <- fisher.test(table(df[[group_col]], carrier))
    if (test$p.value < alpha) {
      results[[allele]] <- list(
        allele = allele,
        p_value = test$p.value,
        table = table(df[[group_col]], carrier)
      )
    }
  }
  results
}

# settings

OUTPUT_DIR <- "results"

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

L <- c("KH180531", "KH180501", "KH180489", "KH180494")

EV <- c("KH180456", "KH180439", "KH180437", "KH180441",
        "KH180444", "KH180453", "KH180421", "KH180452")

LV <- c("KH180416", "KH180423", "KH180428", "KH180429", "KH180450",
        "KH180463", "KH180467", "KH180486", "KH180491", "KH180504",
        "KH180507", "KH180510", "KH180513", "KH180515", "KH180518",
        "KH180521", "KH180523", "KH180530")

# load data

neol_1 <- read.csv("data/combined_HLA_genotypes_russia_class1.csv",
                   sep = ",", na.strings = c("NA", "")) %>%
  filter(Population != "derenburg", group != "Kiel")

neol_1$HLA_A <- paste0(neol_1$HLA_A_1, "-", neol_1$HLA_A_2)
neol_1$HLA_B <- paste0(neol_1$HLA_B_1, "-", neol_1$HLA_B_2)
neol_1$HLA_C <- paste0(neol_1$HLA_C_1, "-", neol_1$HLA_C_2)

neol_1$HLA_A_r <- paste0(neol_1$HLA_A_1_twodigit, "-", neol_1$HLA_A_2_twodigit)
neol_1$HLA_B_r <- paste0(neol_1$HLA_B_1_twodigit, "-", neol_1$HLA_B_2_twodigit)
neol_1$HLA_C_r <- paste0(neol_1$HLA_C_1_twodigit, "-", neol_1$HLA_C_2_twodigit)

neol_2 <- read.csv("data/combined_HLA_genotypes_russia_class2.csv",
                   sep = ",", na.strings = c("NA", "")) %>%
  filter(group != "Kiel", Population != "derenburg")

neol_1_sub <- neol_1 %>%
  mutate(group = case_when(
    Sample_ID %in% L ~ "Lyalovo",
    Sample_ID %in% EV ~ "Early Volosovo",
    Sample_ID %in% LV ~ "Late Volosovo",
    .default = group
  ))

neol_2_sub <- neol_2 %>%
  mutate(group = case_when(
    Sample_ID %in% L ~ "Lyalovo",
    Sample_ID %in% EV ~ "Early Volosovo",
    Sample_ID %in% LV ~ "Late Volosovo",
    .default = group
  ))

neol_DRB <- neol_2 %>%
  select(Sample_ID:group, HLA_DRB1_1_twodigit, HLA_DRB1_2_twodigit) %>%
  drop_na() %>%
  mutate(HLA_DRB = paste0(HLA_DRB1_1_twodigit, "-", HLA_DRB1_2_twodigit)) %>%
  replace_with_na_all(condition = ~.x == "NA") %>%
  replace_with_na_all(condition = ~.x == "NA-NA") %>%
  drop_na() %>%
  select(Sample_ID, Region, Population, group, HLA_DRB)

neol_DQB <- neol_2 %>%
  select(Sample_ID:group, HLA_DQB1_1_twodigit, HLA_DQB1_2_twodigit) %>%
  drop_na() %>%
  mutate(HLA_DQB = paste0(HLA_DQB1_1_twodigit, "-", HLA_DQB1_2_twodigit)) %>%
  replace_with_na_all(condition = ~.x == "NA") %>%
  replace_with_na_all(condition = ~.x == "NA-NA") %>%
  drop_na() %>%
  select(Sample_ID, Region, Population, group, HLA_DQB)

# section 1: pairwise fst and pcoa (pooled sakhtysh)

class1_2f_fst <- neol_1 %>%
  select(Sample_ID, Region, Population, group, HLA_A, HLA_B, HLA_C) %>%
  drop_na()

hf_class1_2f <- genind2hierfstat(df2genind(
  as.data.frame(class1_2f_fst[, c("HLA_A", "HLA_B", "HLA_C"), drop = FALSE]),
  ploidy = 2, ind.names = class1_2f_fst$Sample_ID, pop = class1_2f_fst$group, sep = "-"
))
pair_wc_class1_2f <- pairwise.WCfst(hf_class1_2f)
pair_wc_class1_2f[is.na(pair_wc_class1_2f)] <- 0
pair_wc_class1_2f[pair_wc_class1_2f < 0] <- 0
pco_class1_2f <- plot_pcoa(pair_wc_class1_2f, subtitle = "HLA class-I (second-field)",
                            flip_ax1 = TRUE)

class2_2f_fst <- neol_2 %>%
  select(Sample_ID:group, HLA_DQB1_1:HLA_DRB1_2) %>%
  na.omit() %>%
  mutate(HLA_DRB = paste0(HLA_DRB1_1, "-", HLA_DRB1_2),
         HLA_DQB = paste0(HLA_DQB1_1, "-", HLA_DQB1_2)) %>%
  select(Sample_ID:group, HLA_DRB, HLA_DQB)

hf_class2_2f <- genind2hierfstat(df2genind(
  as.data.frame(class2_2f_fst[, c("HLA_DRB", "HLA_DQB"), drop = FALSE]),
  ploidy = 2, ind.names = class2_2f_fst$Sample_ID, pop = class2_2f_fst$group, sep = "-"
))
pair_wc_class2_2f <- pairwise.WCfst(hf_class2_2f)
pair_wc_class2_2f[is.na(pair_wc_class2_2f)] <- 0
pair_wc_class2_2f[pair_wc_class2_2f < 0] <- 0
pco_class2_2f <- plot_pcoa(pair_wc_class2_2f, subtitle = "HLA class-II (second-field)",
                            flip_ax1 = TRUE)

class1_1f_fst <- neol_1 %>%
  select(Sample_ID, Region, Population, group, HLA_A_r, HLA_B_r, HLA_C_r) %>%
  drop_na()

hf_class1_1f <- genind2hierfstat(df2genind(
  as.data.frame(class1_1f_fst[, c("HLA_A_r", "HLA_B_r", "HLA_C_r"), drop = FALSE]),
  ploidy = 2, ind.names = class1_1f_fst$Sample_ID, pop = class1_1f_fst$group, sep = "-"
))
pair_wc_class1_1f <- pairwise.WCfst(hf_class1_1f)
pair_wc_class1_1f[is.na(pair_wc_class1_1f)] <- 0
pair_wc_class1_1f[pair_wc_class1_1f < 0] <- 0
pco_class1_1f <- plot_pcoa(pair_wc_class1_1f, subtitle = "HLA class-I (first-field)",
                            flip_ax1 = TRUE)

class2_1f_fst <- merge(neol_DRB, neol_DQB)

hf_class2_1f <- genind2hierfstat(df2genind(
  as.data.frame(class2_1f_fst[, c("HLA_DRB", "HLA_DQB"), drop = FALSE]),
  ploidy = 2, ind.names = class2_1f_fst$Sample_ID, pop = class2_1f_fst$group, sep = "-"
))
pair_wc_class2_1f <- pairwise.WCfst(hf_class2_1f)
pair_wc_class2_1f[is.na(pair_wc_class2_1f)] <- 0
pair_wc_class2_1f[pair_wc_class2_1f < 0] <- 0
pco_class2_1f <- plot_pcoa(pair_wc_class2_1f, subtitle = "HLA class-II (first-field)",
                            flip_ax2 = TRUE)

fst_plot <- ggarrange(pco_class1_1f, pco_class2_1f,
                      pco_class1_2f, pco_class2_2f,
                      common.legend = TRUE, legend = "bottom", labels = "AUTO")

ggsave(file.path(OUTPUT_DIR, "fst_pcoa.png"),
       fst_plot, width = 220, height = 120, units = "mm", dpi = 800)
ggsave(file.path(OUTPUT_DIR, "fst_pcoa.pdf"),
       fst_plot, width = 220, height = 120, units = "mm")

# section 2: per-locus subsets with sub-group assignments

neol_A_sub <- neol_1_sub %>%
  select(Sample_ID, Region, Population, group,
         HLA_A_1_twodigit, HLA_A_2_twodigit, HLA_A_r) %>%
  drop_na()

neol_B_sub <- neol_1_sub %>%
  select(Sample_ID, Region, Population, group,
         HLA_B_1_twodigit, HLA_B_2_twodigit, HLA_B_r) %>%
  drop_na()

neol_C_sub <- neol_1_sub %>%
  select(Sample_ID, Region, Population, group,
         HLA_C_1_twodigit, HLA_C_2_twodigit, HLA_C_r) %>%
  drop_na()

neol_A_sub_genind <- df2genind(as.data.frame(neol_A_sub[["HLA_A_r"]]),
  ploidy = 2, ind.names = neol_A_sub$Sample_ID, pop = neol_A_sub$group, sep = "-")
neol_B_sub_genind <- df2genind(as.data.frame(neol_B_sub[["HLA_B_r"]]),
  ploidy = 2, ind.names = neol_B_sub$Sample_ID, pop = neol_B_sub$group, sep = "-")
neol_C_sub_genind <- df2genind(as.data.frame(neol_C_sub[["HLA_C_r"]]),
  ploidy = 2, ind.names = neol_C_sub$Sample_ID, pop = neol_C_sub$group, sep = "-")

neol_A_sub_hf <- genind2hierfstat(neol_A_sub_genind)
neol_B_sub_hf <- genind2hierfstat(neol_B_sub_genind)
neol_C_sub_hf <- genind2hierfstat(neol_C_sub_genind)

neol_DRB_sub <- neol_2_sub %>%
  select(Sample_ID:group, HLA_DRB1_1_twodigit, HLA_DRB1_2_twodigit) %>%
  drop_na() %>%
  mutate(HLA_DRB = paste0(HLA_DRB1_1_twodigit, "-", HLA_DRB1_2_twodigit)) %>%
  replace_with_na_all(condition = ~.x == "NA") %>%
  replace_with_na_all(condition = ~.x == "NA-NA") %>%
  drop_na() %>%
  select(Sample_ID, Region, Population, group, HLA_DRB)

neol_DQB_sub <- neol_2_sub %>%
  select(Sample_ID:group, HLA_DQB1_1_twodigit, HLA_DQB1_2_twodigit) %>%
  drop_na() %>%
  mutate(HLA_DQB = paste0(HLA_DQB1_1_twodigit, "-", HLA_DQB1_2_twodigit)) %>%
  replace_with_na_all(condition = ~.x == "NA") %>%
  replace_with_na_all(condition = ~.x == "NA-NA") %>%
  drop_na() %>%
  select(Sample_ID, Region, Population, group, HLA_DQB)

neol_DRB_sub_genind <- df2genind(as.data.frame(neol_DRB_sub[["HLA_DRB"]]),
  ploidy = 2, ind.names = neol_DRB_sub$Sample_ID, pop = neol_DRB_sub$group, sep = "-")
neol_DQB_sub_genind <- df2genind(as.data.frame(neol_DQB_sub[["HLA_DQB"]]),
  ploidy = 2, ind.names = neol_DQB_sub$Sample_ID, pop = neol_DQB_sub$group, sep = "-")

neol_DRB_sub_hf <- genind2hierfstat(neol_DRB_sub_genind)
neol_DQB_sub_hf <- genind2hierfstat(neol_DQB_sub_genind)

print(allelic.richness(neol_A_sub_hf)$Ar)
print(allelic.richness(neol_B_sub_hf)$Ar)
print(allelic.richness(neol_C_sub_hf)$Ar)
print(allelic.richness(neol_DRB_sub_hf)$Ar)
print(allelic.richness(neol_DQB_sub_hf)$Ar)

basic.stats(neol_A_sub_hf)$Ho
basic.stats(neol_B_sub_hf)$Ho
basic.stats(neol_C_sub_hf)$Ho
basic.stats(neol_DRB_sub_hf)$Ho
basic.stats(neol_DQB_sub_hf)$Ho

basic.stats(neol_A_sub_hf)$n.ind.samp
basic.stats(neol_B_sub_hf)$n.ind.samp
basic.stats(neol_C_sub_hf)$n.ind.samp
basic.stats(neol_DRB_sub_hf)$n.ind.samp
basic.stats(neol_DQB_sub_hf)$n.ind.samp

summary(neol_A_sub_genind)$pop.n.all
summary(neol_B_sub_genind)$pop.n.all
summary(neol_C_sub_genind)$pop.n.all
summary(neol_DRB_sub_genind)$pop.n.all
summary(neol_DQB_sub_genind)$pop.n.all

ev_lv_c1 <- neol_1_sub %>%
  filter(group %in% c("Early Volosovo", "Late Volosovo")) %>%
  select(Sample_ID, Region, Population, group, HLA_A, HLA_B, HLA_C) %>%
  drop_na()

hf_ev_lv_c1 <- genind2hierfstat(df2genind(
  as.data.frame(ev_lv_c1[, c("HLA_A", "HLA_B", "HLA_C"), drop = FALSE]),
  ploidy = 2, ind.names = ev_lv_c1$Sample_ID, pop = ev_lv_c1$group, sep = "-"
))
fst_ev_lv_c1 <- pairwise.WCfst(hf_ev_lv_c1)
fst_ev_lv_c1[is.na(fst_ev_lv_c1)] <- 0
fst_ev_lv_c1[fst_ev_lv_c1 < 0] <- 0
print(fst_ev_lv_c1)

ev_lv_c2 <- neol_2_sub %>%
  filter(group %in% c("Early Volosovo", "Late Volosovo")) %>%
  select(Sample_ID:group, HLA_DQB1_1:HLA_DRB1_2) %>%
  na.omit() %>%
  mutate(HLA_DRB = paste0(HLA_DRB1_1, "-", HLA_DRB1_2),
         HLA_DQB = paste0(HLA_DQB1_1, "-", HLA_DQB1_2)) %>%
  select(Sample_ID:group, HLA_DRB, HLA_DQB)

hf_ev_lv_c2 <- genind2hierfstat(df2genind(
  as.data.frame(ev_lv_c2[, c("HLA_DRB", "HLA_DQB"), drop = FALSE]),
  ploidy = 2, ind.names = ev_lv_c2$Sample_ID, pop = ev_lv_c2$group, sep = "-"
))
fst_ev_lv_c2 <- pairwise.WCfst(hf_ev_lv_c2)
fst_ev_lv_c2[is.na(fst_ev_lv_c2)] <- 0
fst_ev_lv_c2[fst_ev_lv_c2 < 0] <- 0
print(fst_ev_lv_c2)

# section 3: fisher's exact tests (early vs. late volosovo)

fisher_A <- run_fisher_allele_test(neol_1_sub, "group", "HLA_A_1", "HLA_A_2")
for (res in fisher_A) {
  cat("allele:", res$allele, "p =", round(res$p_value, 4), "\n")
  print(res$table)
}

fisher_B <- run_fisher_allele_test(neol_1_sub, "group", "HLA_B_1", "HLA_B_2")
for (res in fisher_B) {
  cat("allele:", res$allele, "p =", round(res$p_value, 4), "\n")
  print(res$table)
}

fisher_C <- run_fisher_allele_test(neol_1_sub, "group", "HLA_C_1", "HLA_C_2")
for (res in fisher_C) {
  cat("allele:", res$allele, "p =", round(res$p_value, 4), "\n")
  print(res$table)
}

fisher_DRB <- run_fisher_allele_test(neol_2_sub, "group", "HLA_DRB1_1", "HLA_DRB1_2")
for (res in fisher_DRB) {
  cat("allele:", res$allele, "p =", round(res$p_value, 4), "\n")
  print(res$table)
}

fisher_DQB <- run_fisher_allele_test(neol_2_sub, "group", "HLA_DQB1_1", "HLA_DQB1_2")
for (res in fisher_DQB) {
  cat("allele:", res$allele, "p =", round(res$p_value, 4), "\n")
  print(res$table)
}

# section 4: lyalovo inclusion vs. exclusion

c1_all <- neol_1 %>%
  select(Sample_ID, Region, Population, group, HLA_A, HLA_B, HLA_C) %>%
  drop_na()

hf_c1_all <- genind2hierfstat(df2genind(
  as.data.frame(c1_all[, c("HLA_A", "HLA_B", "HLA_C"), drop = FALSE]),
  ploidy = 2, ind.names = c1_all$Sample_ID, pop = c1_all$group, sep = "-"
))
fst_c1_all <- pairwise.WCfst(hf_c1_all)
fst_c1_all[is.na(fst_c1_all)] <- 0
fst_c1_all[fst_c1_all < 0] <- 0
pco_c1_all <- plot_pcoa(fst_c1_all, subtitle = "HLA class-I (LYALOVO INCLUDED)")

c2_all <- neol_2 %>%
  select(Sample_ID:group, HLA_DQB1_1:HLA_DRB1_2) %>%
  na.omit() %>%
  mutate(HLA_DRB = paste0(HLA_DRB1_1, "-", HLA_DRB1_2),
         HLA_DQB = paste0(HLA_DQB1_1, "-", HLA_DQB1_2)) %>%
  select(Sample_ID:group, HLA_DRB, HLA_DQB)

hf_c2_all <- genind2hierfstat(df2genind(
  as.data.frame(c2_all[, c("HLA_DRB", "HLA_DQB"), drop = FALSE]),
  ploidy = 2, ind.names = c2_all$Sample_ID, pop = c2_all$group, sep = "-"
))
fst_c2_all <- pairwise.WCfst(hf_c2_all)
fst_c2_all[is.na(fst_c2_all)] <- 0
fst_c2_all[fst_c2_all < 0] <- 0
pco_c2_all <- plot_pcoa(fst_c2_all, subtitle = "HLA class-II (LYALOVO INCLUDED)")

c1_no_lya <- neol_1 %>%
  filter(!Sample_ID %in% L) %>%
  select(Sample_ID, Region, Population, group, HLA_A, HLA_B, HLA_C) %>%
  drop_na()

hf_c1_no_lya <- genind2hierfstat(df2genind(
  as.data.frame(c1_no_lya[, c("HLA_A", "HLA_B", "HLA_C"), drop = FALSE]),
  ploidy = 2, ind.names = c1_no_lya$Sample_ID, pop = c1_no_lya$group, sep = "-"
))
fst_c1_no_lya <- pairwise.WCfst(hf_c1_no_lya)
fst_c1_no_lya[is.na(fst_c1_no_lya)] <- 0
fst_c1_no_lya[fst_c1_no_lya < 0] <- 0
pco_c1_no_lya <- plot_pcoa(fst_c1_no_lya, subtitle = "HLA class-I (LYALOVO EXCLUDED)")

c2_no_lya <- neol_2 %>%
  filter(!Sample_ID %in% L) %>%
  select(Sample_ID:group, HLA_DQB1_1:HLA_DRB1_2) %>%
  na.omit() %>%
  mutate(HLA_DRB = paste0(HLA_DRB1_1, "-", HLA_DRB1_2),
         HLA_DQB = paste0(HLA_DQB1_1, "-", HLA_DQB1_2)) %>%
  select(Sample_ID:group, HLA_DRB, HLA_DQB)

hf_c2_no_lya <- genind2hierfstat(df2genind(
  as.data.frame(c2_no_lya[, c("HLA_DRB", "HLA_DQB"), drop = FALSE]),
  ploidy = 2, ind.names = c2_no_lya$Sample_ID, pop = c2_no_lya$group, sep = "-"
))
fst_c2_no_lya <- pairwise.WCfst(hf_c2_no_lya)
fst_c2_no_lya[is.na(fst_c2_no_lya)] <- 0
fst_c2_no_lya[fst_c2_no_lya < 0] <- 0
pco_c2_no_lya <- plot_pcoa(fst_c2_no_lya, subtitle = "HLA class-II (LYALOVO EXCLUDED)")

lyalovo_comparison <- ggarrange(
  pco_c1_all, pco_c2_all,
  pco_c1_no_lya, pco_c2_no_lya,
  common.legend = TRUE, legend = "bottom",
  labels = "AUTO", nrow = 2, ncol = 2
)

ggsave(file.path(OUTPUT_DIR, "lyalovo_comparison_pcoa.png"),
       lyalovo_comparison, width = 220, height = 120, units = "mm", dpi = 800)
ggsave(file.path(OUTPUT_DIR, "lyalovo_comparison_pcoa.pdf"),
       lyalovo_comparison, width = 220, height = 120, units = "mm")
