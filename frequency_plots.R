library("tidyverse")
library("adegenet")
library("hierfstat")
library("ggpubr")
library("naniar")

# constants

MR_LABEL <- "Russia Nizhny Novgorod, Russians"

STUDY_POPS_FREQ <- c("CEU", "GBR", "FIN", "TSI", "IBS", MR_LABEL,
                     "early_neolithic_farmers", "late_neolithic_farmers", "sakhtysh")

FILL_LIMITS <- c("CEU", "GBR", "FIN", "TSI", "IBS",
                 "Russia Nizhny Novgorod, Russians",
                 "early_neolithic_farmers", "late_neolithic_farmers", "sakhtysh")
FILL_LABELS <- c("CEU", "GBR", "FIN", "TSI", "IBS", "MR", "EF", "LF", "Sakhtysh")
FILL_VALUES <- c("#C8A47A", "#A67C52", "#8A5E35", "#6F4424", "#4E2D14", "#2C1608",
                 "#0072B2", "#CC79A7", "#E69F00")

MAINFIGURE_LIMITS <- c("CEU", "Russia Nizhny Novgorod, Russians",
                      "early_neolithic_farmers", "late_neolithic_farmers", "sakhtysh")
MAINFIGURE_LABELS <- c("CEU", "MR", "EF", "LF", "Sakhtysh")
MAINFIGURE_VALUES <- c("#C8A47A", "#2C1608", "#0072B2", "#CC79A7", "#E69F00")

# functions

extract_freq <- function(genind_obj, locus_idx, locus_col, loci_filter,
                         modern_russia_df) {
  freq_mat <- pop.freq(genind2hierfstat(genind_obj))[[locus_idx]]
  rownames(freq_mat) <- sort(genind_obj@all.names[[locus_col]])
  freq_df <- as.data.frame(freq_mat)
  freq_df <- merge(
    as.data.frame(table(genind_obj@pop)),
    freq_df, by.x = "Var1", by.y = "Var2"
  )
  freq_df <- rbind(
    freq_df,
    modern_russia_df %>% filter(loci == loci_filter) %>% select(-loci)
  )
  freq_df %>% mutate(CI = 1.96 * sqrt((Freq.y * (1 - Freq.y)) / Freq.x))
}

plot_freq_bar <- function(freq_df, x_label, top_n = 5, base_size = 12.5) {
  common_alleles <- freq_df %>%
    arrange(desc(Freq.y)) %>%
    group_by(Var1) %>%
    slice_head(n = top_n) %>%
    pull(x) %>%
    unique()
  plot_df <- freq_df %>%
    filter(x %in% common_alleles) %>%
    mutate(Var1 = factor(Var1, levels = FILL_LIMITS))
  ggplot(plot_df, aes(y = Freq.y, x = x, fill = Var1)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(limits = FILL_LIMITS,
                      labels = FILL_LABELS,
                      values = FILL_VALUES) +
    theme_classic(base_size = base_size) +
    theme(axis.text.x = element_text(size = 8, angle = 25, vjust = 0.7)) +
    xlab(x_label) +
    ylab("frequency") +
    labs(fill = "population") +
    guides(fill = guide_legend(nrow = 1, title = ""))
}

plot_mainfigure <- function(combined_freq_df, alleles, x_order = alleles,
                           y_label = "", base_size = 10.5, label_size = 6.5) {
  df <- combined_freq_df %>%
    filter(Var1 %in% MAINFIGURE_LIMITS, x %in% alleles) %>%
    mutate(Var1 = factor(Var1, levels = MAINFIGURE_LIMITS))
  p <- ggplot(df, aes(y = Freq.y, x = x, fill = Var1)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_x_discrete(limits = x_order) +
    scale_fill_manual(limits = MAINFIGURE_LIMITS,
                      labels = MAINFIGURE_LABELS,
                      values = MAINFIGURE_VALUES) +
    theme_classic(base_size = base_size) +
    xlab("") +
    ylab(y_label) +
    labs(fill = "")
  if (!is.null(label_size)) p <- p + theme(axis.text.x = element_text(size = label_size))
  p
}

# settings

OUTPUT_DIR <- "results"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# load data

modern_russia <- read.csv("data/russia_AFND.csv") %>%
  filter(population == MR_LABEL)
colnames(modern_russia) <- c("x", "loci", "Var1", "Freq.y", "Freq.x")

class1_df <- read.csv("data/combined_HLA_genotypes_class1.csv", sep = ",") %>%
  filter(Population != "derenburg") %>%
  replace_with_na_all(condition = ~.x == "") %>%
  drop_na(c(HLA_A_1_twodigit, HLA_B_1_twodigit, HLA_C_1_twodigit,
            HLA_A_2_twodigit, HLA_B_2_twodigit, HLA_C_2_twodigit))

class1_df$HLA_A <- paste0(class1_df$HLA_A_1, "-", class1_df$HLA_A_2)
class1_df$HLA_B <- paste0(class1_df$HLA_B_1, "-", class1_df$HLA_B_2)
class1_df$HLA_C <- paste0(class1_df$HLA_C_1, "-", class1_df$HLA_C_2)

class1_df <- class1_df %>%
  drop_na() %>%
  select(Sample_ID, Region, Population, Time, group, HLA_A, HLA_B, HLA_C) %>%
  mutate(HLA_A_r = HLA_A, HLA_B_r = HLA_B, HLA_C_r = HLA_C)

class2_df <- read.csv("data/combined_HLA_genotypes_class2.csv", sep = ",") %>%
  filter(Population != "derenburg") %>%
  replace_with_na_all(condition = ~.x == "") %>%
  drop_na(c(HLA_DQB1_1_twodigit, HLA_DQB1_2_twodigit,
            HLA_DRB1_1_twodigit, HLA_DRB1_2_twodigit))

class2_df$HLA_DRB <- paste0(class2_df$HLA_DRB1_1, "-", class2_df$HLA_DRB1_2)
class2_df$HLA_DQB <- paste0(class2_df$HLA_DQB1_1, "-", class2_df$HLA_DQB1_2)

class2_df <- class2_df %>%
  drop_na() %>%
  select(Sample_ID, Region, Population, Time, group, HLA_DRB, HLA_DQB) %>%
  mutate(HLA_DRB_r = HLA_DRB, HLA_DQB_r = HLA_DQB)

recode_population <- function(df) {
  df %>% mutate(Population = case_when(
    group == "late_neolithic_farmers" ~ "late_neolithic_farmers",
    group == "early_neolithic_farmers" ~ "early_neolithic_farmers",
    TRUE ~ Population
  ))
}

class1_freq_df <- recode_population(class1_df)

neol_1_genind <- df2genind(
  as.data.frame(select(class1_freq_df, HLA_A_r, HLA_B_r, HLA_C_r)),
  ploidy = 2,
  ind.names = class1_freq_df$Sample_ID,
  pop = class1_freq_df$Population,
  sep = "-"
)

HLA_A_frequencies <- extract_freq(neol_1_genind, "X1", "HLA_A_r", "A", modern_russia)
HLA_B_frequencies <- extract_freq(neol_1_genind, "X2", "HLA_B_r", "B", modern_russia)
HLA_C_frequencies <- extract_freq(neol_1_genind, "X3", "HLA_C_r", "C", modern_russia)

class2_freq_df <- recode_population(class2_df)

neol_2_genind <- df2genind(
  as.data.frame(select(class2_freq_df, HLA_DRB_r, HLA_DQB_r)),
  ploidy = 2,
  ind.names = class2_freq_df$Sample_ID,
  pop = class2_freq_df$Population,
  sep = "-"
)

HLA_DRB_frequencies <- extract_freq(neol_2_genind, "X1", "HLA_DRB_r", "DRB1", modern_russia)
HLA_DQB_frequencies <- extract_freq(neol_2_genind, "X2", "HLA_DQB_r", "DQB1", modern_russia)

# section 1: frequency bar plots

filter_study_pops <- function(freq_df) freq_df %>% filter(Var1 %in% STUDY_POPS_FREQ)

barplotA <- plot_freq_bar(filter_study_pops(HLA_A_frequencies), "HLA-A alleles")
barplotB <- plot_freq_bar(filter_study_pops(HLA_B_frequencies), "HLA-B alleles")
barplotC <- plot_freq_bar(filter_study_pops(HLA_C_frequencies), "HLA-C alleles")
barplotDRB <- plot_freq_bar(filter_study_pops(HLA_DRB_frequencies), "HLA-DRB1 alleles")
barplotDQB <- plot_freq_bar(filter_study_pops(HLA_DQB_frequencies), "HLA-DQB1 alleles")

combined_barplot <- ggarrange(barplotA, barplotB, barplotC, barplotDRB, barplotDQB,
                              common.legend = TRUE, nrow = 5,
                              align = "v", legend = "bottom")

ggsave(file.path(OUTPUT_DIR,"combined_frequencyplot.png"),
       combined_barplot, units = "mm", width = 195, height = 275, dpi = 800)
ggsave(file.path(OUTPUT_DIR,"combined_frequencyplot.pdf"),
       combined_barplot, units = "mm", width = 195, height = 275)

# section 2: first-field frequency bar plots (Fig. S12)

modern_russia_1f <- read.csv("data/russia_AFND_1field.csv") %>%
  filter(population == MR_LABEL)
colnames(modern_russia_1f) <- c("x", "loci", "Var1", "Freq.y", "Freq.x")

class1_1f_raw <- read.csv("data/combined_HLA_genotypes_class1.csv", sep = ",") %>%
  filter(Population != "derenburg") %>%
  replace_with_na_all(condition = ~.x == "") %>%
  drop_na(c(HLA_A_1_twodigit, HLA_B_1_twodigit, HLA_C_1_twodigit,
            HLA_A_2_twodigit, HLA_B_2_twodigit, HLA_C_2_twodigit))

class1_1f_raw$HLA_A <- paste0(class1_1f_raw$HLA_A_1_twodigit, "-", class1_1f_raw$HLA_A_2_twodigit)
class1_1f_raw$HLA_B <- paste0(class1_1f_raw$HLA_B_1_twodigit, "-", class1_1f_raw$HLA_B_2_twodigit)
class1_1f_raw$HLA_C <- paste0(class1_1f_raw$HLA_C_1_twodigit, "-", class1_1f_raw$HLA_C_2_twodigit)

class1_1f_freq_df <- recode_population(
  class1_1f_raw %>%
    select(Sample_ID, Region, Population, Time, group, HLA_A, HLA_B, HLA_C) %>%
    mutate(HLA_A_r = HLA_A, HLA_B_r = HLA_B, HLA_C_r = HLA_C)
)

neol_1_1f_genind <- df2genind(
  as.data.frame(select(class1_1f_freq_df, HLA_A_r, HLA_B_r, HLA_C_r)),
  ploidy = 2,
  ind.names = class1_1f_freq_df$Sample_ID,
  pop = class1_1f_freq_df$Population,
  sep = "-"
)

HLA_A_1f_frequencies <- extract_freq(neol_1_1f_genind, "X1", "HLA_A_r", "A", modern_russia_1f)
HLA_B_1f_frequencies <- extract_freq(neol_1_1f_genind, "X2", "HLA_B_r", "B", modern_russia_1f)
HLA_C_1f_frequencies <- extract_freq(neol_1_1f_genind, "X3", "HLA_C_r", "C", modern_russia_1f)

class2_1f_raw <- read.csv("data/combined_HLA_genotypes_class2.csv", sep = ",") %>%
  filter(Population != "derenburg") %>%
  replace_with_na_all(condition = ~.x == "") %>%
  drop_na(c(HLA_DQB1_1_twodigit, HLA_DQB1_2_twodigit,
            HLA_DRB1_1_twodigit, HLA_DRB1_2_twodigit))

class2_1f_raw$HLA_DRB <- paste0(class2_1f_raw$HLA_DRB1_1_twodigit, "-", class2_1f_raw$HLA_DRB1_2_twodigit)
class2_1f_raw$HLA_DQB <- paste0(class2_1f_raw$HLA_DQB1_1_twodigit, "-", class2_1f_raw$HLA_DQB1_2_twodigit)

class2_1f_freq_df <- recode_population(
  class2_1f_raw %>%
    select(Sample_ID, Region, Population, Time, group, HLA_DRB, HLA_DQB) %>%
    mutate(HLA_DRB_r = HLA_DRB, HLA_DQB_r = HLA_DQB)
)

neol_2_1f_genind <- df2genind(
  as.data.frame(select(class2_1f_freq_df, HLA_DRB_r, HLA_DQB_r)),
  ploidy = 2,
  ind.names = class2_1f_freq_df$Sample_ID,
  pop = class2_1f_freq_df$Population,
  sep = "-"
)

HLA_DRB_1f_frequencies <- extract_freq(neol_2_1f_genind, "X1", "HLA_DRB_r", "DRB1", modern_russia_1f)
HLA_DQB_1f_frequencies <- extract_freq(neol_2_1f_genind, "X2", "HLA_DQB_r", "DQB1", modern_russia_1f)

barplotA_1f <- plot_freq_bar(filter_study_pops(HLA_A_1f_frequencies), "HLA-A alleles")
barplotB_1f <- plot_freq_bar(filter_study_pops(HLA_B_1f_frequencies), "HLA-B alleles")
barplotC_1f <- plot_freq_bar(filter_study_pops(HLA_C_1f_frequencies), "HLA-C alleles")
barplotDRB_1f <- plot_freq_bar(filter_study_pops(HLA_DRB_1f_frequencies), "HLA-DRB1 alleles")
barplotDQB_1f <- plot_freq_bar(filter_study_pops(HLA_DQB_1f_frequencies), "HLA-DQB1 alleles")

combined_barplot_1f <- ggarrange(barplotA_1f, barplotB_1f, barplotC_1f, barplotDRB_1f, barplotDQB_1f,
                                 common.legend = TRUE, nrow = 5,
                                 align = "v", legend = "bottom")

ggsave(file.path(OUTPUT_DIR, "combined_frequencyplot_1field.png"),
       combined_barplot_1f, units = "mm", width = 195, height = 275, dpi = 800)
ggsave(file.path(OUTPUT_DIR, "combined_frequencyplot_1field.pdf"),
       combined_barplot_1f, units = "mm", width = 195, height = 275)

# section 3: spotlight plots

combined_all <- bind_rows(
  mutate(filter_study_pops(HLA_A_frequencies), locus = "A"),
  mutate(filter_study_pops(HLA_B_frequencies), locus = "B"),
  mutate(filter_study_pops(HLA_C_frequencies), locus = "C"),
  mutate(filter_study_pops(HLA_DRB_frequencies), locus = "DRB1"),
  mutate(filter_study_pops(HLA_DQB_frequencies), locus = "DQB1")
)

pA <- plot_mainfigure(combined_all, c("A*01:01", "A*31:01"), y_label = "frequency", base_size = 11.5, label_size = NULL)
pB <- plot_mainfigure(combined_all, c("B*07:02", "B*27:05"), base_size = 11.5, label_size = NULL)
pC <- plot_mainfigure(combined_all, c("C*07:01", "C*02:02"), base_size = 11.5, label_size = NULL)
pDRB <- plot_mainfigure(combined_all, c("DRB1*15:01", "DRB1*08:01"))
pDQB <- plot_mainfigure(combined_all, c("DQB1*06:02", "DQB1*04:02"))

combined_mainfigure <- ggarrange(pA, pB, pC, pDRB, pDQB,
                                labels = c("A","B","C","D","E"),
                                common.legend = TRUE, nrow = 1,
                                align = "h", legend = "bottom")

ggsave(file.path(OUTPUT_DIR,"combined_mainfigure.png"),
       combined_mainfigure, units = "mm", width = 225, height = 55, dpi = 800)
ggsave(file.path(OUTPUT_DIR,"combined_mainfigure.pdf"),
       combined_mainfigure, units = "mm", width = 225, height = 55)
