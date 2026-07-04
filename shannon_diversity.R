library("tidyverse")
library("stringr")
library("vegan")
library("ggpubr")

# constants

POP_ORDER <- c("CEU", "GBR", "FIN", "TSI", "IBS", "MR",
                "late_neolithic_farmers", "early_neolithic_farmers", "sakhtysh")
POP_LABELS <- c("CEU", "GBR", "FIN", "TSI", "IBS", "MR", "LF", "EF", "Sakhtysh")

# functions

bootstrap_locus <- function(data, col1, col2, min_n, n_iter = 100) {
  pops <- unique(data$group)
  rows <- vector("list", length(pops) * n_iter)
  k <- 1L
  for (pop in pops) {
    alleles <- data %>%
      filter(group == pop) %>%
      select(all_of(c(col1, col2)))
    freq <- table(unlist(alleles))
    freq <- freq / sum(freq)
    for (i in seq_len(n_iter)) {
      si <- diversity(table(sample(names(freq), min_n,
                                                prob = freq, replace = TRUE)))
      rows[[k]] <- data.frame(population = pop, iteration = i, SI = si)
      k <- k + 1L
    }
  }
  do.call(rbind, rows)
}

bootstrap_modern_russia <- function(afnd_data, locus_pattern, min_n, n_iter = 100) {
  locus_data <- afnd_data %>%
    filter(str_detect(allele, locus_pattern))
  rows <- vector("list", n_iter)
  for (i in seq_len(n_iter)) {
    si <- diversity(table(sample(locus_data$allele, min_n,
                                               prob = locus_data$allele_freq,
                                               replace = TRUE)))
    rows[[i]] <- data.frame(population = "MR", iteration = i, SI = si)
  }
  do.call(rbind, rows)
}

plot_shannon <- function(df, locus_name, min_n, resolution) {
  ggplot(df, aes(x = population, y = SI)) +
    geom_boxplot() +
    xlab("population") +
    ylab("Shannon index") +
    ylim(0, max(df$SI) * 1.1) +
    labs(title = paste0(locus_name, " diversity"),
         subtitle = paste0("sample size = ", min_n, "; resolution: ", resolution)) +
    scale_x_discrete(limits = POP_ORDER, labels = POP_LABELS) +
    theme_classic()
}

run_kruskal_dunn <- function(df) {
  print(kruskal.test(SI ~ population, data = df))
  dunn.test::dunn.test(df$SI, df$population,
                       method = "bonferroni", alpha = 0.05 / 5, wrap = TRUE)
}

# settings

OUTPUT_DIR <- "results"
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

STUDY_POPS <- c("CEU", "IBS", "FIN", "TSI", "GBR",
                "early_neolithic_farmers", "late_neolithic_farmers", "sakhtysh")
MR_LABEL <- "Russia Nizhny Novgorod, Russians"
N_ITER <- 100

# load data

neol_1 <- read.csv("data/combined_HLA_genotypes_russia_class1.csv",
                             sep = ",", na.strings = c("NA", ""))
neol_2 <- read.csv("data/combined_HLA_genotypes_russia_class2.csv",
                             sep = ",", na.strings = c("NA", ""))
modern_russia_1f <- read.csv("data/russia_AFND_1field.csv") %>%
  filter(population == MR_LABEL)
modern_russia_2f <- read.csv("data/russia_AFND.csv") %>%
  filter(population == MR_LABEL)

# section 1: first-field

class1_1f <- neol_1 %>%
  select(Sample_ID, group,
         HLA_A_1_twodigit, HLA_A_2_twodigit,
         HLA_B_1_twodigit, HLA_B_2_twodigit,
         HLA_C_1_twodigit, HLA_C_2_twodigit) %>%
  drop_na() %>%
  filter(group %in% STUDY_POPS)

class2_1f_dqb <- neol_2 %>%
  select(group, HLA_DQB1_1_twodigit, HLA_DQB1_2_twodigit) %>%
  drop_na() %>%
  filter(group %in% STUDY_POPS)

class2_1f_drb <- neol_2 %>%
  select(Sample_ID, group, HLA_DRB1_1_twodigit, HLA_DRB1_2_twodigit) %>%
  drop_na() %>%
  filter(group %in% STUDY_POPS)

min_n_1f_class1 <- min(table(class1_1f$group))
min_n_1f_dqb <- min(table(class2_1f_dqb$group))
min_n_1f_drb <- min(table(class2_1f_drb$group))

loci_class1_1f <- list(
  list(name = "HLA-A", col1 = "HLA_A_1_twodigit", col2 = "HLA_A_2_twodigit", pattern = "A\\*"),
  list(name = "HLA-B", col1 = "HLA_B_1_twodigit", col2 = "HLA_B_2_twodigit", pattern = "B\\*"),
  list(name = "HLA-C", col1 = "HLA_C_1_twodigit", col2 = "HLA_C_2_twodigit", pattern = "C\\*")
)

bs_1f <- lapply(loci_class1_1f, function(locus) {
  rbind(
    bootstrap_locus(class1_1f, locus$col1, locus$col2,
                    min_n = min_n_1f_class1, n_iter = N_ITER),
    bootstrap_modern_russia(modern_russia_1f, locus$pattern,
                            min_n = min_n_1f_class1, n_iter = N_ITER)
  )
})
names(bs_1f) <- c("HLA_A", "HLA_B", "HLA_C")

bs_1f_dqb <- rbind(
  bootstrap_locus(class2_1f_dqb, "HLA_DQB1_1_twodigit", "HLA_DQB1_2_twodigit",
                  min_n = min_n_1f_dqb, n_iter = N_ITER),
  bootstrap_modern_russia(modern_russia_1f, "DQB1\\*",
                          min_n = min_n_1f_dqb, n_iter = N_ITER)
)

bs_1f_drb <- rbind(
  bootstrap_locus(class2_1f_drb, "HLA_DRB1_1_twodigit", "HLA_DRB1_2_twodigit",
                  min_n = min_n_1f_drb, n_iter = N_ITER),
  bootstrap_modern_russia(modern_russia_1f, "DRB1\\*",
                          min_n = min_n_1f_drb, n_iter = N_ITER)
)

p_1f_A <- plot_shannon(bs_1f$HLA_A, "HLA-A", min_n_1f_class1, "First-field")
p_1f_B <- plot_shannon(bs_1f$HLA_B, "HLA-B", min_n_1f_class1, "First-field")
p_1f_C <- plot_shannon(bs_1f$HLA_C, "HLA-C", min_n_1f_class1, "First-field")
p_1f_DQB <- plot_shannon(bs_1f_dqb, "HLA-DQB1", min_n_1f_dqb, "First-field")
p_1f_DRB <- plot_shannon(bs_1f_drb, "HLA-DRB1", min_n_1f_drb, "First-field")

combined_1f <- ggarrange(p_1f_A, p_1f_B, p_1f_C, p_1f_DRB, p_1f_DQB,
                         common.legend = TRUE, nrow = 3, ncol = 2,
                         align = "v", legend = "bottom")
ggsave(file.path(OUTPUT_DIR,"combined_SI_1field.png"), combined_1f,
       units = "mm", width = 200, height = 200, dpi = 600)

for (locus in names(bs_1f)) {
  cat(locus, "\n")
  run_kruskal_dunn(bs_1f[[locus]])
}
cat("HLA_DQB\n")
run_kruskal_dunn(bs_1f_dqb)
cat("HLA_DRB\n")
run_kruskal_dunn(bs_1f_drb)

# section 2: second-field

class1_2f <- neol_1 %>%
  select(Sample_ID, group, HLA_A_1, HLA_A_2, HLA_B_1, HLA_B_2, HLA_C_1, HLA_C_2) %>%
  drop_na() %>%
  filter(group %in% STUDY_POPS)

class2_2f_dqb <- neol_2 %>%
  select(group, HLA_DQB1_1, HLA_DQB1_2) %>%
  drop_na() %>%
  filter(group %in% STUDY_POPS)

class2_2f_drb <- neol_2 %>%
  select(Sample_ID, group, HLA_DRB1_1, HLA_DRB1_2) %>%
  drop_na() %>%
  filter(group %in% STUDY_POPS)

min_n_2f_class1 <- min(table(class1_2f$group))
min_n_2f_dqb <- min(table(class2_2f_dqb$group))
min_n_2f_drb <- min(table(class2_2f_drb$group))

loci_class1_2f <- list(
  list(name = "HLA-A", col1 = "HLA_A_1", col2 = "HLA_A_2", pattern = "A\\*"),
  list(name = "HLA-B", col1 = "HLA_B_1", col2 = "HLA_B_2", pattern = "B\\*"),
  list(name = "HLA-C", col1 = "HLA_C_1", col2 = "HLA_C_2", pattern = "C\\*")
)

bs_2f <- lapply(loci_class1_2f, function(locus) {
  rbind(
    bootstrap_locus(class1_2f, locus$col1, locus$col2,
                    min_n = min_n_2f_class1, n_iter = N_ITER),
    bootstrap_modern_russia(modern_russia_2f, locus$pattern,
                            min_n = min_n_2f_class1, n_iter = N_ITER)
  )
})
names(bs_2f) <- c("HLA_A", "HLA_B", "HLA_C")

bs_2f_dqb <- rbind(
  bootstrap_locus(class2_2f_dqb, "HLA_DQB1_1", "HLA_DQB1_2",
                  min_n = min_n_2f_dqb, n_iter = N_ITER),
  bootstrap_modern_russia(modern_russia_2f, "DQB1\\*",
                          min_n = min_n_2f_dqb, n_iter = N_ITER)
)

bs_2f_drb <- rbind(
  bootstrap_locus(class2_2f_drb, "HLA_DRB1_1", "HLA_DRB1_2",
                  min_n = min_n_2f_drb, n_iter = N_ITER),
  bootstrap_modern_russia(modern_russia_2f, "DRB1\\*",
                          min_n = min_n_2f_drb, n_iter = N_ITER)
)

p_2f_A <- plot_shannon(bs_2f$HLA_A, "HLA-A", min_n_2f_class1, "Second-field")
p_2f_B <- plot_shannon(bs_2f$HLA_B, "HLA-B", min_n_2f_class1, "Second-field")
p_2f_C <- plot_shannon(bs_2f$HLA_C, "HLA-C", min_n_2f_class1, "Second-field")
p_2f_DQB <- plot_shannon(bs_2f_dqb, "HLA-DQB1", min_n_2f_dqb, "Second-field")
p_2f_DRB <- plot_shannon(bs_2f_drb, "HLA-DRB1", min_n_2f_drb, "Second-field")

combined_2f <- ggarrange(p_2f_A, p_2f_B, p_2f_C, p_2f_DRB, p_2f_DQB,
                         common.legend = TRUE, nrow = 3, ncol = 2,
                         align = "v", legend = "bottom")
ggsave(file.path(OUTPUT_DIR,"combined_SI_2field.png"), combined_2f,
       units = "mm", width = 200, height = 200, dpi = 600)

for (locus in names(bs_2f)) {
  cat(locus, "\n")
  run_kruskal_dunn(bs_2f[[locus]])
}
cat("HLA_DQB\n")
run_kruskal_dunn(bs_2f_dqb)
cat("HLA_DRB\n")
run_kruskal_dunn(bs_2f_drb)

