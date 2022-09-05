# Check wd
getwd()

library(tidyverse)
library(glue)
library(cowplot)

# Load data ------------------------------------------------------------------------------

py_data <- read_csv(
  file.path("final_results", "final_results.csv")
)

r_data <- read_csv(
  file.path("final_results", "res_10-fold-cv_AroAdm.csv")
)

# Make data compatible -------------------------------------------------------------------

str(r_data)

r_data <- r_data %>% 
  mutate(experiment = "linear") %>% 
  rename("model_type" = model,
         "binary_crossentropy" = bin_ce) %>% 
  mutate(
    grain.weight_cor = ifelse(trait == "grain.weight", cor, NA),
    time.to.flowering.from.sowing_cor = ifelse(trait == "time.to.flowering.from.sowing", cor, NA),
    culm.diameter.1st.internode_accuracy = ifelse(trait == "culm.diameter.1st.internode", accuracy, NA),
    leaf.senescence_accuracy = ifelse(trait == "leaf.senescence", accuracy, NA)
  ) %>%
  mutate(
    culm.diameter.1st.internode_AUC = ifelse(trait == "culm.diameter.1st.internode", AUC, NA),
    leaf.senescence_AUC = ifelse(trait == "leaf.senescence", AUC, NA)
  )

str(py_data)

py_data <- py_data %>% 
  mutate(experiment = "ANN") %>% 
  rename("trait" = traits,
         "MSE" = mse) %>% 
  mutate(
    input_type = recode(input_type,
                        "single-input" = "single.in",
                        "multi-input" = "multi.in"),
    output_type = recode(output_type,
                        "single-output" = "single.out",
                        "multi-output" = "multi.out"),
    model_type = recode(model_type,
                        "top-markers-10k" = "top-m")
  )

# setdiff(
#   colnames(py_data),
#   colnames(r_data)
# )
# 
# intersect(
#   colnames(py_data),
#   colnames(r_data)
# )


# Join data ------------------------------------------------------------------------------

DATA <- bind_rows(py_data, r_data)

# For ANN, complex model_type
DATA <- DATA %>% 
  mutate(model_type = ifelse(
    experiment == "ANN", paste(model_type, architecture, input_type, sep = "_"), model_type
  ))


colnames(DATA)
# Order columns
DATA <- DATA %>% 
  relocate(trait, experiment, model_type, input_type, m_set, part,
           MSE, binary_crossentropy,
           grain.weight_cor, time.to.flowering.from.sowing_cor,
           culm.diameter.1st.internode_accuracy, leaf.senescence_accuracy,
           culm.diameter.1st.internode_AUC, leaf.senescence_AUC,
           accuracy, cor)

# Code all m_set as SNP / all
DATA <- DATA %>% 
  mutate(m_set = recode(
    m_set,
    all.6ker = "all",
    top.all = "all",
    top.SNP = "SNP"
  )) %>% 
  mutate(m_set = factor(m_set, levels = c("SNP", "all")))

# Recode multiple traits
DATA <- DATA %>% 
  mutate(trait = recode(
    trait,
    "['culm.diameter.1st.internode', 'leaf.senescence', 'grain.weight', 'time.to.flowering.from.sowing']" = "multi.out_4.traits",
    "['culm.diameter.1st.internode', 'leaf.senescence']" = "multi.out_2.bin.traits",
    "['grain.weight', 'time.to.flowering.from.sowing']" = "multi.out_2.cont.traits"
  ))


# Order model_type for plot
# levels(factor(DATA$model_type))
DATA$model_type <- factor(DATA$model_type,
                          levels = c("RKHS", "BayesC",
                                     "kerPC_MLP_single.in", "kerPC_MLP_multi.in",
                                     "top-m_MLP_single.in", "top-m_MLP_multi.in",
                                     "top-m_CNN_single.in"))

# Save single_output traits for plotting
single_traits <- DATA %>% 
  filter(output_type == "single.out") %>% 
  select(trait) %>% 
  unique() %>% 
  unlist() %>% 
  unname()
single_traits

# Plots ----------------------------------------------------------------------------------

# TRAIT <- "culm.diameter.1st.internode"
# single_out_var <- "binary_crossentropy"
# multi_out_var <- "culm.diameter.1st.internode_binary_crossentropy"

make_plot <- function(plot_data, TRAIT, single_out_var, multi_out_var) {
  
  # Get y axis range + minimum median of 10-fold-cv ######################################
  
  plot_data <- DATA
  plot_data <- plot_data %>% 
    filter(trait == TRAIT) %>% 
    mutate(trait = "single.out")
  range1 <- range(plot_data[, single_out_var])
  
  median1 <- plot_data %>% 
    filter(part != "AroAdm") %>% 
    group_by(trait, model_type, m_set) %>% 
    summarise(median = median(!!sym(single_out_var))) %>% 
    ungroup()
  
  stat_1 <- plot_data %>% 
    filter(part == "AroAdm") %>% 
    summarise(stat = min(!!sym(single_out_var)))
  
  # ---
  
  plot_data <- DATA
  plot_data <- plot_data %>% 
    filter(! trait %in% single_traits,
           ! is.na(!!sym(multi_out_var)))
  range2 <- range(plot_data[, multi_out_var])
  
  median2 <- plot_data %>% 
    filter(part != "AroAdm") %>% 
    group_by(trait, model_type, m_set) %>% 
    summarise(median = median(!!sym(multi_out_var))) %>% 
    ungroup()
  
  stat_2 <- plot_data %>% 
    filter(part == "AroAdm") %>% 
    summarise(stat = min(!!sym(multi_out_var)))
  
  # ---
  
  
  MEDIAN_table <- bind_rows(median1, median2)
  MEDIAN <- min(MEDIAN_table$median)
  
  stat_table <- bind_rows(stat_1, stat_2)
  stat <- min(stat_table$stat)
  
  
  ylim <- c(0, max(range1[2], range2[2]))
  # mean_median <- mean(MEDIAN_table$median)
  # ylim <- c(0, 2 * mean_median)
  
  # Plot ###################################################################################
  
  title <- ggdraw() + draw_label(TRAIT)
  
  # Single.out
  
  plot_data <- DATA
  plot_data <- plot_data %>% 
    filter(trait == TRAIT) %>% 
    mutate(trait = "single.out")
  
  p1 <- ggplot(plot_data) +
    geom_boxplot(data = filter(plot_data, part != "AroAdm"),
                 aes(x = model_type, y = !!sym(single_out_var), color = m_set)) +
    geom_point(data = filter(plot_data, part == "AroAdm"),
               aes(x = model_type, y = !!sym(single_out_var), fill = m_set),
               position = position_dodge(width = 0.75), shape = 3, size = 4) +
    facet_grid(cols = vars(trait)) +
    ylim(ylim) +
    theme_bw() +
    # theme(legend.position = "none") +
    xlab("") +
    geom_hline(yintercept = MEDIAN, linetype = 2, color = "darkorange2") +
    geom_hline(yintercept = stat, linetype = 2)
  p1
  
  # Multi.out
  
  plot_data <- DATA
  plot_data <- plot_data %>% 
    filter(! trait %in% single_traits,
           ! is.na(!!sym(multi_out_var)))
  
  
  p2 <- ggplot() +
    geom_boxplot(data = filter(plot_data, part != "AroAdm"),
                 aes(x = model_type, y = !!sym(multi_out_var),
                     color = m_set)) +
    geom_point(data = filter(plot_data, part == "AroAdm"),
               aes(x = model_type, y = !!sym(multi_out_var),
                   fill = m_set),
               position = position_dodge(width = 0.75), shape = 3, size = 4) +
    facet_grid(cols = vars(trait)) +
    ylim(ylim) +
    theme_bw() +
    labs(y = single_out_var) +
    xlab("") +
    geom_hline(yintercept = MEDIAN, linetype = 2, color = "darkorange2") +
    geom_hline(yintercept = stat, linetype = 2)
  p2
  
  # Merge
  
  grid <- plot_grid(p1, p2, nrow = 2)
  
  final <- plot_grid(title, grid, nrow = 2, rel_heights = c(0.05, 1))
  print(final)
  
  png(file.path("plots", glue("loss_plot_{TRAIT}.png")),
      units = "in", width = 10, height = 8, res = 300)
  print(final)
  dev.off()
  
}

make_plot(plot_data = DATA,
          TRAIT = "culm.diameter.1st.internode",
          single_out_var = "binary_crossentropy",
          multi_out_var = "culm.diameter.1st.internode_binary_crossentropy")

make_plot(plot_data = DATA,
          TRAIT = "leaf.senescence",
          single_out_var = "binary_crossentropy",
          multi_out_var = "leaf.senescence_binary_crossentropy")

make_plot(plot_data = DATA,
          TRAIT = "grain.weight",
          single_out_var = "MSE",
          multi_out_var = "grain.weight_mse")

# make_plot(plot_data = DATA,
#           TRAIT = "time.to.flowering.from.sowing",
#           single_out_var = "MSE",
#           multi_out_var = "time.to.flowering.from.sowing_mse")

# ----------------------------------------------------------------------------------------

make_plot_max <- function(plot_data, TRAIT, single_out_var, multi_out_var, y_lab, y_lim_0 = TRUE) {
  
  # Get y axis range + minimum median of 10-fold-cv ######################################
  
  plot_data <- DATA
  plot_data <- plot_data %>% 
    filter(trait == TRAIT) %>% 
    mutate(trait = "single.out")
  range1 <- range(plot_data[, single_out_var])
  
  median1 <- plot_data %>% 
    filter(part != "AroAdm") %>% 
    group_by(trait, model_type, m_set) %>% 
    summarise(median = median(!!sym(single_out_var))) %>% 
    ungroup()
  
  stat_1 <- plot_data %>% 
    filter(part == "AroAdm") %>% 
    summarise(stat = max(!!sym(single_out_var)))
  
  # ---
  
  plot_data <- DATA
  plot_data <- plot_data %>% 
    filter(! trait %in% single_traits,
           ! is.na(!!sym(multi_out_var)))
  range2 <- range(plot_data[, multi_out_var])
  
  median2 <- plot_data %>% 
    filter(part != "AroAdm") %>% 
    group_by(trait, model_type, m_set) %>% 
    summarise(median = median(!!sym(multi_out_var))) %>% 
    ungroup()
  
  stat_2 <- plot_data %>% 
    filter(part == "AroAdm") %>% 
    summarise(stat = max(!!sym(multi_out_var)))
  
  # ---
  
  if (y_lim_0) {
    ylim <- c(0, max(range1[2], range2[2]))
  } else {
    ylim <- c(min(range1[1], range2[1]), max(range1[2], range2[2]))
  }
  
  MEDIAN <- bind_rows(median1, median2)
  
  # tmp <- filter(MEDIAN, model_type == "top-m_MLP_single.in")
  # snp <- tmp %>% 
  #   filter(m_set == "SNP") %>% 
  #   select(median)
  # all <- tmp %>% 
  #   filter(m_set == "all") %>% 
  #   select(median)
  # (all / snp - 1) * 100
  
  MEDIAN <- max(MEDIAN$median)
  
    stat_table <- bind_rows(stat_1, stat_2)
  stat <- max(stat_table$stat)
  
  # Plot ###################################################################################
  
  title <- ggdraw() + draw_label(TRAIT)
  
  # Single.out
  
  plot_data <- DATA
  plot_data <- plot_data %>% 
    filter(trait == TRAIT) %>% 
    mutate(trait = "single.out")
  
  p1 <- ggplot(plot_data) +
    geom_boxplot(data = filter(plot_data, part != "AroAdm"),
                 aes(x = model_type, y = !!sym(single_out_var), color = m_set)) +
    geom_point(data = filter(plot_data, part == "AroAdm"),
               aes(x = model_type, y = !!sym(single_out_var), fill = m_set),
               position = position_dodge(width = 0.75), shape = 3, size = 4) +
    facet_grid(cols = vars(trait)) +
    ylim(ylim) +
    theme_bw() +
    # theme(legend.position = "none") +
    xlab("") +
    labs(y = y_lab) +
    geom_hline(yintercept = MEDIAN, linetype = 2, color = "darkorange2") +
    geom_hline(yintercept = stat, linetype = 2)
  p1
  
  # Multi.out
  
  plot_data <- DATA
  plot_data <- plot_data %>% 
    filter(! trait %in% single_traits,
           ! is.na(!!sym(multi_out_var)))
  
  
  p2 <- ggplot() +
    geom_boxplot(data = filter(plot_data, part != "AroAdm"),
                 aes(x = model_type, y = !!sym(multi_out_var),
                     color = m_set)) +
    geom_point(data = filter(plot_data, part == "AroAdm"),
               aes(x = model_type, y = !!sym(multi_out_var),
                   fill = m_set),
               position = position_dodge(width = 0.75), shape = 3, size = 4) +
    facet_grid(cols = vars(trait)) +
    ylim(ylim) +
    theme_bw() +
    xlab("") +
    labs(y = y_lab) +
    geom_hline(yintercept = MEDIAN, linetype = 2, color = "darkorange2") +
    geom_hline(yintercept = stat, linetype = 2)
  p2
  
  # Merge
  
  grid <- plot_grid(p1, p2, nrow = 2)
  
  final <- plot_grid(title, grid, nrow = 2, rel_heights = c(0.05, 1))
  print(final)
  
  png(file.path("plots", glue("{y_lab}_plot_{TRAIT}.png")),
      units = "in", width = 10, height = 8, res = 300)
  print(final)
  dev.off()
  
}

make_plot_max(plot_data = DATA,
          TRAIT = "time.to.flowering.from.sowing",
          single_out_var = "time.to.flowering.from.sowing_cor",
          multi_out_var = "time.to.flowering.from.sowing_cor",
          y_lab = "cor",
          y_lim_0 = FALSE)

make_plot_max(plot_data = DATA,
              TRAIT = "grain.weight",
              single_out_var = "grain.weight_cor",
              multi_out_var = "grain.weight_cor",
              y_lab = "cor",
              y_lim_0 = FALSE)

# ---

make_plot_max(plot_data = DATA,
              TRAIT = "culm.diameter.1st.internode",
              single_out_var = "accuracy",
              multi_out_var = "culm.diameter.1st.internode_accuracy",
              y_lab = "accuracy",
              y_lim_0 = FALSE)

make_plot_max(plot_data = DATA,
              TRAIT = "leaf.senescence",
              single_out_var = "accuracy",
              multi_out_var = "leaf.senescence_accuracy",
              y_lab = "accuracy",
              y_lim_0 = FALSE)

# ---

make_plot_max(plot_data = DATA,
              TRAIT = "culm.diameter.1st.internode",
              single_out_var = "culm.diameter.1st.internode_AUC",
              multi_out_var = "culm.diameter.1st.internode_AUC",
              y_lab = "AUC",
              y_lim_0 = FALSE)

make_plot_max(plot_data = DATA,
              TRAIT = "leaf.senescence",
              single_out_var = "leaf.senescence_AUC",
              multi_out_var = "leaf.senescence_AUC",
              y_lab = "AUC",
              y_lim_0 = FALSE)


# ----

make_plot_custom <- function(plot_data, TRAIT, single_out_var, multi_out_var) {
  
  # Get y axis range + minimum median of 10-fold-cv ######################################
  
  plot_data <- DATA
  plot_data <- plot_data %>% 
    filter(trait == TRAIT) %>% 
    mutate(trait = "single.out")
  range1 <- range(plot_data[, single_out_var])
  
  median1 <- plot_data %>% 
    filter(part != "AroAdm") %>% 
    group_by(trait, model_type, m_set) %>% 
    summarise(median = median(!!sym(single_out_var))) %>% 
    ungroup()
  
  stat_1 <- plot_data %>% 
    filter(part == "AroAdm") %>% 
    summarise(stat = min(!!sym(single_out_var)))
  
  # ---
  
  plot_data <- DATA
  plot_data <- plot_data %>% 
    filter(! trait %in% single_traits,
           ! is.na(!!sym(multi_out_var)))
  range2 <- range(plot_data[, multi_out_var])
  
  median2 <- plot_data %>% 
    filter(part != "AroAdm") %>% 
    group_by(trait, model_type, m_set) %>% 
    summarise(median = median(!!sym(multi_out_var))) %>% 
    ungroup()
  
  stat_2 <- plot_data %>% 
    filter(part == "AroAdm") %>% 
    summarise(stat = min(!!sym(multi_out_var)))
  
  # ---
  
  
  MEDIAN_table <- bind_rows(median1, median2)
  MEDIAN <- min(MEDIAN_table$median)
  
  stat_table <- bind_rows(stat_1, stat_2)
  stat <- min(stat_table$stat)
  
  ylim <- c(0, 0.2)
  # mean_median <- mean(MEDIAN_table$median)
  # ylim <- c(0, 2 * mean_median)
  
  
  # Plot ###################################################################################
  
  title <- ggdraw() + draw_label(TRAIT)
  
  # Single.out
  
  plot_data <- DATA
  plot_data <- plot_data %>% 
    filter(trait == TRAIT) %>% 
    mutate(trait = "single.out")
  
  p1 <- ggplot(plot_data) +
    geom_boxplot(data = filter(plot_data, part != "AroAdm"),
                 aes(x = model_type, y = !!sym(single_out_var), color = m_set)) +
    geom_point(data = filter(plot_data, part == "AroAdm"),
               aes(x = model_type, y = !!sym(single_out_var), fill = m_set),
               position = position_dodge(width = 0.75), shape = 3, size = 4) +
    facet_grid(cols = vars(trait)) +
    ylim(ylim) +
    theme_bw() +
    # theme(legend.position = "none") +
    xlab("") +
    geom_hline(yintercept = MEDIAN, linetype = 2, color = "darkorange2") +
    geom_hline(yintercept = stat, linetype = 2)
  p1
  
  # Multi.out
  
  plot_data <- DATA
  plot_data <- plot_data %>% 
    filter(! trait %in% single_traits,
           ! is.na(!!sym(multi_out_var)))
  
  
  p2 <- ggplot() +
    geom_boxplot(data = filter(plot_data, part != "AroAdm"),
                 aes(x = model_type, y = !!sym(multi_out_var),
                     color = m_set)) +
    geom_point(data = filter(plot_data, part == "AroAdm"),
               aes(x = model_type, y = !!sym(multi_out_var),
                   fill = m_set),
               position = position_dodge(width = 0.75), shape = 3, size = 4) +
    facet_grid(cols = vars(trait)) +
    ylim(ylim) +
    theme_bw() +
    xlab("") +
    labs(y = single_out_var) +
    geom_hline(yintercept = MEDIAN, linetype = 2, color = "darkorange2") +
    geom_hline(yintercept = stat, linetype = 2)
  p2
  
  # Merge
  
  grid <- plot_grid(p1, p2, nrow = 2)
  
  final <- plot_grid(title, grid, nrow = 2, rel_heights = c(0.05, 1))
  print(final)
  
  png(file.path("plots", glue("loss_plot_{TRAIT}.png")),
      units = "in", width = 10, height = 8, res = 300)
  print(final)
  dev.off()
  
}


make_plot_custom(plot_data = DATA,
          TRAIT = "time.to.flowering.from.sowing",
          single_out_var = "MSE",
          multi_out_var = "time.to.flowering.from.sowing_mse")

