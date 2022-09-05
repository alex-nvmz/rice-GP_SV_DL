here::i_am(file.path("GP-linear", "01a_do-param-grid.R"))

library(here)
library(tidyverse)
library(glue)
library(cowplot)
library(ungeviz)


# 1. Parameters --------------------------------------------------------------------------
project_folder <- "GP-linear"

par_files <- c("pars_10-fold-cv_RKHS.csv",
               "pars_10-fold-cv_BayesC.csv",
               "pars_AroAdm_RKHS.csv",
               "pars_AroAdm_BayesC.csv")


# 2. Get data ----------------------------------------------------------------------------

pars_full <- NULL
for (par_file in par_files) {
  tmp <- read_csv(here(
    project_folder, "parameters", par_file
  ))
  pars_full <- bind_rows(pars_full, tmp)
}

DATA <- NULL

for (i in 1:nrow(pars_full)) {
  pars <- pars_full[i, ]
  
  model <- pars[["model"]]
  part <- pars[["part"]]
  trait <- pars[["trait"]]
  
  if (model == "BayesC") {
    m_set <- pars[["m_set"]]
    n_mark <- pars[["n_mark"]]
    p0 <- pars[["p0"]]
    pi0 <- pars[["pi0"]]
    
    out_folder <- glue("out_{model}_{part}_{trait}_{m_set}-{n_mark}_p0-{p0}_pi0-{pi0}")
    
    tmp <- read_csv(here(
      project_folder, out_folder, "results.csv"
    ))
    
    tmp <- cbind(tmp, model, part, trait, m_set, n_mark, p0, pi0)
    DATA <- bind_rows(DATA, tmp)
    
    
  } else if (model == "RKHS") {
    m_set <- pars[["m_set"]]
    
    out_folder <- glue("out_{model}_{part}_{trait}_{m_set}")
    
    tmp <- read_csv(here(
      project_folder, out_folder, "results.csv"
    ))
    tmp <- cbind(tmp, model, part, trait, m_set)
    
    DATA <- bind_rows(DATA, tmp)
    
  } else {
    print(glue("Unknown BGLR model: {model}"))
    quit(save = "no")
  }
}

DATA <- filter(DATA, pi0 == 0.01 | is.na(pi0))

DATA <- DATA %>% 
  relocate(model, part, trait, m_set, pi0, MSE, bin_ce) %>% 
  arrange(trait)

write_csv(DATA, here(
  project_folder, "summary", "res_10-fold-cv_AroAdm.csv"
))


# DATA summaries

DATA_sum_10.cv <- DATA %>% 
  filter(part != "AroAdm") %>% 
  group_by(model, trait, m_set) %>% 
  summarise(n = n(),
            MSE_mean = mean(MSE),
            cor_mean = mean(cor),
            bin.ce_mean = mean(bin_ce),
            accuracy_mean = mean(accuracy),
            MSE_sd = sd(MSE),
            cor_sd = sd(cor),
            bin.ce_sd = sd(bin_ce),
            accuracy_sd = sd(accuracy)) %>% 
  ungroup()

write_csv(DATA_sum_10.cv, here(
  project_folder, "summary", "summary_10-fold-cv.csv"
))

DATA_sum_all <- DATA %>% 
  group_by(model, trait, m_set) %>% 
  summarise(n = n(),
            MSE_mean = mean(MSE),
            cor_mean = mean(cor),
            bin.ce_mean = mean(bin_ce),
            accuracy_mean = mean(accuracy),
            MSE_sd = sd(MSE),
            cor_sd = sd(cor),
            bin.ce_sd = sd(bin_ce),
            accuracy_sd = sd(accuracy)) %>% 
  ungroup()

write_csv(DATA_sum_all, here(
  project_folder, "summary", "summary_10-fold-cv_AroAdm.csv"
))


DATA_sum_AroAdm <- DATA %>% 
  filter(part == "AroAdm") %>% 
  group_by(model, trait, m_set) %>% 
  summarise(n = n(),
            MSE_mean = mean(MSE),
            cor_mean = mean(cor),
            bin.ce_mean = mean(bin_ce),
            accuracy_mean = mean(accuracy),
            MSE_sd = sd(MSE),
            cor_sd = sd(cor),
            bin.ce_sd = sd(bin_ce),
            accuracy_sd = sd(accuracy)) %>% 
  ungroup()

# 3. Bar plot ----------------------------------------------------------------------------

# plot_data <- DATA_sum_10.cv
# plot_data <- filter(plot_data, ! is.na(MSE_mean) & model == "RKHS")
# 
# 
# title <- ggdraw() + draw_label(glue("10-fold-cv_{unique(DATA$model)}"))
# 
# p1 <- ggplot(plot_data) +
#   geom_col(aes(x = m_set, y = MSE_mean), position = "dodge", fill = "lightskyblue") +
#   geom_errorbar(aes(x = m_set, ymin = MSE_mean - MSE_sd, ymax = MSE_mean + MSE_sd),
#                 width = 0.3, color = "darkorange", size = 0.9) +
#   facet_grid(cols = vars(trait),
#              rows = vars(model))
# p1
# 
# p2 <- ggplot(plot_data) +
#   geom_col(aes(x = m_set, y = cor_mean), position = "dodge", fill = "lightskyblue") +
#   geom_errorbar(aes(x = m_set, ymin = cor_mean - cor_sd, ymax = cor_mean + cor_sd),
#                 width = 0.3, color = "darkorange", size = 0.9) +
#   ylim(0, 1) +
#   facet_grid(cols = vars(trait),
#              rows = vars(model))
# p2
# 
# # ---
# 
# plot_data <- DATA_sum_10.cv
# plot_data <- filter(plot_data, ! is.na(bin.ce_mean) & model == "RKHS")
# 
# p3 <- ggplot(plot_data) +
#   geom_col(aes(x = m_set, y = bin.ce_mean), position = "dodge", fill = "lightskyblue") +
#   geom_errorbar(aes(x = m_set, ymin = bin.ce_mean - bin.ce_sd, ymax = bin.ce_mean + bin.ce_sd),
#                 width = 0.3, color = "darkorange", size = 0.9) +
#   facet_grid(cols = vars(trait),
#              rows = vars(model))
# p3
# 
# grid <- plot_grid(p1, p2, p3, nrow = 3)
# 
# 
# png(here(
#   project_folder, "plots", glue("res_10-fold-cv_{unique(DATA$model)}.png")
# ), units = "in", width = 13, height = 9, res = 300)
# plot_grid(title, grid, nrow = 2, rel_heights = c(0.1, 1))
# dev.off()

# 4. Dot plot ----------------------------------------------------------------------------

# 4.1. RKHS analyses #####################################################################

# 4.1.1. 10-fold-cv ======================================================================
# tag <- "10-fold-cv"
# 
# # Continuous traits
# plot_data <- DATA
# 
# plot_data <- plot_data %>% 
#   filter(model == "RKHS") %>%  
#   filter(part != "AroAdm") %>% 
#   filter(! is.na(MSE))
# 
# plot_data$m_set <- factor(plot_data$m_set, levels = c("SNP", "all.6ker"))
# 
# # binwidth = 1/35 * diff(range(plot_data$MSE)
# # rows = vars(model)
# 
# title <- ggdraw() + draw_label(glue("{tag}_{unique(plot_data$model)}"))
# 
# p1 <- ggplot(plot_data) +
#   geom_dotplot(aes(x = m_set, y = MSE, fill = m_set),
#                binaxis = "y", stackdir = "center",
#                binwidth = 1/35 * diff(range(plot_data$MSE))) +
#   stat_summary(aes(x = m_set, y = MSE), geom = "hpline", fun = "mean", width = 0.4) +
#   ylim(0, range(plot_data$MSE)[2]) +
#   facet_grid(cols = vars(trait)) +
#   theme_bw() +
#   theme(axis.ticks.x = element_blank(),
#         axis.text.x = element_blank(),
#         legend.position="none") +
#   xlab("")
# p1
# 
# p2 <- ggplot(plot_data) +
#   geom_dotplot(aes(x = m_set, y = cor, fill = m_set),
#                binaxis = "y", stackdir = "center",
#                binwidth = 1/30 * diff(range(plot_data$cor))) +
#   stat_summary(aes(x = m_set, y = cor), geom = "hpline", fun = "mean", width = 0.4) +
#   ylim(0, 1) +
#   facet_grid(cols = vars(trait)) +
#   theme_bw() +
#   theme(axis.ticks.x = element_blank(),
#         axis.text.x = element_blank()) +
#   xlab("")
# p2
# 
# # Binary traits
# plot_data <- DATA
# 
# plot_data <- plot_data %>% 
#   filter(model == "RKHS") %>%  
#   filter(part != "AroAdm") %>% 
#   filter(! is.na(bin_ce))
# 
# plot_data$m_set <- factor(plot_data$m_set, levels = c("SNP", "all.6ker"))
# 
# p3 <- ggplot(plot_data) +
#   geom_dotplot(aes(x = m_set, y = bin_ce, fill = m_set),
#                binaxis = "y", stackdir = "center",
#                binwidth = 1/20 * diff(range(plot_data$bin_ce))) +
#   stat_summary(aes(x = m_set, y = bin_ce), geom = "hpline", fun = "mean", width = 0.4) +
#   ylim(0, range(plot_data$bin_ce)[2]) +
#   facet_grid(cols = vars(trait)) +
#   theme_bw() +
#   theme(axis.ticks.x = element_blank(),
#         axis.text.x = element_blank()) +
#   xlab("")
# p3
# 
# grid <- plot_grid(p1, p3, ncol = 2)
# grid
# 
# final <- plot_grid(title, grid, nrow = 2, rel_heights = c(0.05, 1))
# final
# 
# png(here(
#   project_folder, "plots", glue("dotplot_{tag}_{unique(plot_data$model)}.png")
# ), units = "in", width = 13, height = 9, res = 300)
# print(final)
# dev.off()
# 
# png(here(
#   project_folder, "plots", glue("dotplot_{tag}_{unique(plot_data$model)}_cor.png")
# ), units = "in", width = 13, height = 9, res = 300)
# plot_grid(title, p2, nrow = 2, rel_heights = c(0.05, 1))
# dev.off()

# 4.1.2. All (11) folds ==================================================================
tag <- "10-fold-cv_AroAdm"

# Continuous traits
plot_data <- DATA

plot_data <- plot_data %>% 
  filter(model == "RKHS") %>%  
  filter(! is.na(MSE))

plot_data$m_set <- factor(plot_data$m_set, levels = c("SNP", "all.6ker"))

title <- ggdraw() + draw_label(glue("{tag}_{unique(plot_data$model)}"))

p1 <- ggplot(plot_data) +
  geom_dotplot(aes(x = m_set, y = MSE, fill = m_set),
               binaxis = "y", stackdir = "center",
               binwidth = 1/35 * diff(range(plot_data$MSE))) +
  stat_summary(aes(x = m_set, y = MSE), geom = "hpline", fun = "mean", width = 0.4) +
  geom_text(aes(x = m_set, y = MSE, label = ifelse(part == "AroAdm", part, "")),
            hjust = -0.2, vjust = 0) +
  ylim(0, range(plot_data$MSE)[2]) +
  facet_grid(cols = vars(trait)) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position="none") +
  xlab("")
p1

p2 <- ggplot(plot_data) +
  geom_dotplot(aes(x = m_set, y = cor, fill = m_set),
               binaxis = "y", stackdir = "center",
               binwidth = 1/30 * diff(range(plot_data$cor))) +
  stat_summary(aes(x = m_set, y = cor), geom = "hpline", fun = "mean", width = 0.4) +
  geom_text(aes(x = m_set, y = cor, label = ifelse(part == "AroAdm", part, "")),
            hjust = -0.2, vjust = 0) +
  ylim(0, 1) +
  facet_grid(cols = vars(trait)) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  xlab("")
p2

# Binary traits
plot_data <- DATA

plot_data <- plot_data %>% 
  filter(model == "RKHS") %>%  
  filter(! is.na(bin_ce))

plot_data$m_set <- factor(plot_data$m_set, levels = c("SNP", "all.6ker"))

p3 <- ggplot(plot_data) +
  geom_dotplot(aes(x = m_set, y = bin_ce, fill = m_set),
               binaxis = "y", stackdir = "center",
               binwidth = 1/20 * diff(range(plot_data$bin_ce))) +
  stat_summary(aes(x = m_set, y = bin_ce), geom = "hpline", fun = "mean", width = 0.4) +
  geom_text(aes(x = m_set, y = bin_ce, label = ifelse(part == "AroAdm", part, "")),
            hjust = -0.2, vjust = 0) +
  ylim(0, range(plot_data$bin_ce)[2]) +
  facet_grid(cols = vars(trait)) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  xlab("")
p3

grid <- plot_grid(p1, p3, ncol = 2, rel_widths = c(0.85, 1))
grid

final <- plot_grid(title, grid, nrow = 2, rel_heights = c(0.05, 1))
final

png(here(
  project_folder, "plots", glue("dotplot_{tag}_{unique(plot_data$model)}.png")
), units = "in", width = 13, height = 9, res = 300)
print(final)
dev.off()

png(here(
  project_folder, "plots", glue("dotplot_{tag}_{unique(plot_data$model)}_cor.png")
), units = "in", width = 13, height = 9, res = 300)
plot_grid(title, p2, nrow = 2, rel_heights = c(0.05, 1))
dev.off()

# 4.1.3. AroAdm ==========================================================================
# tag <- "AroAdm"
# 
# # Continuous traits
# plot_data <- DATA
# 
# plot_data <- plot_data %>% 
#   filter(model == "RKHS") %>%  
#   filter(part == "AroAdm") %>% 
#   filter(! is.na(MSE))
# 
# plot_data$m_set <- factor(plot_data$m_set, levels = c("SNP", "all.6ker"))
# 
# title <- ggdraw() + draw_label(glue("{tag}_{unique(plot_data$model)}"))
# 
# p1 <- ggplot(plot_data) +
#   geom_dotplot(aes(x = m_set, y = MSE, fill = m_set),
#                binaxis = "y", stackdir = "center",
#                binwidth = 1/35 * diff(range(plot_data$MSE))) +
#   # stat_summary(aes(x = m_set, y = MSE), geom = "hpline", fun = "mean", width = 0.4) +
#   ylim(0, range(plot_data$MSE)[2]) +
#   facet_grid(cols = vars(trait)) +
#   theme_bw() +
#   theme(axis.ticks.x = element_blank(),
#         axis.text.x = element_blank(),
#         legend.position="none") +
#   xlab("")
# p1
# 
# p2 <- ggplot(plot_data) +
#   geom_dotplot(aes(x = m_set, y = cor, fill = m_set),
#                binaxis = "y", stackdir = "center",
#                binwidth = 1/30 * diff(range(plot_data$cor))) +
#   # stat_summary(aes(x = m_set, y = cor), geom = "hpline", fun = "mean", width = 0.4) +
#   ylim(0, 1) +
#   facet_grid(cols = vars(trait)) +
#   theme_bw() +
#   theme(axis.ticks.x = element_blank(),
#         axis.text.x = element_blank()) +
#   xlab("")
# p2
# 
# # Binary traits
# plot_data <- DATA
# 
# plot_data <- plot_data %>% 
#   filter(model == "RKHS") %>%  
#   filter(part == "AroAdm") %>% 
#   filter(! is.na(bin_ce))
# 
# plot_data$m_set <- factor(plot_data$m_set, levels = c("SNP", "all.6ker"))
# 
# p3 <- ggplot(plot_data) +
#   geom_dotplot(aes(x = m_set, y = bin_ce, fill = m_set),
#                binaxis = "y", stackdir = "center",
#                binwidth = 1/3 * diff(range(plot_data$bin_ce))) +
#   # stat_summary(aes(x = m_set, y = bin_ce), geom = "hpline", fun = "mean", width = 0.4) +
#   ylim(0, range(plot_data$bin_ce)[2]) +
#   facet_grid(cols = vars(trait)) +
#   theme_bw() +
#   theme(axis.ticks.x = element_blank(),
#         axis.text.x = element_blank()) +
#   xlab("")
# p3
# 
# grid <- plot_grid(p1, p3, ncol = 2)
# grid
# 
# final <- plot_grid(title, grid, nrow = 2, rel_heights = c(0.05, 1))
# final
# 
# png(here(
#   project_folder, "plots", glue("dotplot_{tag}_{unique(plot_data$model)}.png")
# ), units = "in", width = 13, height = 9, res = 300)
# print(final)
# dev.off()
# 
# png(here(
#   project_folder, "plots", glue("dotplot_{tag}_{unique(plot_data$model)}_cor.png")
# ), units = "in", width = 13, height = 9, res = 300)
# plot_grid(title, p2, nrow = 2, rel_heights = c(0.05, 1))
# dev.off()


# 4.2. BayesC analyses #####################################################################

tag <- "10-fold-cv_AroAdm"

# Continuous traits
plot_data <- DATA

plot_data <- plot_data %>% 
  filter(model == "BayesC") %>%  
  filter(! is.na(MSE))

plot_data$m_set <- factor(plot_data$m_set, levels = c("top.SNP", "top.all"))

title <- ggdraw() + draw_label(glue("{tag}_{unique(plot_data$model)}"))

p1 <- ggplot(plot_data) +
  geom_dotplot(aes(x = m_set, y = MSE, fill = m_set),
               binaxis = "y", stackdir = "center",
               binwidth = 1/35 * diff(range(plot_data$MSE))) +
  stat_summary(aes(x = m_set, y = MSE), geom = "hpline", fun = "mean", width = 0.4) +
  geom_text(aes(x = m_set, y = MSE, label = ifelse(part == "AroAdm", part, "")),
            hjust = -0.2, vjust = 0) +
  ylim(0, range(plot_data$MSE)[2]) +
  facet_grid(cols = vars(trait)) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        legend.position="none") +
  xlab("")
p1

p2 <- ggplot(plot_data) +
  geom_dotplot(aes(x = m_set, y = cor, fill = m_set),
               binaxis = "y", stackdir = "center",
               binwidth = 1/30 * diff(range(plot_data$cor))) +
  stat_summary(aes(x = m_set, y = cor), geom = "hpline", fun = "mean", width = 0.4) +
  geom_text(aes(x = m_set, y = cor, label = ifelse(part == "AroAdm", part, "")),
            hjust = -0.2, vjust = 0) +
  ylim(range(plot_data$cor)[1], 1) +
  facet_grid(cols = vars(trait)) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  xlab("")
p2

# Binary traits
plot_data <- DATA

plot_data <- plot_data %>% 
  filter(model == "BayesC") %>%  
  filter(! is.na(bin_ce))

plot_data$m_set <- factor(plot_data$m_set, levels = c("top.SNP", "top.all"))

p3 <- ggplot(plot_data) +
  geom_dotplot(aes(x = m_set, y = bin_ce, fill = m_set),
               binaxis = "y", stackdir = "center",
               binwidth = 1/20 * diff(range(plot_data$bin_ce))) +
  stat_summary(aes(x = m_set, y = bin_ce), geom = "hpline", fun = "mean", width = 0.4) +
  geom_text(aes(x = m_set, y = bin_ce, label = ifelse(part == "AroAdm", part, "")),
            hjust = -0.2, vjust = 0) +
  ylim(0, range(plot_data$bin_ce)[2]) +
  facet_grid(cols = vars(trait)) +
  theme_bw() +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  xlab("")
p3

grid <- plot_grid(p1, p3, ncol = 2, rel_widths = c(0.85, 1))
grid

final <- plot_grid(title, grid, nrow = 2, rel_heights = c(0.05, 1))
final

png(here(
  project_folder, "plots", glue("dotplot_{tag}_{unique(plot_data$model)}.png")
), units = "in", width = 13, height = 9, res = 300)
print(final)
dev.off()

png(here(
  project_folder, "plots", glue("dotplot_{tag}_{unique(plot_data$model)}_cor.png")
), units = "in", width = 13, height = 9, res = 300)
plot_grid(title, p2, nrow = 2, rel_heights = c(0.05, 1))
dev.off()


# Stats

DATA %>% 
  filter(part != "AroAdm") %>% 
  summarise(test = mean(size_test),
            train = mean(size_train))
DATA %>% 
  filter(part == "AroAdm") %>% 
  summarise(test = mean(size_test),
            train = mean(size_train))


mean(DATA$size_train)
mean(DATA$size_test)
