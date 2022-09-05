here::i_am(file.path("GP-linear", "01a_do-param-grid.R"))

library(here)
library(tidyverse)
library(glue)
library(cowplot)


# 1. Parameters --------------------------------------------------------------------------
project_folder <- "GP-linear"

par_file <- "pars_AroAdm_BayesC.csv"



pars_full <- read_csv(here(
  project_folder, "parameters", par_file
))
pars_full


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

DATA <- DATA %>% 
  relocate(model, part, trait, m_set, pi0, MSE, bin_ce) %>% 
  arrange(trait)

plot_data <- DATA
plot_data <- filter(plot_data, ! is.na(MSE))

plot_data

title <- ggdraw() + draw_label(glue("{unique(plot_data$part)}_{unique(plot_data$model)}"))

p1 <- ggplot(plot_data) +
  geom_col(aes(x = m_set, y =  MSE, fill = factor(pi0)), position = "dodge") +
  facet_grid(cols = vars(trait)) +
  labs(fill = "pi0")

p2 <- ggplot(plot_data) +
  geom_col(aes(x = m_set, y =  cor, fill = factor(pi0)), position = "dodge") +
  ylim(range(plot_data$cor)[1], 1) +
  facet_grid(cols = vars(trait)) +
  labs(fill = "pi0")


# ---

plot_data <- DATA
plot_data <- filter(plot_data, ! is.na(bin_ce))

plot_data
p3 <- ggplot(plot_data) +
  geom_col(aes(x = m_set, y = bin_ce, fill = factor(pi0)), position = "dodge") +
  facet_grid(cols = vars(trait)) +
  labs(fill = "pi0")

grid <- plot_grid(p1, p2, p3, nrow = 3)


png(here(
  project_folder, "plots", glue("res_{unique(plot_data$model)}_{unique(plot_data$part)}.png")
), units = "in", width = 13, height = 8, res = 300)
plot_grid(title, grid, nrow = 2, rel_heights = c(0.1, 1))
dev.off()

