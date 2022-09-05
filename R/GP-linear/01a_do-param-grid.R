here::i_am(file.path("GP-linear", "01a_do-param-grid.R"))

library(here)
library(tidyverse)
library(glue)

# ----------------------------------------------------------------------------------------
project_folder <- "GP-linear"


pheno <- read_csv(here("data", "main", "pheno_original.csv"))
parts <- read_csv(here("data", "main", "partitions.csv"))

# Partition names
parts <- colnames(parts)[-1]
parts
# Choose traits
traits <- colnames(pheno)[c(4, 10, 11, 13)]
traits

# 1. Hyperparameter tuning (test with AroAdm partition) ----------------------------------

# It also serves as prediction for AroAdm scenario

# RKHS

par_grid <- expand.grid(
  model = c("RKHS"),
  part = parts[11],
  trait = traits,
  m_set = c("SNP", "all.6ker")
)
par_grid

write_csv(par_grid, here(
  project_folder, "parameters", glue("pars_AroAdm_RKHS.csv")
))


# BayesC

par_grid <- expand.grid(
  model = c("BayesC"),
  part = parts[11],
  trait = traits,
  m_set = c("top.SNP", "top.all"),
  n_mark = c(10000),
  p0 = c(5),
  pi0 = c(0.01, 0.1)
)
par_grid

write_csv(par_grid, here(
  project_folder, "parameters", glue("pars_AroAdm_BayesC.csv")
))


# 2. Prediction --------------------------------------------------------------------------

# RKHS

par_grid <- expand.grid(
  model = c("RKHS"),
  part = parts[1:10],
  trait = traits,
  m_set = c("SNP", "all.6ker")
)
par_grid

write_csv(par_grid, here(
  project_folder, "parameters", glue("pars_10-fold-cv_RKHS.csv")
))

# BayesC

par_grid <- expand.grid(
  model = c("BayesC"),
  part = parts[1:10],
  trait = traits,
  m_set = c("top.SNP", "top.all"),
  n_mark = c(10000),
  p0 = c(5),
  pi0 = c(0.01)
)
par_grid

write_csv(par_grid, here(
  project_folder, "parameters", glue("pars_10-fold-cv_BayesC.csv")
))
