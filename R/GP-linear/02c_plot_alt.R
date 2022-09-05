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

# ---
pheno <- read_csv(here(
  "data", "main", "pheno_original.csv"
))

traits <- pheno[, c(11, 13)]
traits

var(traits$grain.weight, na.rm = TRUE)
var(traits$time.to.flowering.from.sowing, na.rm = TRUE)






