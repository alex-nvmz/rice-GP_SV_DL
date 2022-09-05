here::i_am(file.path("GP-linear", "01a_do-param-grid.R"))

library(here)
library(tidyverse)
library(glue)
library(BGLR)
library(pROC)

# 1. Parameters --------------------------------------------------------------------------
project_folder <- "GP-linear"

par_files <- c("pars_10-fold-cv_RKHS.csv",
               "pars_10-fold-cv_BayesC.csv",
               "pars_AroAdm_RKHS.csv",
               "pars_AroAdm_BayesC.csv")


# 2. Get all analyses list  --------------------------------------------------------------

pars_full <- NULL
for (par_file in par_files) {
  tmp <- read_csv(here(
    project_folder, "parameters", par_file
  ))
  pars_full <- bind_rows(pars_full, tmp)
}
pars_full

# 3. For every analysis, compute new metric if appropriate --------------------------------

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
    
  } else if (model == "RKHS") {
    m_set <- pars[["m_set"]]
    
    out_folder <- glue("out_{model}_{part}_{trait}_{m_set}")
    
  } else {
    print(glue("Unknown BGLR model: {model}"))
    quit(save = "no")
  }
  
  # Load data post BGLR analysis #########################################################
  load(here(
    project_folder, out_folder, "env.RData"
  ))
  
  # Variable type
  var_type
  
  if (var_type == "binary") {
    
    # Compute accuracy
    do_accuracy <- function(y_true, y_prob) {
      y_pred <- ifelse(y_prob > 0.5, 1, 0)
      true_count <- sum(y_true == y_pred)  
      accuracy <- true_count / length(y_true)
      return(accuracy)
    }
    
    accuracy <- do_accuracy(y_true = y_true, y_prob = fm$probs[test, 2])
    
    # Compute AUC
    y_prob <- fm$probs[test, 2]
    y_pred <- ifelse(y_prob > 0.5, 1, 0)
    
    auc_res <- auc(y_true, y_pred)
    AUC <- as.numeric(auc_res)
    
    # Save results
    RES[, "accuracy"] <- accuracy
    RES[, "AUC"] <- AUC
    
    # Write
    write_csv(RES, here(
      project_folder, out_folder, "results.csv"
    ))
  }
}

