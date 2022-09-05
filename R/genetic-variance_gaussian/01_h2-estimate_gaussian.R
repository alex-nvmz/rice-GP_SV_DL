here::i_am(file.path("genetic-variance_gaussian", "01_h2-estimate_gaussian.R"))

library(here)
library(tidyverse)
library(glue)
library(BGLR)

# 1. Parameters --------------------------------------------------------------------------
project_folder <- "genetic-variance_gaussian"

# HARD-CODED ###

n_iter <- 100000
burn_in <- 500


# 2. Load data ---------------------------------------------------------------------------

# Traits
pheno <- read_csv(here(
  "data", "main", "pheno_original.csv"
))
# Choose traits
traits <- colnames(pheno)[- c(1:3)]
traits

# Variable type (used later)
var_table <- tribble(
  ~trait, ~var_type,
  "culm.diameter.1st.internode", "binary",
  "culm.strength", "binary",
  "flag.leaf.angle", "categorical",
  "grain.length", "continuous",
  "grain.width", "continuous",
  "leaf.length", "categorical",
  "leaf.senescence", "binary",
  "grain.weight", "continuous",
  "salt.injury.at.EC12", "binary",
  "time.to.flowering.from.sowing", "continuous",
  "panicle.threshability", "binary"
)

# Features: kernel combinations
K_df <- data.frame(K_1 = "SNP",
                   K_2 = c("RLX-RIX", "MITE-DTX", "DEL", "DUP", "INV"))
K_df

# 3. Perform analysis --------------------------------------------------------------------

# For every kernel combination
for (i in 1:nrow(K_df)) {
  kernel_names <- K_df[i, ]
  print(kernel_names)
  
  # Load kernels
  K_list <- list()
  for (kernel_name in kernel_names) {
    
    K <- read.csv(here(
      "data", "kernels", glue("kernel_{kernel_name}.csv")
    ), row.names = 1) %>% 
      as.matrix()
    
    K_list[[kernel_name]] <- K
  }
  
  # Construct linear predictor
  ETA <- list()
  for (kernel_name in kernel_names) {
    ETA[[kernel_name]] <- list(K = K_list[[kernel_name]],
                               model = "RKHS")
  }
  
  
  # Output folder
  out_folder <- glue("out_RKHS_{kernel_names[1]}_{kernel_names[2]}")
  out_path <- here(project_folder, out_folder)
  if (! dir.exists(out_path)) {dir.create(out_path)}
  
  # ======================================================================================
  
  # Results for this kernel combination
  RES <- NULL
  
  # For every trait
  for (trait in traits) {
    print(trait)
    
    # Get data, variable type and response type
    y <- pheno[, trait, drop = TRUE]
    
    var_type <- var_table[which(var_table$trait == trait), "var_type", drop = TRUE]
    
    response_type <- "gaussian"
    
    # if (var_type == "continuous") {
    #   response_type <- "gaussian"
    #   
    # } else {
    #   response_type <- "ordinal"
    # }
    
    # No scaling
    # if (var_type == "continuous") {
    #   y <- scale(y)
    # }
    
    # ====================================================================================
    
    # Fit the model
    
    cat("\n\n")
    t_start <- Sys.time()
    
    fm <- BGLR(
      y = y,
      response_type = response_type,
      ETA = ETA,
      nIter = n_iter,
      burnIn = burn_in,
      saveAt = here(
        project_folder, out_folder, glue("{trait}_")
      ),
      verbose = FALSE
    )
    
    t_end <- Sys.time()
    print(t_end - t_start)
    cat("\n\n")
    
    # ====================================================================================
    
    # Extract variances and compute heritabilities
    
    # Residual variance (environmental + genetic not captured by markers)
    varE <- fm$varE
    
    # Genetic variance for each marker set
    varU_list <- list()
    for (kernel_name in kernel_names) {
      varU_list[[kernel_name]] <- fm$ETA[[kernel_name]]$varU
    }
    # All genetic variance
    varG <- sum(as_vector(varU_list))
    
    
    # Compute heritabilities
    RES_tmp <- NULL
    RES_tmp[["trait"]] <- trait
    
    for (kernel_name in kernel_names) {
      h2 <- varU_list[[kernel_name]] / (varE + varG)
      
      RES_tmp[[kernel_name]] <- h2
    }
    RES_tmp <- data.frame(RES_tmp)
    
    # Save results
    RES <- rbind(RES, RES_tmp)
    

    # Trace plots ####
    # 
    # varE_t <- scan(here(
    #   project_folder, out_folder, glue("{trait}_varE.dat")
    # ))
    # 
    # varU_t_list <- list()
    # varG_t <- numeric(length = length(varE_t))
    # for (kernel_name in kernel_names) {
    #   varU_t <- scan(here(
    #     project_folder, out_folder, glue("{trait}_ETA_{kernel_name}_varU.dat")
    #   ))
    #   
    #   varU_t_list[[kernel_name]] <- varU_t
    #   varG_t <- varG_t + varU_t
    # }
    # 
    # 
    # for (kernel_name in kernel_names) {
    #   plot(varU_t_list[[kernel_name]],
    #        type = "o", cex = 0.5, col = 4,
    #        ylab = glue("varU {kernel_name}"))
    #   abline(h = mean(varU_t_list[[kernel_name]]), lty = 2)
    # }
    # 
    # plot(varG_t,
    #      type = "o", cex = 0.5, col = 4)
    # abline(h = mean(varG_t), lty = 2)
    # 
    # 
    # plot(varE_t,
    #      type = "o", cex = 0.5, col = 4)
    # abline(h = mean(varE_t), lty = 2)
  
  }
  
  # Save heritabilities for this kernel combination
  
  write_csv(RES, here(
    project_folder, out_folder, "results.csv"
  ))
  
  # Save environment to come back later if necessary
  save.image(here(
    project_folder, out_folder, "env.RData"
  ))
  
}