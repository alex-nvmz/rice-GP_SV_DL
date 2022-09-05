here::i_am(file.path("GP-linear", "01a_do-param-grid.R"))

library(here)
library(tidyverse)
library(glue)
library(BGLR)

# 1. Parameters --------------------------------------------------------------------------
project_folder <- "GP-linear"

# HARD-CODED ###

n_iter <- 100000
burn_in <- 500

# VARIABLE ###

# Parse parameters file and row from command line ---
args <- commandArgs(trailingOnly = T)

par_file <- args[1]
par_row <- as.numeric(args[2])

# ---
# par_file <- "pars_AroAdm_RKHS.csv"
# par_file <- "pars_AroAdm_BayesC.csv"
# par_row <- 5
# par_row <- 7
# ---

# Get pars
pars <- read_csv(here(
  project_folder, "parameters", par_file
))
pars <- pars[par_row, ]
pars

# Extract pars
model <- pars[["model"]]
part <- pars[["part"]]
trait <- pars[["trait"]]

# Model-specific pars
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

# Out folder check
out_path <- here(project_folder, out_folder)
if (! dir.exists(out_path)) {dir.create(out_path)}

# Variable type
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

var_type <- var_table[which(var_table$trait == trait), "var_type", drop = TRUE]

# Response type for BGLR
if (var_type == "continuous") {
  response_type <- "gaussian"
  
} else if (var_type == "binary") {
  response_type <- "ordinal"
  
} else {
  print(glue("Script not implemented for var_type: {var_type}"))
  quit(save = "no")
}

# 2. Load data ---------------------------------------------------------------------------

# Traits
pheno <- read_csv(here(
  "data", "main", "pheno_original.csv"
))

# Partitions
parts <- read_csv(here(
  "data", "main", "partitions.csv"
))
part_guide <- parts[, part]

# Features
# For BayesC, genotype matrices of top markers
if (model == "BayesC") {
  tmp <- if (m_set == "top.SNP") {"SNP"} else if (m_set == "top.all") {"all"} 
  
  X <- read_csv(here(
    "data", "top-markers", glue("geno_top-{n_mark}-{tmp}_{part}_{trait}.csv")
  ))
  X <- as.matrix(X[, - 1])

# For RKHS, marker kernels 
} else if (model == "RKHS") {
  K_list <- list()
  
  # 1 SNP kernel
  if (m_set == "SNP") {
    kernel_names <- c("SNP")
  # 6 kernels
  } else if (m_set == "all.6ker") {
    kernel_names <- c("SNP", "RLX-RIX", "MITE-DTX", "DEL", "DUP", "INV")
  }
    
  for (kernel_name in kernel_names) {
    K <- read.csv(here(
      "data", "kernels", glue("kernel_{kernel_name}.csv")
    ), row.names = 1) %>% 
      as.matrix()
    
    K_list[[kernel_name]] <- K
  }
}


# 3. Process data ------------------------------------------------------------------------

# Trait vector
y <- pheno[, trait, drop = TRUE]

# Train / test masks
train <- which(part_guide == "train" & ! is.na(y))
test <- which(part_guide == "test" & ! is.na(y))

# Scale features
if (model == "BayesC") {
  # Scale genotypes using all train / test data because genomic prediction problem
  X <- scale(X, center = TRUE, scale = TRUE)
}

# Scale continuous traits, just using data of training
if (var_type == "continuous") {
  y <- (y - mean(y[train])) / sd(y[train])
}


# Introduce NA in test accessions for prediction
y_NA <- y
y_NA[test] <- NA

# Construct linear predictor
ETA <- list()
if (model == "BayesC") {
  ETA[[1]] <- list(X = X,
                   model = model,
                   probIn = pi0,
                   counts = p0)

} else if (model == "RKHS") {
  for (kernel_name in kernel_names) {
    ETA[[kernel_name]] <- list(K = K_list[[kernel_name]],
                               model = model)
  }
}


# 4. Prediction --------------------------------------------------------------------------

cat("\n\n")
print(pars)
cat("\n\n")

t_start <- Sys.time()

fm <- BGLR(
  y = y_NA,
  response_type = response_type,
  ETA = ETA,
  nIter = n_iter,
  burnIn = burn_in,
  saveAt = here(
    project_folder, out_folder, ""
  ),
  verbose = FALSE
)

t_end <- Sys.time()
print(t_end - t_start)


# 5. Record results ----------------------------------------------------------------------

RES <- NULL

# Continuous variable
if (var_type == "continuous") {
  # Performance metrics
  cor <- cor(y[test], fm$yHat[test])
  MSE <- mean((y[test] - fm$yHat[test])**2)
  
  # Partition metrics
  size_train <- length(y[train])
  size_test <- length(y[test])
  
  RES <- data.frame(cor, MSE, size_train, size_test)
  
  # Plot predicted vs. observed
  pdf_path <- here(
    project_folder, out_folder, "cor_plot.pdf"
  )
  pdf(pdf_path)
  plot(fm$yHat[test] ~ y[test],
       xlab = "Observed (scaled)",
       ylab = "Predicted (scaled)",
       main = glue("{out_folder}, cor = {round(cor, 2)}"))
  abline(lm(fm$yHat[test] ~ y[test]))
  dev.off()

# Binary variable
} else if (var_type == "binary") {
  
  # If variable is 1 / 2 coded, recode to 0 / 1
  if (all(levels(factor(y[test])) == c("1", "2"))) {
    y_true <- recode(y[test], `2` = 1, `1` = 0)
  } else {
    y_true <- y[test]
  }
  
  # Calculate binary cross-entropy
  do_bin_crossentropy <- function(y_true, y_prob) {
    y_prob_corrected <- ifelse(y_true == 1, y_prob, 1 - y_prob)
    y_prob_corrected_log <- log(y_prob_corrected, base = exp(1))
    bin_ce <- sum(y_prob_corrected_log) / (- length(y_prob_corrected_log))
    return(bin_ce)
  }
  
  bin_ce <- do_bin_crossentropy(y_true = y_true, y_prob = fm$probs[test, 2])
  bin_ce
  
  # Partition metrics
  size_train <- length(y[train])
  size_test <- length(y[test])
  fraction_class2_train <- (table(y[train]) / length(y[train]))[2]
  fraction_class2_test <- (table(y[test]) / length(y[test]))[2]
  
  RES <- data.frame(bin_ce, size_train, size_test, fraction_class2_train, fraction_class2_test)
}

# Write
write_csv(RES, here(
  project_folder, out_folder, "results.csv"
))

# Save environment
save.image(here(
  project_folder, out_folder, "env.RData"
))
