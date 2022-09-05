here::i_am(file.path("data", "01_accessions_phenotypes.R"))

library(here)
library(tidyverse)
library(glue)

# 1. Parameters --------------------------------------------------------------------------

marker_sets <- c("INV", "DUP", "DEL", "RLX-RIX", "MITE-DTX", "SNP")
out_folder <- "kernels"


# 2. Load kernels into list --------------------------------------------------------------

DATA <- lapply(marker_sets, function(x) {
  read.csv(here(
    "data", "kernels", glue("kernel_{x}.csv")
  ), row.names = 1)
})
names(DATA) <- marker_sets


# 3. Check kernels are positive definite and do PCA --------------------------------------

for (i in 1:length(DATA)) {
  data_name <- names(DATA[i])
  print(data_name)
  
  kernel <- DATA[[i]]
  kernel_scaled <- scale(kernel, center = TRUE, scale = TRUE)
  
  # Eigendecomposition
  res_eigen <- eigen(kernel_scaled)
  values <- res_eigen$values
  vectors <- res_eigen$vectors
  
  # Check if positive definite
  print(all(values >= 0))
  print(min(values))
  cat("\n")
  # DEL, MITE-DTX and SNP are not positive definite
  
  # Add small constant (0.0001) to covariance matrix and repeat
  diag(kernel_scaled) <- diag(kernel_scaled) * (1 + 1e-4)
  
  res_eigen <- eigen(kernel_scaled)
  values <- res_eigen$values
  vectors <- res_eigen$vectors
  
  print(all(values >= 0))
  print(min(values))
  cat("\n\n")
  
  # Prepare artificial neural network input
  ann_input <- vectors
  ann_input <- as.data.frame(ann_input)
  colnames(ann_input) <- glue("PC{1:ncol(ann_input)}")
  row.names(ann_input) <- row.names(kernel_scaled)
  ann_input <- rownames_to_column(ann_input, "Accession")
  
  # Save
  write_csv(ann_input, here(
    "data", "kernels", glue("PC_{data_name}.csv")
  ))
}
