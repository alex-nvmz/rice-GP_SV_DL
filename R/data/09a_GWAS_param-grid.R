here::i_am(file.path("data", "01_accessions_phenotypes.R"))

library(here)
library(tidyverse)

# Generate params defining each GWAS study (trait and data partition)

pheno <- read_csv(here("data", "main", "pheno_original.csv"))
parts <- read_csv(here("data", "main", "partitions.csv"))

# Choose traits
traits <- colnames(pheno)[c(4, 10, 11, 13)]
# traits <- colnames(pheno)[c(4, 8, 10, 13)]
# traits <- colnames(pheno)[c(11)]

par_grid <- expand.grid(
  trait = traits,
  part = colnames(parts)[-1]
)

par_grid
write_csv(par_grid, here("data", "utility", "pars_GWAS.csv"))
