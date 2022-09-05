
library(tidyverse)
library(here)

data <- read_csv(here(
  "data", "top-markers", "geno_top-10000-all_AroAdm_grain.weight.csv"
))


data_sub <- data[, 1:1001]

write_csv(data_sub, here(
  "DL", "data", "geno_test.csv"
))
