here::i_am(file.path("data", "01_accessions_phenotypes.R"))

library(here)
library(tidyverse)
library(glue)

# Count how many accessions are not NA for different trait combinations ------------------

pheno <- read_csv(here("data", "main", "pheno_original.csv"))
traits <- colnames(pheno)[c(4, 10, 11, 13)]

combn_size <- 1:length(traits)

# For each combination size
for (i in combn_size){
  # Combinations list
  combinations <- combn(traits, m = i, simplify = FALSE)
  # Results dataframe
  RES <- NULL
  
  # For each combination
  for (j in 1:length(combinations)) {
    comb <- combinations[[j]]
    # Final boolean vector
    bool_res <- NULL
    
    # For each trait
    for (k in 1:length(comb)) {
      trait <- comb[[k]]
      print(trait)
      
      # Not NA
      y <- pheno[trait]
      bool <- !is.na(y)
      
      if (k == 1) {
        bool_res <- bool
      } else {
        bool_res <- bool_res & bool
      }
    
      # size <- sum(bool)
      # print(size)
    }
    # Count
    size <- sum(bool_res)
    print(size)
    
    # Save results
    col <- size
    RES <- cbind(RES, col)
    col_name <- paste(comb, collapse = "_")
    colnames(RES)[j] <- col_name
    RES <- as.data.frame(RES)
    
    cat("\n")
  }
  
  # Write for this combination size
  write_csv(RES, here(
    "data", "utility", glue("not-NA-count_comb-{i}-traits.csv")
  ))
}
