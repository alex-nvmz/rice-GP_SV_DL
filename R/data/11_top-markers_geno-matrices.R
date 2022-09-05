here::i_am(file.path("data", "01_accessions_phenotypes.R"))

library(here)
library(tidyverse)
library(glue)


# 1. Parameters --------------------------------------------------------------------------

# HARD-CODED
n_mark <- 10000
out_folder <- "top-markers"
in_folder <- "top-markers"

# VARIABLE
pars <- read_csv(here("data", "utility", "pars_GWAS.csv"))
# Marker sets to get top markers from
top_markers <- c("all", "SNP")
pars <- expand.grid(trait = unique(pars$trait),
                    part = unique(pars$part),
                    top_marker = top_markers)


# 2. Load genotype matrices list ---------------------------------------------------------

# Marker data in 1 list
marker_sets <- c("INV", "DUP", "DEL", "RLX-RIX", "MITE-DTX", "SNP")
rds_path <- here(
  "data", "utility", glue("markers-list_{paste(marker_sets, collapse = '_')}.rds")
)

# If no R object, create; else load
if (! file.exists(rds_path)) {
  markers <- lapply(marker_sets, function(x){
    read_csv(here("data", "main", glue("geno_{x}.csv")))
  })
  names(markers) <- marker_sets
  saveRDS(markers, rds_path)
  
} else {
  markers <- readRDS(rds_path)
}

# 3. Obtain genotype matrix for each partition, trait and top-marker_set -----------------


for (i in 1:nrow(pars)){
  print(pars[i, ])
  # Extract parameters
  trait <- pars[[i, "trait"]]
  part <- pars[[i, "part"]]
  top_marker <- pars[[i, "top_marker"]]
  
  # Load P-value table
  table <- read_csv(here(
    "data", in_folder, glue("Pvals_top-{n_mark}-{top_marker}_{part}_{trait}.csv")
  ))
  
  # Resulting genotype matrix
  RESULTS <- markers[[1]][, "Accession"]
  
  # For every marker set
  for (j in 1:length(markers)) {
    # Get geno matrix
    geno <- markers[[j]]
    # Select columns in P-value table
    RESULTS <- bind_cols(RESULTS,
                         geno[, which(colnames(geno) %in% table$marker)])
  }
  
  # Sort mixed geno matrix by marker positions
  cols <- colnames(RESULTS)[-1]
  pos <-
    str_split(cols, "_") %>% 
    map(function(x) {
      x[-1] %>% 
        paste(collapse = "_")
    }) %>% 
    unlist()
  # Columns in the sorted order
  cols <- cols[match(str_sort(pos, numeric = TRUE), pos)]
  
  RESULTS <- RESULTS[, c("Accession", cols)]
  
  # Write
  write_csv(RESULTS, here(
    "data", out_folder, glue("geno_top-{n_mark}-{top_marker}_{part}_{trait}.csv")
  ))
}
