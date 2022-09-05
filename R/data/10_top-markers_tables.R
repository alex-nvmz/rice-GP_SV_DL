here::i_am(file.path("data", "01_accessions_phenotypes.R"))

library(here)
library(tidyverse)
library(glue)


# 1. Parameters --------------------------------------------------------------------------

# HARD-CODED
n_mark <- 10000
out_folder <- "top-markers"
in_folder <- "output_GWAS"

# VARIABLE
pars <- read_csv(here("data", "utility", "pars_GWAS.csv"))
# Marker sets to get top markers from
top_markers <- c("all", "SNP")
pars <- expand.grid(trait = unique(pars$trait),
                    part = unique(pars$part),
                    top_marker = top_markers)


# 2. Get filtered P-value tables ---------------------------------------------------------

for (i in 1:nrow(pars)){
  print(pars[i, ])
  # Extract parameters
  trait <- pars[[i, "trait"]]
  part <- pars[[i, "part"]]
  top_marker <- pars[[i, "top_marker"]]
  
  # Load P-value table
  table <- read_csv(here(
    "data", in_folder, glue("Pvals_{part}_{trait}.csv")
  ))
  
  # Filter for a given marker set if necessary
  if (top_marker != "all") {
    table <- filter(table, marker_set == top_marker)
  }
  
  # Table with top markers
  table_top <- table %>% 
    arrange(Pvalue) %>%   # Sort by P-value
    slice(1:n_mark) %>%   # Take top n
    mutate(               # Create position field
      pos = {
        marker %>% 
          str_split("_") %>% 
          map(function(x) {
            x[-1] %>% 
              paste(collapse = "_")
          }) %>% 
          as.character()
      }
    ) %>% 
    relocate(pos) %>% 
    # Sort markers by position (create factor with custom order of levels)
    mutate(pos = factor(pos, levels = str_sort(unique(pos), numeric = TRUE))) %>% 
    arrange(pos)
  
  # Write
  write_csv(table_top, here(
    "data", out_folder, glue("Pvals_top-{n_mark}-{top_marker}_{part}_{trait}.csv")
  ))
}

