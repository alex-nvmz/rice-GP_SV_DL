here::i_am(file.path("data", "01_accessions_phenotypes.R"))

library(here)
library(tidyverse)
library(glue)
library(ggfortify)
library(cowplot)

# 1. Parameters --------------------------------------------------------------------------

marker_sets <- c("INV", "DUP", "DEL", "RLX-RIX", "MITE-DTX", "SNP")
out_folder <- "plots"


# 2. Load genotype matrices into list ----------------------------------------------------

DATA <- lapply(marker_sets, function(x) {
  read_csv(here(
    "data", "main", glue("geno_{x}.csv")
  ))
})
names(DATA) <- marker_sets

# Add insertions
INS <- read_csv(here("data", "utility", "geno_INS.csv"))
INS <- list("INS" = INS)
DATA <- append(DATA, INS)
rm(INS)

# 3. Integrate population data -----------------------------------------------------------

# Load population data
pop <- read_csv(here(
  "data", "main", "pheno_original.csv"
  ), col_select = c("Accession", "Population", "Pedigree"))


# Add population info as variable
add_population <- function(geno, pop) {
  geno %>% 
    left_join(pop[, c("Accession", "Population")], by = "Accession") %>% 
    relocate(Accession, Population)
}

DATA <- lapply(DATA, function(x) {
  add_population(x, pop)
})


# 4. Do PCA and save plots + PCA results -------------------------------------------------
PCA <- list()

# For every marker set
for (i in 1:length(DATA)) {
  data_name <- names(DATA[i])
  
  res_pca <- prcomp(DATA[[i]][, - c(1, 2)], scale. = TRUE)
  # Save PCA into list
  PCA[[data_name]] <- res_pca
  
  # Plot
  pdf(here(
    "data", "plots", glue("PCA_{data_name}.pdf")
  ))
  print(
    autoplot(res_pca, data = DATA[[i]], colour = "Population") +
      theme_bw() +
      labs(title = glue("{data_name}"))
  )
  dev.off()
  
}

# Save R object
saveRDS(PCA, here(
  "data", "utility", "PCA_list.rds"
))





# ----------------------------------------------------------------------------------------

# Make fancy plot


marker_sets <- c("SNP", "DEL", "MITE-DTX", "RLX-RIX", "DUP", "INV")
out_folder <- "plots"


# Load data
PCA <- readRDS(here(
  "data", "utility", "PCA_list.rds"
))

pop <- read_csv(here(
  "data", "main", "pheno_original.csv"
), col_select = c("Accession", "Population", "Pedigree"))

# Get plots ----
# plot_PCA <- function(marker_set, pca_list, pop_data) {
#   autoplot(PCA[[marker_set]], data = pop_data, colour = "Population") +
#     theme_bw() +
#     labs(title = marker_set,
#          colour = "Group")
# }
# 
# plot_PCA(marker_set = "SNP",
#          pca_list = PCA,
#          pop_data = pop)
# ----

plot_list <- list()

for (marker_set in marker_sets) {
  print(marker_set)
  
  tmp <- autoplot(PCA[[marker_set]], data = pop, colour = "Population", alpha = 0.5) +
    theme_bw() +
    labs(title = marker_set,
         colour = "Group")
  
  plot_list[[marker_set]] <- tmp
  
  rm(tmp)
}

length(plot_list)
rm(PCA)


grid <- plot_grid(plotlist = plot_list, ncol = 2, labels = "AUTO")
grid

png(here(
  "data", out_folder, "PCA_all.png"
), units = "in", width = 10, height = 11, res = 300)
print(grid)
dev.off()
