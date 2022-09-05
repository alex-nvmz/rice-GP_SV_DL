here::i_am(file.path("data", "01_accessions_phenotypes.R"))

library(here)
library(tidyverse)
library(glue)
library(AGHmatrix)

# 1. Parameters --------------------------------------------------------------------------

out_folder <- "kernels"
marker_sets <- c("INV", "DUP", "DEL", "RLX-RIX", "MITE-DTX", "SNP")
# marker_sets <- c("all", "DUP-INV")

# 2. Load data ---------------------------------------------------------------------------

DATA <- lapply(marker_sets, function(x) {
  read_csv(here(
    "data", "main", glue("geno_{x}.csv")
  ))
})
names(DATA) <- marker_sets


# 3. Change format to matrix, as AGH package requires ------------------------------------

to_matrix <- function(data) {
  accs <- data$Accession
  data <- as.matrix(data[, -1])
  row.names(data) <- accs
  return(data)
}

DATA <- lapply(DATA, to_matrix)
DATA$INV[1:5, 1:10]


# 4. Recode binary markers ---------------------------------------------------------------

# Recode genotype matrices, 1 to 2 -> The ones which are not SNPs
# Because AGHmatrix is meant for SNPs coded as 0,1,2; so it calculates allelic frequencies
# regarding that

recode_binary_markers <- function(geno) {
  geno.recode <- apply(geno, 2, function(x) recode(as.numeric(x), `1` = 2))
  rownames(geno.recode) <- rownames(geno)
  return(geno.recode)
}

for (m_set in marker_sets) {
  # If geno matrix does not contain SNPs, recode
  if(! str_detect(m_set, "SNP")) {
    DATA[[m_set]] <- recode_binary_markers(DATA[[m_set]])
  }
}
DATA$INV[1:5, 1:10]

# For geno matrices mixed with SNPs ------------------------------------------------------
# Figure out which columns are SNPs, so as not to recode them
# ncol(SNP)
# identical(all[, 1:228871], SNP)
# identical(SNP.TIP[, 1:228871], SNP)
# 
# all[, - c(1:228871)] <- recode_binary_markers(all[, - c(1:228871)])
# SNP.TIP[, - c(1:228871)] <- recode_binary_markers(SNP.TIP[, - c(1:228871)])
# ----------------------------------------------------------------------------------------


# 5. Compute relationship matrices -------------------------------------------------------

DATA_Gm <- lapply(DATA, function(x) {
  Gmatrix(SNPmatrix = x, method = "VanRaden", missingValue = NA, maf = 0.01)
})

# 6. Write -------------------------------------------------------------------------------

for (m_set in marker_sets) {
  write.csv(DATA_Gm[[m_set]], here(
    "data", out_folder, glue("kernel_{m_set}.csv")
  ))
}
