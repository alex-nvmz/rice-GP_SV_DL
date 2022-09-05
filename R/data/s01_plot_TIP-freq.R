here::i_am(file.path("data", "01_accessions_phenotypes.R"))

library(here)
library(tidyverse)
library(glue)

# Compute allelic frequencies of TIPs (no MAF filtering)

# Load genotype matrices
dtx <- read_delim(here("data", "source_TIP-SNP-pheno", "DTX_matrix_IRIS.txt"))
mite <- read_delim(here("data", "source_TIP-SNP-pheno", "MITE_matrix_IRIS.txt"))
rix <- read_delim(here("data", "source_TIP-SNP-pheno", "RIX_matrix_IRIS.txt"))
rlx <- read_delim(here("data", "source_TIP-SNP-pheno","RLX_matrix_IRIS.txt"))


get_frequency_table <- function(geno, maf) {
  f1 <- apply(geno[, -1], 1, function(x) sum(x == 1) / sum(! is.na(x)))
  f0 <- 1 - f1
  MAF <- pmin(f1, f0)
  
  out <- data.frame(f1 = f1, f0 = f0, MAF = MAF, row.names = geno$chrom_start_end)
  return(out)
}


# Get allele frequency tables
dtx_f <- get_frequency_table(dtx)
mite_f <- get_frequency_table(mite)
rix_f <- get_frequency_table(rix)
rlx_f <- get_frequency_table(rlx)

# Plots f1
hist(dtx_f$f1,
     breaks = 40,
     xlab = "Frequency of presence allele",
     ylab = "Count",
     main = "DTX")

hist(mite_f$f1,
     breaks = 40,
     xlab = "Frequency of presence allele",
     ylab = "Count",
     main = "MITE")

hist(rix_f$f1,
     breaks = 40,
     xlab = "Frequency of presence allele",
     ylab = "Count",
     main = "RIX")

hist(rlx_f$f1,
     breaks = 40,
     xlab = "Frequency of presence allele",
     ylab = "Count",
     main = "RLX")

