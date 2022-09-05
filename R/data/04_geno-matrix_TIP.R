here::i_am(file.path("data", "01_accessions_phenotypes.R"))

library(here)
library(tidyverse)
library(glue)

# 1. Parameters --------------------------------------------------------------------------
# MAF to filter by
maf <- 0.01


# 2. Load data ---------------------------------------------------------------------------

DTX <- read_delim(here("data", "source_TIP-SNP-pheno", "DTX_matrix_IRIS.txt"))
MITE <- read_delim(here("data", "source_TIP-SNP-pheno", "MITE_matrix_IRIS.txt"))
RIX <- read_delim(here("data", "source_TIP-SNP-pheno", "RIX_matrix_IRIS.txt"))
RLX <- read_delim(here("data", "source_TIP-SNP-pheno", "RLX_matrix_IRIS.txt"))


# 3. Modify ------------------------------------------------------------------------------

# Transpose original TIP matrices to have Accessions in the rows

transpose_data <- function(DATA) {
  return(
    DATA %>% 
      column_to_rownames(var = "chrom_start_end") %>% 
      t() %>%
      as.data.frame() %>% 
      rownames_to_column(var = "Accession") %>% 
      as_tibble()
  )
}

DTX <- transpose_data(DTX)
MITE <- transpose_data(MITE)
RIX <- transpose_data(RIX)
RLX <- transpose_data(RLX)

# Fix coding of accessions
# Include marker type in marker IDs
DTX <- mutate(DTX, Accession = str_replace(
  str_replace_all(Accession, "_", "-"), "-", "_"
  ))
colnames(DTX)[-1] <- glue("DTX_{colnames(DTX)[-1]}")

MITE <- mutate(MITE, Accession = str_replace(
  str_replace_all(Accession, "_", "-"), "-", "_"
  ))
colnames(MITE)[-1] <- glue("MITE_{colnames(MITE)[-1]}")

RIX <- mutate(RIX, Accession = str_replace(
  str_replace_all(Accession, "_", "-"), "-", "_"
  ))
colnames(RIX)[-1] <- glue("RIX_{colnames(RIX)[-1]}")

RLX <- mutate(RLX, Accession = str_replace(
  str_replace_all(Accession, "_", "-"), "-", "_"
  ))
colnames(RLX)[-1] <- glue("RLX_{colnames(RLX)[-1]}")

# Check order of accessions
accs <- read_csv(here("data", "utility", "accessions.txt"))
accs <- accs$Accession
identical(DTX$Accession, accs)


# 4. Filter by MAF -----------------------------------------------------------------------

filter_binary_by_MAF <- function(geno, maf) {
  f1 <- apply(geno[, -1], 2, function(x) sum(x == 1) / sum(! is.na(x)))
  f0 <- 1 - f1
  MAF <- pmin(f1, f0)
  MAF.filt <- MAF[MAF > maf]
  geno.filt <- geno[, c("Accession", names(MAF.filt))]
  return(geno.filt)
}

DTX.filt <- filter_binary_by_MAF(DTX, maf = maf)
MITE.filt <- filter_binary_by_MAF(MITE, maf = maf)
RIX.filt <- filter_binary_by_MAF(RIX, maf = maf)
RLX.filt <- filter_binary_by_MAF(RLX, maf = maf)

# Count how many remaining
paste(ncol(DTX[, -1]), ncol(DTX.filt[, -1]))
paste(ncol(MITE[, -1]), ncol(MITE.filt[, -1]))
paste(ncol(RIX[, -1]), ncol(RIX.filt[, -1]))
paste(ncol(RLX[, -1]), ncol(RLX.filt[, -1]))
print("")
ncol(RIX.filt[, -1]) + ncol(RLX.filt[, -1])
ncol(MITE.filt[, -1]) + ncol(DTX.filt[, -1])

# Markers already sorted by position
identical(colnames(DTX.filt)[-1], str_sort(colnames(DTX.filt)[-1], numeric = TRUE))
identical(colnames(MITE.filt)[-1], str_sort(colnames(MITE.filt)[-1], numeric = TRUE))
identical(colnames(RIX.filt)[-1], str_sort(colnames(RIX.filt)[-1], numeric = TRUE))
identical(colnames(RLX.filt)[-1], str_sort(colnames(RLX.filt)[-1], numeric = TRUE))

# 5. Write matrices ----------------------------------------------------------------------

write_csv(DTX.filt, here("data", "main", "geno_DTX.csv"))
write_csv(MITE.filt, here("data", "main", "geno_MITE.csv"))
write_csv(RIX.filt, here("data", "main", "geno_RIX.csv"))
write_csv(RLX.filt, here("data", "main", "geno_RLX.csv"))
