here::i_am(file.path("data", "01_accessions_phenotypes.R"))

library(here)
library(tidyverse)

# ----------------------------------------------------------------------------------------
# Produce genotype matrices of mixed marker types

# Load data
SNP <- read_csv(here("data", "main", "geno_SNP.csv"))

RIX <- read_csv(here("data", "main", "geno_RIX.csv"))
RLX <- read_csv(here("data", "main", "geno_RLX.csv"))
DTX <- read_csv(here("data", "main", "geno_DTX.csv"))
MITE <- read_csv(here("data", "main", "geno_MITE.csv"))

DEL <- read_csv(here("data", "main", "geno_DEL.csv"))
DUP <- read_csv(here("data", "main", "geno_DUP.csv"))
INV <- read_csv(here("data", "main", "geno_INV.csv"))

# Join datasets
all.markers <- bind_cols(SNP,
                         RIX[, -1],
                         RLX[, -1],
                         DTX[, -1],
                         MITE[, -1],
                         DEL[, -1],
                         DUP[, -1],
                         INV[, -1])

RLX.RIX <- bind_cols(RLX,
                     RIX[, -1])

MITE.DTX <- bind_cols(MITE,
                      DTX[, -1])

SNP.TIP <- bind_cols(SNP,
                     RLX[, -1],
                     RIX[, -1],
                     MITE[, -1],
                     DTX[, -1])

MITE.DTX.RLX.RIX.DEL <- bind_cols(MITE,
                                  DTX[, -1],
                                  RLX[, -1],
                                  RIX[, -1],
                                  DEL[, -1])

DUP.INV <- bind_cols(DUP,
                     INV[, -1])

# Write
write_csv(all.markers, here("data", "main", "geno_all.csv"))
write_csv(RLX.RIX, here("data", "main", "geno_RLX-RIX.csv"))
write_csv(MITE.DTX, here("data", "main", "geno_MITE-DTX.csv"))
write_csv(SNP.TIP, here("data", "main", "geno_SNP-TIP.csv"))
write_csv(MITE.DTX.RLX.RIX.DEL, here("data", "main", "geno_MITE-DTX-RLX-RIX-DEL.csv"))
write_csv(DUP.INV, here("data", "main", "geno_DUP-INV.csv"))
