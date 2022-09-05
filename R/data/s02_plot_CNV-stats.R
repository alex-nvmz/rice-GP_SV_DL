here::i_am(file.path("data", "01_accessions_phenotypes.R"))

library(here)
library(tidyverse)

# ----------------------------------------------------------------------------------------
# Plot missing rate of CNV markers and Accessions

miss.site <- read_delim(here("data", "output_CNV", "CNVnator_Q10_goodRD_noCN1-3_accs.recode.lmiss"))
miss.indv <- read_delim(here("data", "output_CNV", "CNVnator_Q10_goodRD_noCN1-3_accs.recode.imiss"))

# Per CNV marker
head(miss.site)
nrow(miss.site)

length(miss.site$F_MISS)
mean(miss.site$F_MISS)

hist(miss.site$F_MISS,
     breaks = 20,
     main = "CNVs (207,927 SVs)",
     xlab = "Missing rate per CNV",
     ylab = "Count")

# Per accession
head(miss.indv)
nrow(miss.indv)
# 661 samples from 738
hist(miss.indv$F_MISS,
     breaks = 20,
     main = "CNVs (661 accessions)",
     xlab = "Missing rate per accession",
     ylab = "Count")

