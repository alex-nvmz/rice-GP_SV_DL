here::i_am(file.path("data", "01_accessions_phenotypes.R"))

library(here)
library(tidyverse)
library(glue)

# Transfrom SNPs from PLINK format to genotype matrix

# 1. Load data ---------------------------------------------------------------------------
# Accessions
accs <- read_csv(here("data", "utility", "accessions.txt"))
accs <- accs$Accession

# Genotype matrix from .bed file
library(BEDMatrix)
geno.snp <- BEDMatrix(here("data", "source_TIP-SNP-pheno", "snps.bed"), simple_names = TRUE)
geno.snp[1:5, 1:5]

# SNP info
info.snp <- read.table(here("data", "source_TIP-SNP-pheno", "snps.bim"), sep = "\t")


# 2. Tidy up and checking ----------------------------------------------------------------

# Get SNP names coded as VarType_chrX_position
head(info.snp)
info.snp <- mutate(info.snp, chr_position = glue("SNP_chr{sprintf('%02d', V1)}_{V4}"))
colnames(geno.snp) <- info.snp$chr_position

# Check accessions in same order as reference
identical(rownames(geno.snp), accs)

# Check same matrix as Ioanna
load(here("data", "source_TIP-SNP-pheno", "snps_matrix_Ioanna.RData"))
SNP.i <- snps
rm(snps)
identical(unname(as.matrix(geno.snp)), unname(SNP.i))


# 3. Check if data was filtered by MAF ---------------------------------------------------
maf <- 0.01

f1 <- apply(geno.snp, 2, function(x) sum(x) / (sum(! is.na(x)) * 2))
f0 <- 1 - f1
MAF <- pmin(f1, f0)
MAF.filt <- MAF[MAF > maf]
geno.filt <- geno.snp[, names(MAF.filt)]

identical(dim(geno.filt), dim(geno.snp))


# 4. Final tidy up -----------------------------------------------------------------------

geno.snp.final <- geno.snp %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Accession") %>% 
  as_tibble()
geno.snp.final[1:5, 1:5]

# Markers already sorted by position
identical(colnames(geno.snp.final[-1]),
          str_sort(colnames(geno.snp.final)[-1], numeric = TRUE))

write_csv(geno.snp.final, here("data", "main", "geno_SNP.csv"))
