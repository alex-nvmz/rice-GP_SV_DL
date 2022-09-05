here::i_am(file.path("data", "01_accessions_phenotypes.R"))

library(here)
library(tidyverse)


# 1. Load accs and original SV files -----------------------------------------------------

accs <- read_csv(here("data", "utility", "accessions.txt"))
accs <- accs$Accession

INS <- read_delim(here("data", "source_SV", "NB_INS_mergesam_clustered.txt"),
                  col_names = c("chrom", "start", "end", "type_length", "Accession", "cluster"))
DEL <- read_delim(here("data", "source_SV", "NB_DEL_mergesam_clustered.txt"),
                  col_names = c("chrom", "start", "end", "type_length", "Accession", "cluster"))
DUP <- read_delim(here("data", "source_SV", "NB_DUP_mergesam_clustered.txt"),
                  col_names = c("chrom", "start", "end", "type_length", "Accession", "cluster"))
INV <- read_delim(here("data", "source_SV", "NB_INV_mergesam_clustered.txt"),
                  col_names = c("chrom", "start", "end", "type_length", "Accession", "cluster"))


# 2. Check number of genotyped accessions ------------------------------------------------

head(INS)
length(unique(INS$Accession))
# 562 different accessions appear -> The ones mentioned in the paper
# Detected in 562 high-coverage samples.
INS <- filter(INS, Accession %in% accs)
INS_accs <- unique(INS$Accession)
length(unique(INS$Accession))
# 390 accessions for INS

head(DEL)
length(unique(DEL$Accession))
# 3023 accessions
DEL <- filter(DEL, Accession %in% accs)
DEL_accs <- unique(DEL$Accession)
length(unique(DEL$Accession))
# 738 (all) accessions for DEL

head(DUP)
length(unique(DUP$Accession))
# 3023 accessions
DUP <- filter(DUP, Accession %in% accs)
DUP_accs <- unique(DUP$Accession)
length(unique(DUP$Accession))
# 738 (all) accessions for DUP

head(INV)
length(unique(INV$Accession))
# 3023 accessions
INV <- filter(INV, Accession %in% accs)
INV_accs <- unique(INV$Accession)
# 738 (all) accessions for INV
length(unique(INV$Accession))

# 3. Tidy up 1 ---------------------------------------------------------------------------

# Make "chrom_start_end" field to sum up position
INS <- mutate(INS, "chrom_start_end" = paste(chrom, start, end, sep = "_"))
DEL <- mutate(DEL, "chrom_start_end" = paste(chrom, start, end, sep = "_"))
DUP <- mutate(DUP, "chrom_start_end" = paste(chrom, start, end, sep = "_"))
INV <- mutate(INV, "chrom_start_end" = paste(chrom, start, end, sep = "_"))

# Select relevant variables
INS <- INS %>% 
  select(cluster, Accession, chrom_start_end, type_length)
DEL <- DEL %>% 
  select(cluster, Accession, chrom_start_end, type_length)
DUP <- DUP %>% 
  select(cluster, Accession, chrom_start_end, type_length)
INV <- INV %>% 
  select(cluster, Accession, chrom_start_end, type_length)

# 4. Extract necessary information -------------------------------------------------------

# For each SV type, recode info by each cluster (unique SV) and compute relevant fields

# Insertions
INS_ <- INS %>% 
  arrange(cluster) %>% 
  # Group by unique SV (cluster)
  group_by(cluster) %>% 
  mutate(
    # Save number of calls of a given cluster/SV (times it has been detected in the different samples)
    n_calls = n(),
    # Save unique accession names
    accessions = paste0(unique(Accession), collapse = ","),
    # Save number of unique accessions (frequency of the SV)
    n_accessions = length(unique(Accession)),
    # Save positions at which the SV can be located
    positions = paste0(unique(chrom_start_end), collapse = ","),
    # Save the sequences which might be alleles (unique ones)
    # Not sure if 1 sequence per call, or several (separated by ",")
    
    # SEPARATE BY ";"
    maybe_alleles = paste0(unique(type_length), collapse = ";"),
    # Number of (maybe) alleles
    n_maybe_alleles = length(unique(type_length))
  ) %>%
  # Select variables
  select(cluster, n_accessions, accessions, positions, maybe_alleles, n_maybe_alleles, n_calls) %>% 
  # Remove repeated rows
  distinct() %>% 
  ungroup() %>% 
  # Create frequency of presence of SV
  mutate(frequency = n_accessions / length(INS_accs)) %>% 
  arrange(cluster)

head(INS_)
# Number of SVs
nrow(INS_)
length(INS_accs)
# Plot distribution of frequencies 
hist(INS_$frequency,
     breaks = 40,
     xlab = "Frequency of presence allele",
     ylab = "Count",
     main = "INS (351,333 SVs, 390 accessions)")

# More calls than accessions, I don't quite understand
plot(INS_$n_accessions, INS_$n_calls)
colSums(INS_[,2])
colSums(INS_[,7])

# Deletions
DEL_ <- DEL %>% 
  arrange(cluster) %>% 
  group_by(cluster) %>% 
  mutate(
    n_calls = n(),
    accessions = paste0(unique(Accession), collapse = ","),
    n_accessions = length(unique(Accession)),
    positions = paste0(unique(chrom_start_end), collapse = ","),
    # SEPARATE BY ","
    maybe_alleles = paste0(unique(type_length), collapse = ","),
    n_maybe_alleles = length(unique(type_length))
  ) %>% 
  select(cluster, n_accessions, accessions, positions, maybe_alleles, n_maybe_alleles, n_calls) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(frequency = n_accessions / length(DEL_accs)) %>% 
  arrange(cluster)


head(DEL_)
# Number of SVs
nrow(DEL_)
# Plot distribution of frequencies 
hist(DEL_$frequency,
     breaks = 40,
     xlab = "Frequency of presence allele",
     ylab = "Count",
     main = "DEL (466,236 SVs, 738 accessions)")

# Duplications
DUP_ <- DUP %>% 
  arrange(cluster) %>% 
  group_by(cluster) %>% 
  mutate(
    n_calls = n(),
    accessions = paste0(unique(Accession), collapse = ","),
    n_accessions = length(unique(Accession)),
    positions = paste0(unique(chrom_start_end), collapse = ","),
    # SEPARATE BY ","
    maybe_alleles = paste0(unique(type_length), collapse = ","),
    n_maybe_alleles = length(unique(type_length))
  ) %>% 
  select(cluster, n_accessions, accessions, positions, maybe_alleles, n_maybe_alleles, n_calls) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(frequency = n_accessions / length(DUP_accs)) %>% 
  arrange(cluster)


head(DUP_)
# Number of SVs
nrow(DUP_)
# Plot distribution of frequencies 
hist(DUP_$frequency,
     breaks = 40,
     xlab = "Frequency of presence allele",
     ylab = "Count",
     main = "DUP (51,693 SVs, 738 accessions)")

# Inversions
INV_ <- INV %>% 
  arrange(cluster) %>% 
  group_by(cluster) %>% 
  mutate(
    n_calls = n(),
    accessions = paste0(unique(Accession), collapse = ","),
    n_accessions = length(unique(Accession)),
    positions = paste0(unique(chrom_start_end), collapse = ","),
    # SEPARATE BY ","
    maybe_alleles = paste0(unique(type_length), collapse = ","),
    n_maybe_alleles = length(unique(type_length))
  ) %>% 
  select(cluster, n_accessions, accessions, positions, maybe_alleles, n_maybe_alleles, n_calls) %>% 
  distinct() %>% 
  ungroup() %>% 
  mutate(frequency = n_accessions / length(INV_accs)) %>% 
  arrange(cluster)


head(INV_)
# Number of SVs
nrow(INV_)
# Plot distribution of frequencies 
hist(INV_$frequency,
     breaks = 200,
     xlab = "Frequency of presence allele",
     ylab = "Count",
     main = "INV (116,190 SVs, 738 accessions)")


# 5. Write information files

# Full dataset
write_csv(INS_, here("data", "utility", "info_INS.csv"))
write_csv(DEL_, here("data", "utility", "info_DEL.csv"))
write_csv(DUP_, here("data", "utility", "info_DUP.csv"))
write_csv(INV_, here("data", "utility", "info_INV.csv"))

# Keep INS accessions
write_csv(data.frame(Accession = INS_accs), here("data", "utility", "accessions_INS.txt"))
