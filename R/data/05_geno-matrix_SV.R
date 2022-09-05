here::i_am(file.path("data", "01_accessions_phenotypes.R"))

library(here)
library(tidyverse)
library(glue)

# 1. Parameters --------------------------------------------------------------------------
# MAF to filter by
maf <- 0.01


# 2. Load data ---------------------------------------------------------------------------

INS <- read_delim(here("data", "utility", "info_INS.csv"))
DEL <- read_delim(here("data", "utility", "info_DEL.csv"))
DUP <- read_delim(here("data", "utility", "info_DUP.csv"))
INV <- read_delim(here("data", "utility", "info_INV.csv"))

# All 738 accessions
accs <- read_table(here("data", "utility", "accessions.txt"))
accs <- accs$Accession

# INS in 390 accessions
INS.accs <- read_table(here("data", "utility", "accessions_INS.txt"))
INS.accs <- INS.accs$Accession
# Sort by the order of 738 accessions
INS.accs <- INS.accs[order(match(INS.accs, accs))]


# 3. Filter by MAF -----------------------------------------------------------------------

filter_by_MAF <- function(DATA, maf) {
  DATA.filt <- DATA %>% 
    mutate(f0 = 1- frequency) %>% 
    mutate(MAF = pmin(frequency, f0)) %>% 
    filter(MAF > maf)
  return(DATA.filt)
}

INS.filt <- filter_by_MAF(INS, maf = maf)
DEL.filt <- filter_by_MAF(DEL, maf = maf)
DUP.filt <- filter_by_MAF(DUP, maf = maf)
INV.filt <- filter_by_MAF(INV, maf = maf)

# Count how many remaining
paste(nrow(INS), nrow(INS.filt))
paste(nrow(DEL), nrow(DEL.filt))
paste(nrow(DUP), nrow(DUP.filt))
paste(nrow(INV), nrow(INV.filt))


# 4. Save info tables of SVs filtered by MAF ---------------------------------------------

write_csv(INS.filt, here("data", "utility", "info_INS-filtered.csv"))
write_csv(DEL.filt, here("data", "utility", "info_DEL-filtered.csv"))
write_csv(DUP.filt, here("data", "utility", "info_DUP-filtered.csv"))
write_csv(INV.filt, here("data", "utility", "info_INV-filtered.csv"))


# 5. Obtain single identifier for SV position --------------------------------------------

# Get minimum and maximum position from the possible positions
get_single_positions <- function(DATA) {
  out <- DATA %>% 
    mutate(
      positions_vector = str_split(positions, ",")
    ) %>% 
    mutate(
      chrs = map(positions_vector,
                 function(x) {
                   x %>% 
                     str_split("_") %>% 
                     map(function(y) y[1])
                 }),
      poss = map(positions_vector,
                 function(x) {
                   x %>% 
                     str_split("_") %>% 
                     map(function(y) y[c(2,3)]) %>% 
                     as_vector()
                 })
    ) %>% 
    rowwise() %>% 
    mutate(
      chr_start_end = paste(chrs[1], min(poss), max(poss), sep = "_")
    ) %>% 
    pull(chr_start_end)
  return(out)
}

chr_start_end.ins <- get_single_positions(INS.filt)
chr_start_end.del <- get_single_positions(DEL.filt)
chr_start_end.dup <- get_single_positions(DUP.filt)
chr_start_end.inv <- get_single_positions(INV.filt)

# Add marker type to the code
chr_start_end.ins <- paste0("INS_", chr_start_end.ins)
chr_start_end.del <- paste0("DEL_", chr_start_end.del)
chr_start_end.dup <- paste0("DUP_", chr_start_end.dup)
chr_start_end.inv <- paste0("INV_", chr_start_end.inv)


# 6. Initialize genotype matrices --------------------------------------------------------

initialize_genotype_matrix <- function(accessions, positions) {
  m.fill <- matrix(nrow = length(accessions), ncol = length(positions))
  colnames(m.fill) <- positions
  row.names(m.fill) <- accessions
  return(m.fill)
}

ins <- initialize_genotype_matrix(INS.accs, chr_start_end.ins)
del <- initialize_genotype_matrix(accs, chr_start_end.del)
dup <- initialize_genotype_matrix(accs, chr_start_end.dup)
inv <- initialize_genotype_matrix(accs, chr_start_end.inv)


# 7. Fill genotype matrices (biallelic) --------------------------------------------------

fill_geno_matrix <- function(m.fill, source) {
  # For every column (SV) in the dataframe to fill
  for (i in 1:ncol(m.fill)) {
    print(colnames(m.fill)[i])
    # Save accessions string from the source, turn into vector
    accs.vec <- str_split(source[i, "accessions"], ",")[[1]]
    # If row in accs.vec, return 1. Insert in column
    m.fill[, i] <- ifelse(row.names(m.fill) %in% accs.vec, 1, 0)
  }
  return(m.fill)
}

ins <- fill_geno_matrix(ins, INS.filt)
del <- fill_geno_matrix(del, DEL.filt)
dup <- fill_geno_matrix(dup, DUP.filt)
inv <- fill_geno_matrix(inv, INV.filt)


# 8. Tidy up -----------------------------------------------------------------------------

format_geno_matrix <- function(m.geno) {
  m.geno <- m.geno %>% 
    data.frame() %>% 
    rownames_to_column(var = "Accession") %>% 
    as_tibble()
  return(m.geno)
}

ins <- format_geno_matrix(ins)
del <- format_geno_matrix(del)
dup <- format_geno_matrix(dup)
inv <- format_geno_matrix(inv)

# 9. Sort markers by position ------------------------------------------------------------

ins <- ins[, c("Accession", str_sort(colnames(ins)[-1], numeric = TRUE))]
del <- del[, c("Accession", str_sort(colnames(del)[-1], numeric = TRUE))]
dup <- dup[, c("Accession", str_sort(colnames(dup)[-1], numeric = TRUE))]
inv <- inv[, c("Accession", str_sort(colnames(inv)[-1], numeric = TRUE))]


# 10. Save --------------------------------------------------------------------------------

write_csv(ins, here("data", "utility", "geno_INS.csv")) # INS to utility

write_csv(del, here("data", "main", "geno_DEL.csv"))
write_csv(dup, here("data", "main", "geno_DUP.csv"))
write_csv(inv, here("data", "main", "geno_INV.csv"))


# ----------------------------------------------------------------------------------------
# Modification check (now they are sorted differently)
inv_final <- read_csv(here("data", "main", "geno_INV.csv"))
colnames(inv_final) <- 1:ncol(inv_final)

inv_previous_script <- read_csv("/home/anavartinez/CRAG_not_tracked/riceDL_main/alex/data/Final/geno_INV.csv")
colnames(inv_previous_script) <- 1:ncol(inv_previous_script)

identical(inv_final[, -1], inv_previous_script[, -1])
all.equal(inv_final[, -1], inv_previous_script[, -1])
