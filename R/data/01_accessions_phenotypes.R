# Make sure you are in the correct project
here::i_am(file.path("data", "01_accessions_phenotypes.R"))
# To be able to use relative paths from the project's root (where .Rproj is located)
library(here)
library(tidyverse)


# 1. Load files --------------------------------------------------------------------------

pheno <- read_csv(here("data", "source_TIP-SNP-pheno", "all.phenotypes.csv"),
                  col_select = -1)
pedigree <- read_csv(here("data", "source_TIP-SNP-pheno", "iris_pedigree.csv"))


# 2. Check accessions of both files are in the same order --------------------------------

pheno <- pheno %>% 
  rename("Accession" = Variety) %>% 
  mutate(Accession = str_replace(Accession, " ", "_"))

pedigree <- pedigree %>% 
  rename("Accession" = X3K_DNA_IRIS_UNIQUE_ID) %>% 
  mutate(Accession = str_replace(
    str_replace_all(Accession, "_", "-"), "-", "_"
  )) 

identical(pheno$Accession, pedigree$Accession)

# Record accession list
accs <- pheno[, "Accession"]
write_csv(accs, here("data", "utility", "accessions.txt"))


# 3. Tidy data ---------------------------------------------------------------------------

pedigree <- pedigree %>% 
  rename("Pedigree" = Status_with_Pedigree.plus.knowledge)
pedigree_short <- pedigree[, c("Accession", "Pedigree")]

pheno <- pheno %>% 
  rename("Population" = SNPsubsp,
         "culm.diameter.1st.internode" = Culm.diameter...1st.internode,
         "grain.weight" = Grain.weight,
         "salt.injury.at.EC12" = Salt.injury.at.EC12,
         "time.to.flowering.from.sowing" = Time.to.flowering..from.sowing)

pheno <- pheno %>% 
  left_join(pedigree_short, by = "Accession") %>% 
  relocate(Accession, Population, Pedigree)

# File with all info about accessions
write_csv(pheno, here("data", "main", "pheno_original_extra-info.csv"))  

pheno_short <- pheno %>% 
  select(- c(Centre, K9.group, country, region, readepth, Filtered.insertions))

# File with most relevant fields
write_csv(pheno_short, here("data", "main", "pheno_original.csv"))

