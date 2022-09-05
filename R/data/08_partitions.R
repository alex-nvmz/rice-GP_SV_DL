here::i_am(file.path("data", "01_accessions_phenotypes.R"))

library(here)
library(tidyverse)
library(glue)

# 1. Parameters --------------------------------------------------------------------------
# Number of folds for k-fold cross validation
k <- 10
# Percentage of data to include for testing in each fold
p <- 1 / k


# 2. Load data ---------------------------------------------------------------------------

pheno <- read_csv(here("data", "main", "pheno_original.csv"))
head(pheno)



# 3. Create k-fold partitions ------------------------------------------------------------

# Results with indices of each partition 
res <- vector(mode = "list", length = k)
names(res) <- glue("Fold.{1:k}")

# All indices
indices <- seq(1:nrow(pheno))

# List with indices already used for test
test_list <- NULL

# For each fold
for (i in 1:k){
  # Indices to sample from not already used for test
  x <- indices[! indices %in% test_list]
  # Size of the sample (p % of the total)
  size <- round(length(indices) * p)
  
  # If population is smaller than sample size, take all population
  if (length(x) < size) {
    size <- length(x)
  }
  
  # Indices for each case
  test <- sample(x = x,
                 size = size,
                 replace = FALSE)
  training <- indices[- test]
  
  test_list <- c(test_list, test)
  
  # Append to list
  res[[i]] <- list(training = training, test = test)
}

# Check sizes of folds
lapply(res, function(x) lapply(x, length))


# 4. Generate results table --------------------------------------------------------------

# Final table with partitions
PARTS <- data.frame(Accession = pheno$Accession)

# For each fold
for (i in 1:length(res)) {
  # column vector
  col <- vector(mode = "character", length = length(PARTS$Accession))
  col[res[[i]]$test] <- "test"
  col[res[[i]]$training] <- "train"
  
  PARTS <- cbind(PARTS, col)
  colnames(PARTS)[1 + i] <- names(res[i])
}
PARTS

# Add AroAdm column
col <- vector(mode = "character", length = length(PARTS$Accession))
test <- which(pheno$Population %in% c("ARO", "ADM"))
col[test] <- "test"
col[- test] <- "train"

PARTS <- cbind(PARTS, AroAdm = col)

# Write
write_csv(PARTS, here("data", "utility", "partitions.csv"))


##########################################################################################
##########      MOVE TO "main" FOLDER TO AVOID RISK OF OVERWRITING    ####################
##########################################################################################


# 5. Generate table for IInd -------------------------------------------------------------

col <- vector(mode = "character", length = length(PARTS$Accession))
test <- which(pheno$Population == "IND" & pheno$Pedigree == "I")
col[test] <- "test"
col[- test] <- "train"

IIND <- data.frame(Accession = pheno$Accession, IInd = col)

write_csv(IIND, here("data", "utility", "partitions_IInd.csv"))




# ----------------------------------------------------------------------------------------

# Attempt with caret
# For some reason, folds are of different sizes

# library(caret)
# testIndex <- createFolds(pheno$Accession,
#                          k = 10,
#                          list = TRUE,
#                          returnTrain = FALSE)
# 
# lapply(testIndex, length)
# round(738 / 10)
# 
# lapply(testIndex, function(x){
#   test <- pheno[x, ]
#   train <- pheno[-x, ]
#   return(list(test = test, train = train))
# })

