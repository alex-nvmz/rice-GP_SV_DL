here::i_am(file.path("data", "01_accessions_phenotypes.R"))

library(here)
library(tidyverse)
library(glue)

# 1. Parameters --------------------------------------------------------------------------

# HARD-CODED ----

# Number of PCs to correct by for population structure in regression
n_PC <- 4
# Separation with "-" is important
marker_sets <- c("INV", "DUP", "DEL", "RLX-RIX", "MITE-DTX", "SNP")

out_folder <- "output_GWAS"

# VARIABLE ----

# ----
# Read from parameter file
par_file <- "pars_GWAS.csv"
par_table <- read_csv(here("data", "utility", par_file))
# Get row number from command line
args <- commandArgs(trailingOnly = T)
par_row <- as.numeric(args[1])
pars <- par_table[par_row, ]
# ----

# Extract pars
trait <- pars[["trait"]]
part <- pars[["part"]]

# Variable type
var_table <- tribble(
  ~trait, ~var_type,
  "culm.diameter.1st.internode", "binary",
  "culm.strength", "binary",
  "flag.leaf.angle", "categorical",
  "grain.length", "continuous",
  "grain.width", "continuous",
  "leaf.length", "categorical",
  "leaf.senescence", "binary",
  "grain.weight", "continuous",
  "salt.injury.at.EC12", "binary",
  "time.to.flowering.from.sowing", "continuous",
  "panicle.threshability", "binary"
)

var_type <- var_table[which(var_table$trait == trait), "var_type", drop = TRUE]

# 2. Load data ---------------------------------------------------------------------------

pheno <- read_csv(here("data", "main", "pheno_original.csv"))
parts <- read_csv(here("data", "main", "partitions.csv"))

# Marker data in 1 list
rds_path <- here(
  "data", "utility", glue("markers-list_{paste(marker_sets, collapse = '_')}.rds")
)

# If no R object, create; else load
if (! file.exists(rds_path)) {
  markers <- lapply(marker_sets, function(x){
    read_csv(here("data", "main", glue("geno_{x}.csv")))
  })
  names(markers) <- marker_sets
  saveRDS(markers, rds_path)
  
} else {
  markers <- readRDS(rds_path)
}


# 3. GWAS --------------------------------------------------------------------------------

# Mask to select training data
train_mask <- parts[, part, drop = TRUE] == "train"

# Extract target variable
y <- pheno[, trait, drop = TRUE]
y_sub <- y[train_mask]

# PDF for plots
pdf_path <- here("data", out_folder, glue("GWAS_{part}_{trait}.pdf"))
pdf(pdf_path, onefile = TRUE)
par(mfrow = c(2,1))

# Final table with each marker and Pvalue
RESULTS <- NULL

cat("\n\n\n")
print(pars)
cat("\n")

t_start_global <- Sys.time()
# For every marker set
for (i in 1:length(markers)) {
  print(names(markers)[i])
  t_start <- Sys.time()
  
  # Extract regressor variables
  geno <- markers[[i]]
  geno_sub <- geno[train_mask, -1]
  rm(geno)
  # Make sure regressors don't have sd = 0 after applying mask, because can't scale or PCA
  geno_sub <- geno_sub[, which(apply(geno_sub, 2, sd) != 0)]
  
  # Initialize temporary dataframe for marker set
  TMP <- data.frame(
    marker = colnames(geno_sub),
    marker_set = names(markers)[i],
    Pvalue = NA
  )
  
  # Extract PCs to correct for population structure in regression
  res.pca <- prcomp(geno_sub, center = TRUE, scale = TRUE)
  PC <- res.pca$x[, 1:n_PC]
  
  # For each marker
  for (j in 1:ncol(geno_sub)) {
    # Extract marker
    geno_sub_marker <- geno_sub[, j]
    
    # Scale predictor variables: faster convergence
    geno_sub_marker <- scale(geno_sub_marker)
    PC <- scale(PC)
    
    # Regression for continuous or binary variables
    if (var_type == "continuous") {
      # Scale only continuous targets
      y_sub <- scale(y_sub)
      # Linear regression
      m <- summary(lm(y_sub ~ geno_sub_marker + PC))
      
      if ("geno_sub_marker" %in% rownames(m$coefficients)) {
        pvalue <- m$coefficients["geno_sub_marker", "Pr(>|t|)"]
      } else {
        pvalue <- NA
      }
    } else if (var_type == "binary") {
      # Logistic regression
      m <- summary(glm(as.factor(y_sub) ~ geno_sub_marker + PC, family = "binomial"))
      
      if ("geno_sub_marker" %in% rownames(m$coefficients)) {
        pvalue <- m$coefficients["geno_sub_marker", "Pr(>|z|)"]
      } else {
        pvalue <- NA
      }
    }
    # Add Pvalue to tmp results
    TMP[j, "Pvalue"] <- pvalue
  }
  # Add data of marker set to global data
  RESULTS <- rbind(RESULTS, TMP)
  
  # Plots ----
  
  # Manhattan plot
  # Significance threshold of 0.05 with Bonferroni correction
  plot(-log10(TMP$Pvalue),
       col = "grey", cex = 0.6, pch = 19,
       main = c(
        glue("{names(markers)[i]}_{trait}"),
        glue("Pvalues < {0.05 / ncol(geno_sub)}:"),
        length(which(TMP$Pvalue < 0.05 / ncol(geno_sub)))
      )
    )
  abline(h = -log10(0.05 / ncol(geno_sub)), lty = 2, col = "red")
  
  # QQ-plot
  qqplot(rexp(length(TMP$Pvalue), rate = log(10)),
         -log10(TMP$Pvalue),
         xlab = "Expected quantile",
         pch = 19,
         cex = 0.4,
         main = c(
           glue("{names(markers)[i]}_{trait}"),
           "Pvalues < 0.05:",
           length(which(TMP$Pvalue < 0.05))
         )
        )
  abline(0, 1, col = "red")
  # ---------
  
  t_end <- Sys.time()
  print(t_end - t_start)
}
dev.off()

# dim(RESULTS)

write_csv(RESULTS, here("data", out_folder, glue("Pvals_{part}_{trait}.csv")))

t_end_global <- Sys.time()
cat("\nGlobal time\n")
print(t_end_global - t_start_global)
