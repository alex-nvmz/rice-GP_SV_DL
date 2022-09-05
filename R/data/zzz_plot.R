
library(tidyverse)
library(here)
library(cowplot)


DATA <- read_csv(here(
  "data", "utility", "all.final.results.csv"
))

DATA$architecture <- factor(DATA$architecture, levels = c(
  "BayesC", "RKHS", "CNN", "MLP", "MLP_PCA", "MLP_PCA_MULTIPLE"
))
DATA$marker <- factor(DATA$marker, levels = c(
  "ALL", "SNPS", "SIX_Kernels_ONE_INPUT", "SNP_Kernel", "PCA_SNPS", "SIX_PCAS_ONE_INPUT",
  "PCAS_SIX_INPUTS"
))

# ---

pheno <- pheno <- read_csv(here(
  "data", "main", "pheno_original.csv"
))

traits <- pheno[, c(11, 13)]
traits

# Artificial Neural Network (ANN) data
# Traits were not scaled, so approximate results by dividing metrics by variance of the trait

DATA_ANN <- DATA %>% 
  filter(! architecture %in% c("BayesC", "RKHS")) %>% 
  mutate(MSE_GRAIN_WEIGHT = MSE_GRAIN_WEIGHT / (var(traits$grain.weight, na.rm = TRUE)),
         MSE_TIME_TO_FLOWERING = MSE_TIME_TO_FLOWERING /
           (var(traits$time.to.flowering.from.sowing, na.rm = TRUE)))

# Linear models
# Traits were scaled, so not change

DATA_linear <- DATA %>% 
  filter(architecture %in% c("BayesC", "RKHS"))

DATA <- bind_rows(DATA_linear, DATA_ANN)

# ---

# 10-fold-cv plot

plot_data <- DATA
plot_data <- plot_data %>% 
  filter(Partition != "ARO_ADM")

plot_data <- plot_data %>% 
  pivot_longer(cols = c(BIN_CULM_DIAMETER, BIN_LEAF_SENESCENCE),
               values_to = "bin_ce",
               names_to = "bin_traits")

plot_data
range1 <- c(0, range(plot_data$bin_ce)[2])
p1 <- ggplot(plot_data) +
  geom_boxplot(aes(x = architecture, y = bin_ce, color = marker)) +
  facet_grid(cols = vars(bin_traits)) +
  ylim(range1)
p1

plot_data <- DATA
plot_data <- plot_data %>% 
  filter(Partition != "ARO_ADM")
plot_data <- plot_data %>% 
  pivot_longer(cols = c(MSE_GRAIN_WEIGHT, MSE_TIME_TO_FLOWERING),
               values_to = "MSE",
               names_to = "cont_traits")

range2 <- c(0, range(plot_data$MSE)[2])
p2 <- ggplot(plot_data) +
  geom_boxplot(aes(x = architecture, y = MSE, color = marker)) +
  facet_grid(cols = vars(cont_traits)) +
  ylim(range2)
p2

final <- plot_grid(p1, p2, nrow = 2)

png(here(
  "data", "plots", "prediction_alternative_fix.png"
), units = "in", width = 15, height = 9, res = 300)
print(final)
dev.off()

# ----

# AroAdm plot

plot_data <- DATA
plot_data <- plot_data %>% 
  filter(Partition == "ARO_ADM")
plot_data <- plot_data %>% 
  pivot_longer(cols = c(BIN_CULM_DIAMETER, BIN_LEAF_SENESCENCE),
               values_to = "bin_ce",
               names_to = "bin_traits") 

plot_data
p1 <- ggplot(plot_data) +
  geom_point(aes(x = architecture, y = bin_ce, color = marker),
               position = position_dodge(width = 0.5), size = 3) +
  facet_grid(cols = vars(bin_traits)) +
  ylim(range1)
p1

plot_data <- DATA
plot_data <- plot_data %>% 
  filter(Partition == "ARO_ADM")
plot_data <- plot_data %>% 
  pivot_longer(cols = c(MSE_GRAIN_WEIGHT, MSE_TIME_TO_FLOWERING),
             values_to = "MSE",
             names_to = "cont_traits")


p2 <- ggplot(plot_data) +
  geom_point(aes(x = architecture, y = MSE, color = marker),
             position = position_dodge(width = 0.5), size = 3) +
  facet_grid(cols = vars(cont_traits)) +
  ylim(range2)
p2

final <- plot_grid(p1, p2, nrow = 2)

png(here(
  "data", "plots", "prediction_alternative_AroAdm_fix.png"
), units = "in", width = 15, height = 9, res = 300)
print(final)
dev.off()

