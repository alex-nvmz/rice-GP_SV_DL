here::i_am(file.path("data", "01_accessions_phenotypes.R"))

library(here)
library(tidyverse)
library(glue)
library(cowplot)

# ----------------------------------------------------------------------------------------

# Prepare data
out_folder <- "plots"
pheno <- read_csv(here("data", "main", "pheno_original.csv"))

traits <- pheno[, c(2, 4, 10, 11, 13)]

traits$culm.diameter.1st.internode <- factor(traits$culm.diameter.1st.internode)
traits$leaf.senescence <- factor(traits$leaf.senescence)

# Make plots

make_plot_binary <- function(var_name) {
  p1 <- ggplot(traits) +
    geom_histogram(
      data = filter(traits, ! is.na(!!sym(var_name))),
      aes(x = !!sym(var_name)), stat = "count", width = 0.5
    ) +
    theme_bw() +
    xlab("") +
    ylab("Count") +
    labs(title = var_name)
  
  p2 <- ggplot(traits) +
    geom_bar(
      data = filter(traits, ! is.na(!!sym(var_name))),
      aes(x = !!sym(var_name), y = ..prop.., group = 1, fill = Population),
    ) +
    facet_grid(cols = vars(Population)) +
    theme_bw() +
    labs(fill = "Group") +
    xlab("") +
    ylab("Proportion") +
    labs(title = var_name)
  
  return(list(p1, p2))
}

culm <- make_plot_binary(var_name = "culm.diameter.1st.internode")

leaf <- make_plot_binary(var_name = "leaf.senescence")

make_plot_cont <- function(var_name, units, binwidth) {
  p1 <- ggplot(traits) +
    geom_histogram(aes(x = !!sym(var_name)), binwidth = binwidth) +
    theme_bw() +
    xlab(units) +
    ylab("Count") +
    labs(title = var_name)
  
  p2 <- ggplot(traits) +
    geom_boxplot(aes(x = Population, y = !!sym(var_name), fill = Population),
                 alpha = 0.7) +
    theme_bw() +
    xlab("") +
    ylab(units) +
    labs(fill = "Group",
         title = var_name)
  
  return(list(p1, p2))
}


time <- make_plot_cont(var_name = "time.to.flowering.from.sowing",
                       units = "log( Days )",
                       binwidth = 0.06)

weight <- make_plot_cont(var_name = "grain.weight",
                       units = "Grams",
                       binwidth = 0.2)

# Join plots

plot_list <- time %>% 
  append(weight) %>% 
  append(culm) %>% 
  append(leaf)

grid <- plot_grid(plotlist = plot_list, ncol = 2, labels = "AUTO")
grid

png(here(
  "data", out_folder, "trait_distributions.png"
), units = "in", width = 10, height = 11, res = 300)
print(grid)
dev.off()


# ----------------------------------------------------------------------------------------

# Trait plots, correlations

# Load and prepare data

pheno <- read_csv(here("data", "main", "pheno_original.csv"))
traits <- pheno[, c(2, 4, 10, 11, 13)]
traits <- rename(traits,
       "culm.diameter" = culm.diameter.1st.internode,
       "time.to.flowering" = time.to.flowering.from.sowing)

# PCA

pca_mat <- traits[, -1] %>% 
  na.omit() %>% 
  as.matrix() %>% 
  scale(center = TRUE, scale = TRUE)

res_pca <- prcomp(pca_mat, center = FALSE, scale. = FALSE)

library(factoextra)
library(viridis)

pca_plot <- fviz_pca_var(res_pca,
             col.var = "contrib",
             gradient.cols = viridis(n = 3, option = "viridis"),
             repel = TRUE) +
  theme_bw() +
  labs(color = "Contribution (%)",
       title = element_blank()) +
  xlab(glue("PC1 ({round(((res_pca$sdev^2)[1] / sum(res_pca$sdev^2)) * 100, 2)}%)")) +
  ylab(glue("PC2 ({round(((res_pca$sdev^2)[2] / sum(res_pca$sdev^2)) * 100, 2)}%)"))
pca_plot

png(here(
  "data", out_folder, "trait_PCA.png"
), units = "in", width = 8, height = 8, res = 300)
print(pca_plot)
dev.off()

# Observations plotted on PCs of traits, not interesting

library(ggfortify)
autoplot(res_pca, data = na.omit(traits), colour = "Population", alpha = 0.7) +
  theme_bw()

# Cor plot

library(corrplot)

data <- na.omit(traits[, -1])
cor_mat <- cor(data)

corrplot(cor_mat, method = "color",
   type = "upper", order = "hclust",
   addCoef.col = "black",
   tl.srt = 45, tl.col = "black",
   number.cex = 1.2,
   diag = FALSE)

png(here(
  "data", out_folder, "trait_cor.png"
), units = "in", width = 7, height = 7, res = 300)

corrplot(cor_mat, method = "color",
         type = "upper", order = "hclust",
         addCoef.col = "black",
         tl.srt = 45, tl.col = "black",
         number.cex = 1.2,
         diag = FALSE)

dev.off()

plot(traits$grain.weight, traits$leaf.senescence)
plot(traits$grain.weight, traits$time.to.flowering)
plot(traits$culm.diameter, traits$time.to.flowering)
plot(traits$culm.diameter, traits$time.to.flowering)
plot(traits$culm.diameter, traits$leaf.senescence)
