here::i_am(file.path("data", "01_accessions_phenotypes.R"))

library(here)
library(tidyverse)
library(glue)
library(viridis)


# 1. Parameters --------------------------------------------------------------------------

# HARD-CODED
n_mark <- 10000
out_folder <- "plots"
in_folder <- "top-markers"

# VARIABLE
pars <- read_csv(here("data", "utility", "pars_GWAS.csv"))
# Marker sets to get top markers from
top_markers <- c("all")
pars <- expand.grid(trait = unique(pars$trait),
                    part = unique(pars$part),
                    top_marker = top_markers)


# 2. Plot --------------------------------------------------------------------------------


for (PART in unique(pars$part)) {
  pars_sub <- filter(pars, part == PART)
  
  plot_data <- NULL
  
  for (i in 1:nrow(pars_sub)) {
    trait <- pars_sub[[i, "trait"]]
    part <- pars_sub[[i, "part"]]
    top_marker <- pars_sub[[i, "top_marker"]]
    
    table <- read_csv(here(
      "data", in_folder, glue("Pvals_top-{n_mark}-{top_marker}_{part}_{trait}.csv")
    ))
    
    tmp <- table %>% 
      group_by(marker_set) %>% 
      summarise("{trait}" := n())
    
    if (i == 1) {
      plot_data <- tmp
    } else {
      plot_data <- left_join(plot_data, tmp)
    }
  }
  plot_data
  
  plot_data <- plot_data %>% 
    pivot_longer(
      cols = 2:last_col(),
      values_to = "count",
      names_to = "trait"
    )
  
  plot_data$trait <- factor(plot_data$trait, levels = c("culm.diameter.1st.internode",
                                                        "leaf.senescence",
                                                        "grain.weight",
                                                        "time.to.flowering.from.sowing"))
  
  p1 <- ggplot(plot_data) +
    geom_col(aes(x = trait, y = count, fill = marker_set)) +
    scale_fill_viridis(discrete = TRUE, option = "turbo") +
    xlab("Trait") +
    ylab("No. of markers") +
    labs(fill = "Marker",
         title = glue("{n_mark / 1000}k most associated markers, {part}")) +
    theme_bw()
  
  # PNG
  png(here(
    "data", out_folder, glue("top-{n_mark}_{part}.png")
  ), units = "in", width = 12, height = 7.2, res = 300)
  print(p1)
  dev.off()
}
