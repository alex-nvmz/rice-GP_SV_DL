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

# Load data

DATA <- NULL

for (i in 1:nrow(pars)) {
  trait <- pars[[i, "trait"]]
  part <- pars[[i, "part"]]
  top_marker <- pars[[i, "top_marker"]]
  
  DATA_tmp <- read_csv(here(
    "data", in_folder, glue("Pvals_top-{n_mark}-{top_marker}_{part}_{trait}.csv")
  ))
  
  DATA_tmp <- DATA_tmp %>% 
    group_by(marker_set) %>% 
    summarise("count" = n())
  DATA_tmp$trait <- trait
  DATA_tmp$part <- part
  
  DATA <- bind_rows(DATA, DATA_tmp)
}

DATA

# Plot

plot_data <- DATA

plot_data$trait <- factor(plot_data$trait, levels = c(
  "time.to.flowering.from.sowing",
  "grain.weight",
  "culm.diameter.1st.internode",
  "leaf.senescence"
))

plot_data$marker_set <- factor(
  plot_data$marker_set, levels = c("SNP", "MITE-DTX", "RLX-RIX", "DEL", "DUP", "INV")
)

n_mark
p1 <- ggplot(plot_data) +
  geom_boxplot(aes(
    x = marker_set, y = (count / n_mark), color = marker_set
  ), alpha = 1) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  facet_wrap(vars(trait)) +
  xlab("Marker set") +
  ylab("Proportion") +
  theme_bw()
p1

# PNG
png(here(
  "data", out_folder, glue("top-{n_mark}_across-parts.png")
), units = "in", width = 9, height = 6, res = 300)
print(p1)
dev.off()

# labs(title = glue(
#   "Distributions of marker types of the top {n_mark / 1000}k markers across 10-fold partitions + AroAdm partition"
# )) +


p2 <- ggplot(plot_data) +
  geom_boxplot(aes(
    x = marker_set, y = (count / n_mark)
  ), fill = "coral", alpha = 0.6) +
  facet_wrap(vars(trait)) +
  xlab("Marker set") +
  ylab("Fraction") +
  labs(title = glue(
    "Distributions of marker types of the top {n_mark / 1000}k markers across 10-fold partitions + AroAdm partition"
  )) +
  theme_bw()
p2

png(here(
  "data", out_folder, glue("top-{n_mark}_across-parts_fraction.png")
), units = "in", width = 12, height = 7.2, res = 300)
print(p2)
dev.off()


# ----

# Marker fraction from total available markers

marker_sets <- unique(DATA$marker_set)

total_count <- c()
for (m in marker_sets) {
  count <- read_csv(here(
    "data", "main", glue("geno_{m}.csv")
    ), n_max = 1, col_select = -1) %>%
    ncol()
  
  names(count) <- m
  total_count <- c(total_count, count)
}
total_count

frac_from_total <- apply(plot_data, 1, function(x){
  as.numeric(x["count"]) / total_count[x["marker_set"]]
})

plot_data$frac_from_total <- frac_from_total
plot_data

p3 <- ggplot(plot_data) +
  geom_boxplot(aes(
    x = marker_set, y = frac_from_total * 100, color = marker_set
  )) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  facet_wrap(vars(trait)) +
  xlab("Marker set") +
  ylab("% from total available markers") +
  theme_bw()
p3

png(here(
  "data", out_folder, glue("top-{n_mark}_across-parts_from-total.png")
), units = "in", width = 12, height = 7.2, res = 300)
print(p3)
dev.off()


# labs(title = glue(
#   "Distributions of % from total available markers for the different marker types
#     in the top {n_mark / 1000}k markers, across 10-fold partitions + AroAdm partition"
# )) +