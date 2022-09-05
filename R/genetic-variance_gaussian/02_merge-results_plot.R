here::i_am(file.path("genetic-variance_gaussian", "01_h2-estimate_gaussian.R"))

library(here)
library(tidyverse)
library(glue)

# Merge results --------------------------------------------------------------------------
project_folder <- "genetic-variance_gaussian"
out_folder <- "results"


# Features: kernel combinations
K_df <- data.frame(K_1 = "SNP",
                   K_2 = c("RLX-RIX", "MITE-DTX", "DEL", "DUP", "INV"))
K_df


# Merge
RES <- NULL
SNP_df <- NULL
# For every kernel combination
for (i in 1:nrow(K_df)) {
  kernel_names <- K_df[i, ]
  
  # Load results
  in_folder <- glue("out_RKHS_{kernel_names[1]}_{kernel_names[2]}")
  RES_tmp <- read_csv(here(
    project_folder, in_folder, "results.csv"
  ))
  RES_tmp
  
  # Join
  if (i == 1) {
    RES <- RES_tmp[, c(1, 3)]
    SNP_df <- RES_tmp[, 2]
    
  } else {
    RES <- bind_cols(RES, RES_tmp[, 3])
    SNP_df <- bind_cols(SNP_df, RES_tmp[, 2])
  }
}

# Compute mean of SNP heritabilites
RES$SNP_mean <- apply(SNP_df, 1, mean)
RES$SNP_sd <- apply(SNP_df, 1, sd)

# Save 
write_csv(RES, here(
  project_folder, out_folder, "results.csv"
))


# Plots good -----------------------------------------------------------------------------
plot_data <- RES
plot_data <- plot_data %>%
  pivot_longer(cols = RLX.RIX:SNP_mean,
               names_to = "marker_set",
               values_to = "h2") %>% 
  mutate(var_set = ifelse(marker_set == "SNP_mean", "SNP", "SV")) %>% 
  relocate(trait, var_set, marker_set, h2)

plot_data$marker_set <- factor(
  plot_data$marker_set, levels = c("SNP_mean", "MITE.DTX", "RLX.RIX", "DEL", "DUP", "INV")
)

plot_data$trait <- factor(
  plot_data$trait, levels = c(
    "time.to.flowering.from.sowing", "grain.weight", "culm.diameter.1st.internode",
    "leaf.senescence", "culm.strength", "flag.leaf.angle", "grain.length", "grain.width",
    "leaf.length", "salt.injury.at.EC12", "panicle.threshability"
  )
)

library(viridis)
library(ggh4x)
design <- "
AABBCC
DDEEFF
GGHHII
#JJKK#
"

# Column good
p1 <- ggplot(plot_data) +
  geom_col(aes(x = var_set, y = h2, fill = marker_set), position = "dodge") +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  facet_manual(vars(trait), design = design) +
  # facet_wrap(vars(trait)) +
  theme_bw()
p1

png(here(
  project_folder, out_folder, "h2_plot.png"
), units = "in", width = 12, height = 11, res = 300)
print(p1)
dev.off()

# ----

# Points
ggplot(plot_data) +
  geom_point(aes(x = var_set, y = h2, color = marker_set, shape = marker_set),
             size = 3) +
  scale_color_viridis(discrete = TRUE, option = "turbo") +
  facet_wrap(vars(trait)) +
  theme_bw()


# Just 4 traits ---
plot_data <- RES
plot_data <- filter(plot_data, trait %in% c(
  "culm.diameter.1st.internode",
  "leaf.senescence",
  "grain.weight",
  "time.to.flowering.from.sowing"
))

plot_data <- plot_data %>%
  pivot_longer(cols = RLX.RIX:SNP_mean,
               names_to = "marker_set",
               values_to = "h2") %>% 
  mutate(var_set = ifelse(marker_set == "SNP_mean", "SNP", "SV")) %>% 
  relocate(trait, var_set, marker_set, h2)

plot_data$marker_set <- factor(
  plot_data$marker_set, levels = c("SNP_mean", "MITE.DTX", "RLX.RIX", "DEL", "DUP", "INV")
)

plot_data$trait <- factor(
  plot_data$trait, levels = c(
    "time.to.flowering.from.sowing", "grain.weight", "culm.diameter.1st.internode",
    "leaf.senescence"
  )
)

p2 <- ggplot(plot_data) +
  geom_col(aes(x = var_set, y = h2, fill = marker_set), position = "dodge") +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  # facet_manual(vars(trait), design = design) +
  facet_wrap(vars(trait)) +
  theme_bw()
p2

png(here(
  project_folder, out_folder, "h2_plot_4-traits.png"
), units = "in", width = 8, height = 6, res = 300)
print(p2)
dev.off()

# Plot bad -------------------------------------------------------------------------------
plot_data <- RES
plot_data <- plot_data %>% 
  pivot_longer(cols = RLX.RIX:SNP_mean,
               names_to = "marker_set",
               values_to = "h2") %>% 
  relocate(trait, marker_set, h2)

plot_data$marker_set <- factor(
  plot_data$marker_set, levels = c("SNP_mean", "MITE.DTX", "RLX.RIX", "DEL", "DUP", "INV")
)

ggplot(plot_data) +
  geom_col(aes(x = marker_set, y = h2, fill = marker_set)) +
  scale_fill_viridis(discrete = TRUE, option = "turbo") +
  facet_wrap(vars(trait)) +
  theme_bw()

