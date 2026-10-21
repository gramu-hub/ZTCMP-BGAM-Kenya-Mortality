# =====================================================================
# STEP 2: EXPLORATORY DATA ANALYSIS
# Probabilistic Forecasting of Mortality in Kenya (ZTCMP-BGAM)
# ---------------------------------------------------------------------
# Two figures describe the data before modelling: the mortality surface
# over the age-time plane, and the distribution of the death counts.
#
# Input  : mortality_long.rds
# Output : fig1_mortality_surface.png
#          fig2_count_distribution.png
# =====================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

dat <- readRDS("mortality_long.rds")

thm <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

# =====================================================================
# FIGURE 1: Mortality surface (age x time)
# ---------------------------------------------------------------------
# Log mortality rate over the full age-time plane: the surface the
# model estimates. 
# =====================================================================
p1 <- ggplot(dat, aes(x = Year, y = Age, fill = log(Rate))) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma") +
  scale_y_continuous(breaks = seq(0, 100, by = 5),
                     labels = c(as.character(seq(0, 95, by = 5)), "100+")) +
  scale_x_continuous(breaks = seq(1950, 2020, by = 5)) +
  labs(x = "Year",
       y = "Age (single years, 100+ open interval)",
       fill = "log(rate)") +
  thm

ggsave("fig1_mortality_surface.png", p1, width = 10, height = 6, dpi = 300)

# =====================================================================
# FIGURE 2: Distribution of death counts
# ---------------------------------------------------------------------
# The death counts span 2 to 78,730, a range of roughly 39,000-fold.
# This dynamic range is why the CMP normalising constant is summed over
# a window rather than from zero. A log x-axis is used, since a linear
# scale would compress almost all cells into the first bin.
# =====================================================================
p2 <- ggplot(dat, aes(x = Deaths)) +
  geom_histogram(bins = 50, fill = "steelblue", colour = "white") +
  scale_x_log10(labels = comma) +
  labs(x = "Death count (log scale)",
       y = "Number of age-year cells") +
  thm

ggsave("fig2_count_distribution.png", p2, width = 8, height = 5, dpi = 300)

cat("Done. Two EDA figures written:\n",
    "  fig1_mortality_surface.png\n",
    "  fig2_count_distribution.png\n")