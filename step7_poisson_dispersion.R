# =====================================================================
# STEP 7: DISPERSION DIAGNOSIS FROM THE FITTED POISSON (POISSON ONLY)
# Probabilistic Forecasting of Mortality in Kenya (ZTCMP-BGAM)
# ---------------------------------------------------------------------
# Run in a fresh session. Reads the saved Poisson fit; refits nothing.
#
# The Poisson assumes equidispersion, so it cannot model departures from
# it, but its residuals reveal them: the mean squared Pearson residual
# within each age is ~1 under equidispersion, >1 over, <1 under. This is
# a model-agnostic diagnosis, compared later against the ZTCMP nu. (The
# ZTP diagnosis is identical, truncation being inert, so it is not
# duplicated.)
#
# Classification threshold (same scale as 1/nu in the ZTCMP step):
#   over  if dispersion ratio > 1.053  (nu < 0.95)
#   under if dispersion ratio < 0.962  (nu > 1.04)
#   equi  otherwise
#
# Input  : poisson_bgam_fit.rds, design_matrices.rds
# Output : fig5_dispersion_ratio_by_age.png,
#          table5_poisson_dispersion_counts.xlsx,
#          appendix_poisson_dispersion_by_age.xlsx
# =====================================================================

library(openxlsx)
library(ggplot2)

design <- readRDS("design_matrices.rds")
pois   <- readRDS("poisson_bgam_fit.rds")

Y    <- design$Y
Age  <- design$Age
mu   <- exp(pois$log_mu)

# Pearson residual: (Y - mu) / sqrt(mu). Under a correct Poisson its
# variance is ~1; departures from 1 signal mis-fit dispersion.
pearson <- (Y - mu) / sqrt(mu)

OVER_T  <- 1.053   # ratio above this = overdispersed
UNDER_T <- 0.962   # ratio below this = underdispersed

# ---- Per-age dispersion ratio ----
# Mean squared Pearson residual within each age = empirical VMR.
disp_ratio <- tapply(pearson^2, Age, mean)
ages    <- as.integer(names(disp_ratio))
verdict <- ifelse(disp_ratio > OVER_T, "over",
                  ifelse(disp_ratio < UNDER_T, "under", "equi"))

# ---- Figure 5: per-age ratio (log scale, line at 1) ----
dfrat <- data.frame(Age = ages, ratio = as.numeric(disp_ratio),
                    verdict = verdict)
p5 <- ggplot(dfrat, aes(Age, ratio, colour = verdict)) +
  geom_point(size = 1.6) +
  geom_hline(yintercept = 1,       colour = "grey30", linetype = "dashed") +
  geom_hline(yintercept = OVER_T,  colour = "grey70", linetype = "dotted") +
  geom_hline(yintercept = UNDER_T, colour = "grey70", linetype = "dotted") +
  scale_y_log10(breaks = 10^seq(-1, 3, by = 0.5),
                labels = as.character(seq(-1, 3, by = 0.5))) +
  scale_x_continuous(breaks = seq(0, 100, by = 5),
                     labels = c(as.character(seq(0, 95, by = 5)), "100+")) +
  scale_colour_manual(values = c(over = "firebrick", equi = "grey50",
                                 under = "steelblue")) +
  labs(x = "Age (single years, 100+ open interval)",
       y = "dispersion ratio (log10 scale)", colour = NULL) +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave("fig5_dispersion_ratio_by_age.png", p5, width = 12, height = 6, dpi = 300)

# ---- Table 5: classification counts ----
table5 <- data.frame(
  Dispersion = c("Overdispersed (ratio > 1.053)",
                 "Equidispersed (0.962 - 1.053)",
                 "Underdispersed (ratio < 0.962)",
                 "TOTAL"),
  Count = c(sum(verdict == "over"),
            sum(verdict == "equi"),
            sum(verdict == "under"),
            length(verdict))
)
cat("============= TABLE 5: POISSON DISPERSION DIAGNOSIS =============\n")
print(table5, row.names = FALSE)
cat("\n ratio range:", round(range(disp_ratio), 3), "\n\n")
write.xlsx(table5, "table5_poisson_dispersion_counts.xlsx")

# ---- Appendix: full per-age ratios ----
# The terminal age (100) is the open interval 100+.
age_label <- as.character(ages)
age_label[ages == max(ages)] <- paste0(max(ages), "+")
appendix <- data.frame(
  Age = age_label,
  Dispersion_ratio = round(as.numeric(disp_ratio), 4),
  Classification = verdict
)
write.xlsx(appendix, "appendix_poisson_dispersion_by_age.xlsx")

cat("Saved: fig5_dispersion_ratio_by_age.png,",
    "table5_poisson_dispersion_counts.xlsx,",
    "appendix_poisson_dispersion_by_age.xlsx\n")