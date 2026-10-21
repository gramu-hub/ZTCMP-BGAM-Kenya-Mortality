# =====================================================================
# STEP 11: FORECASTING 2024-2053 (SELECTED MODEL: ZTCMP-BGAM) - Objective 4
# Probabilistic Forecasting of Mortality in Kenya (ZTCMP-BGAM)
# ---------------------------------------------------------------------
# Run in a fresh session. Reads the saved ZTCMP fit and rebuilds the
# spline basis at future years. Refits nothing.
#
# Extrapolation limitation. Forecasting evaluates the time smooth f2(t)
# and the time margin of f3(x,t) beyond the observed range (1950-2023).
# Penalised spline bases are not designed to extrapolate: beyond the data
# the basis continues its last local trend with no guarantee it stays
# demographically sensible, and the uncertainty grows rapidly. The
# forecasts are therefore reported with their widening intervals and
# interpreted with increasing caution over the horizon - a genuine
# limitation relative to, e.g., a Lee-Carter model with an explicit
# random-walk-with-drift. It is reported, not hidden.
#
# Age 100 is the open interval 100+; the life-expectancy calculation
# closes the terminal age accordingly.
#
# Uncertainty. The mean-structure parameters (alpha, beta_smooth) are
# drawn from the Laplace-approximate posterior covariance (sdreport), the
# linear predictor is rebuilt for each draw, and the predictive quantiles
# of the rate are read off. The rate is exp(eta) and does not contain nu,
# so nu is not itself sampled; the age-specific dispersion instead shapes
# this covariance - overdispersed ages are less informative about the
# mean, so their parameters are less precisely estimated - so the fan
# widths differ by age.
#
# Input  : ztcmp_bgam_fit.rds (with $sdr), design_matrices.rds,
#          mortality_long.rds
# Output : fig11_fan_charts.png, fig12_forecast_surface.png,
#          fig13_life_expectancy.png, table11_forecast_summary.xlsx,
#          table12_life_expectancy.xlsx, appendix_full_forecast_2024_2053.xlsx
# =====================================================================

library(mgcv)
library(openxlsx)
library(ggplot2)
library(MASS)                    # mvrnorm
set.seed(2024)

design <- readRDS("design_matrices.rds")
ztcmp  <- readRDS("ztcmp_bgam_fit.rds")
dat    <- readRDS("mortality_long.rds")
if (is.null(ztcmp$sdr))
  stop("ztcmp fit has no $sdr - rerun 5d with keep_sdr = TRUE.")

# ---- 1. Rebuild the basis at future years ----
# Fit a throwaway gam with the same formula and data as Step 3 to get a
# predict() method; its bases are identical, so predict(type = "lpmatrix")
# gives the design matrix at any Age x Year, including future years.
# Coefficients are overwritten with the fitted ZTCMP betas; nothing is
# refitted for the forecast.
shell <- gam(Deaths ~ s(Age,  k = 20, bs = "ps") +
               s(Year, k = 20, bs = "ps") +
               ti(Age, Year, k = c(10, 10), bs = "cr"),
             data = dat, family = poisson, offset = log(dat$Exposure))

# ---- 2. Future prediction grid ----
# Future exposure is unknown, so forecast RATES: set the offset aside and
# report mu/E = exp(eta). This is what life expectancy needs and does not
# require projecting exposure.
ages_all   <- 0:100
yrs_future <- 2024:2053
grid <- expand.grid(Age = ages_all, Year = yrs_future)

Xf <- predict(shell, newdata = grid, type = "lpmatrix")   # basis at future cells

# ztcmp betas in the same [intercept, smooths] order as the basis columns.
beta_hat <- c(ztcmp$par_full[names(ztcmp$par_full) == "alpha"],
              ztcmp$par_full[names(ztcmp$par_full) == "beta_smooth"])
stopifnot(length(beta_hat) == ncol(Xf))

eta_f      <- as.vector(Xf %*% beta_hat)     # log rate at future cells
log_rate_f <- eta_f                          # offset 0 => log(mu/E)

# ---- 3. Posterior draws for uncertainty ----
# Draw the mean-structure parameters from the Laplace posterior covariance,
# rebuild the future log-rate for each draw, and take quantiles.
sdr         <- ztcmp$sdr
Vfixed      <- sdr$cov.fixed
fixed_names <- rownames(Vfixed)
idx_beta    <- which(fixed_names %in% c("alpha", "beta_smooth"))
mu_beta     <- sdr$par.fixed[idx_beta]
V_beta      <- Vfixed[idx_beta, idx_beta]

B     <- 500
draws <- MASS::mvrnorm(B, mu = mu_beta, Sigma = V_beta)   # B x length(beta)

eta_draws <- Xf %*% t(draws)                               # (n_future) x B
qs <- t(apply(eta_draws, 1, quantile, probs = c(0.025, 0.5, 0.975)))
grid$lo  <- exp(qs[, 1]); grid$med <- exp(qs[, 2]); grid$hi <- exp(qs[, 3])
grid$rate_point <- exp(log_rate_f)

# ---- 4. Figure 11: fan charts for representative ages ----
sel_ages <- c(0, 5, 30, 60, 85)
dffan <- subset(grid, Age %in% sel_ages)
dffan$Age <- factor(paste("Age", dffan$Age), levels = paste("Age", sel_ages))
library(ggh4x)

scales_list <- list(
  "Age 0"  = scale_y_continuous(breaks = seq(0.002, 0.010, by = 0.002),
                                limits = c(0.002, 0.010)),
  "Age 5" = scale_y_continuous(breaks = seq(0.0002, 0.0008, by = 0.0002),
                               limits = c(0.0002, 0.0008),
                               labels = scales::label_number(accuracy = 0.0001)),
  "Age 30" = scale_y_continuous(breaks = seq(0.001, 0.003, by = 0.001),
                                limits = c(0.001, 0.003)),
  "Age 60" = scale_y_continuous(breaks = seq(0.006, 0.018, by = 0.006),
                                limits = c(0.006, 0.018)),
  "Age 85" = scale_y_continuous(breaks = seq(0.04, 0.12, by = 0.04),
                                limits = c(0.04, 0.12))
)

p11 <- ggplot(dffan, aes(Year, rate_point)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.25) +
  geom_line(colour = "steelblue") +
  facet_wrap(~ Age, scales = "free_y") +
  facetted_pos_scales(y = scales_list) +
  labs(x = "Year", y = "mortality rate") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave("fig11_fan_charts.png", p11, width = 9, height = 6, dpi = 300)

# ---- 5. Figure 12: forecast surface ----
p12 <- ggplot(grid, aes(Year, Age, fill = log(rate_point))) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma") +
  scale_y_continuous(breaks = seq(0, 100, by = 5),
                     labels = c(as.character(seq(0, 95, by = 5)), "100+")) +
  scale_x_continuous(breaks = seq(2025, 2050, by = 5)) +
  labs(x = "Year", y = "Age (single years, 100+ open interval)", fill = "log(rate)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
ggsave("fig12_forecast_surface.png", p12, width = 12, height = 7, dpi = 300)

# ---- 6. Life expectancy at birth by year ----
# Period life table from forecast rates m(x). n = 1 single-year ages;
# a(x) = 0.5 except the infant correction a(0) = 0.3; terminal age 100 is
# the open interval 100+, closed with e(100) = 1/m(100).
life_exp_at_birth <- function(m) {          # m = rates, ages 0..100
  ax <- rep(0.5, length(m)); ax[1] <- 0.3   # infant correction
  qx <- m / (1 + (1 - ax) * m)
  qx[length(qx)] <- 1                       # open interval closes
  px <- 1 - qx
  lx <- c(1, cumprod(px)[-length(px)])
  dx <- lx * qx
  Lx <- lx - (1 - ax) * dx
  Lx[length(Lx)] <- ifelse(m[length(m)] > 0, lx[length(lx)] / m[length(m)], 0)
  sum(Lx) / lx[1]
}

le_point <- sapply(yrs_future, function(y) {
  r <- grid$rate_point[grid$Year == y][order(ages_all)]
  life_exp_at_birth(r)
})

# uncertainty: life expectancy for each posterior draw
le_draws <- matrix(NA, length(yrs_future), B)
for (b in seq_len(B)) {
  rate_b <- exp(matrix(eta_draws[, b], nrow = length(ages_all)))   # ages x years
  for (j in seq_along(yrs_future)) {
    le_draws[j, b] <- life_exp_at_birth(rate_b[, j])
  }
}
le_lo <- apply(le_draws, 1, quantile, 0.025)
le_hi <- apply(le_draws, 1, quantile, 0.975)
dfle  <- data.frame(Year = yrs_future, LE = le_point, lo = le_lo, hi = le_hi)

# ---- 7. Figure 13: life expectancy trajectory ----
p13 <- ggplot(dfle, aes(Year, LE)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "darkgreen", alpha = 0.2) +
  geom_line(colour = "darkgreen", linewidth = 0.8) +
  scale_x_continuous(breaks = seq(2025, 2050, by = 5)) +
  scale_y_continuous(breaks = seq(60, 90, by = 5)) +
  labs(x = "Year", y = "life expectancy at birth (years)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold"))
ggsave("fig13_life_expectancy.png", p13, width = 10, height = 6, dpi = 300)

# ---- 8. Table 11: forecast summary at 2030 / 2040 / 2053 ----
sel_yrs <- c(2030, 2040, 2053)
t11 <- subset(grid, Year %in% sel_yrs & Age %in% sel_ages,
              select = c(Age, Year, rate_point, lo, hi))
t11[c("rate_point", "lo", "hi")] <- lapply(t11[c("rate_point", "lo", "hi")],
                                           function(z) signif(z, 3))
t11 <- t11[order(t11$Age, t11$Year), ]
names(t11) <- c("Age", "Year", "Forecast_rate", "CI_lower", "CI_upper")
cat("=============== TABLE 11: FORECAST RATES (sel ages/years) ===============\n")
print(t11, row.names = FALSE)
write.xlsx(t11, "table11_forecast_summary.xlsx")

# Full forecast (all ages x all years) as the appendix behind Table 11.
full_forecast <- grid[order(grid$Age, grid$Year),
                      c("Age", "Year", "rate_point", "lo", "hi")]
full_forecast[c("rate_point", "lo", "hi")] <-
  lapply(full_forecast[c("rate_point", "lo", "hi")], function(z) signif(z, 4))
names(full_forecast) <- c("Age", "Year", "Forecast_rate", "CI_lower", "CI_upper")
full_forecast$Age <- ifelse(full_forecast$Age == 100, "100+",
                            as.character(full_forecast$Age))   # terminal open interval
write.xlsx(full_forecast, "appendix_full_forecast_2024_2053.xlsx")
cat(" (full", nrow(full_forecast), "-row forecast saved to appendix)\n")

# ---- 9. Table 12: life expectancy by year ----
t12 <- data.frame(Year = yrs_future,
                  Life_expectancy = round(le_point, 2),
                  CI_lower = round(le_lo, 2),
                  CI_upper = round(le_hi, 2))
cat("\n=============== TABLE 12: FORECAST LIFE EXPECTANCY ===============\n")
print(t12[t12$Year %in% c(2024, 2030, 2040, 2053), ], row.names = FALSE)
write.xlsx(t12, "table12_life_expectancy.xlsx")

cat("\nSaved: fig11-13, table11_forecast_summary.xlsx,",
    "table12_life_expectancy.xlsx, appendix_full_forecast_2024_2053.xlsx\n")
cat("\nReminder: forecasts rely on spline extrapolation beyond 2023 -",
    "report the widening intervals as a limitation.\n")