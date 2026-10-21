# =====================================================================
# STEP 10: UNCERTAINTY QUANTIFICATION (SELECTED MODEL: ZTCMP-BGAM)
# Probabilistic Forecasting of Mortality in Kenya (ZTCMP-BGAM)
# ---------------------------------------------------------------------
# Run in a fresh session. Reads the saved ZTCMP fit (which stored its
# sdreport via keep_sdr = TRUE in 5d) and the Poisson fit. Refits nothing.
#
# sdreport() forms the Laplace approximation to the posterior: the mode
# (the MAP estimates) plus the posterior covariance from the Hessian at
# the optimum. This makes the "Bayesian" description concrete - priors
# combined with a Laplace-approximate posterior, not MCMC - and supplies
# the uncertainty the forecast needs. The positive-definite Hessian
# confirmed at the estimation-validation stage is what makes these
# standard errors valid.
#
# nu at age 0 is fixed, so it has no standard error and is reported as
# fixed throughout.
#
# Input  : ztcmp_bgam_fit.rds (with $sdr), poisson_bgam_fit.rds,
#          design_matrices.rds
# Output : table10_parameter_estimates.xlsx,
#          appendix_nu_credible_intervals.xlsx,
#          fig9_interval_widths.png, fig10_nu_credible_bands.png
# =====================================================================

library(openxlsx)
library(ggplot2)

design <- readRDS("design_matrices.rds")
ztcmp  <- readRDS("ztcmp_bgam_fit.rds")
pois   <- readRDS("poisson_bgam_fit.rds")

if (is.null(ztcmp$sdr))
  stop("ztcmp_bgam_fit.rds has no $sdr - rerun 5d with keep_sdr = TRUE.")
sdr  <- ztcmp$sdr
ages <- 0:100

# ---- nu with standard errors ----
# nu = exp(gamma_nu). The interval is built on the gamma_nu (log) scale
# and exponentiated, so it is guaranteed positive without a floor. The
# ADREPORTed nu gives the point estimate; gamma_nu and its SE give the
# interval. The fixed age (0) has no SE.
ad      <- summary(sdr, "report")           # ADREPORTed nu (point estimate)
nu_rows <- grep("^nu$", rownames(ad))
nu_est  <- ad[nu_rows, "Estimate"]

fe        <- summary(sdr, "fixed")
g_rows    <- grep("^gamma_nu$", rownames(fe))
gamma_hat <- fe[g_rows, "Estimate"]
gamma_se  <- fe[g_rows, "Std. Error"]

# gamma_nu covers only the estimated ages; age 0 is held fixed (mapped
# out), so it has no row. Align the SEs and interval back to all 101 ages.
nu_se  <- rep(NA_real_, length(ages))
nu_lo  <- rep(NA_real_, length(ages))
nu_hi  <- rep(NA_real_, length(ages))
fixed0 <- ages == 0
est_ages <- !fixed0
nu_se[est_ages] <- nu_est[est_ages] * gamma_se          # delta-method SE of nu
nu_lo[est_ages] <- exp(gamma_hat - 1.96 * gamma_se)     # positive by construction
nu_hi[est_ages] <- exp(gamma_hat + 1.96 * gamma_se)

# ---- Figure 10: nu by age with credible bands ----
dfnu <- data.frame(Age = ages, nu = nu_est, lo = nu_lo, hi = nu_hi, fixed = fixed0)
p10 <- ggplot(dfnu, aes(Age, nu)) +
  geom_ribbon(data = subset(dfnu, !fixed), aes(ymin = lo, ymax = hi),
              fill = "steelblue", alpha = 0.25) +
  geom_hline(yintercept = 1, colour = "grey30", linetype = "dashed") +
  geom_line(colour = "steelblue") +
  geom_point(data = subset(dfnu, fixed), aes(Age, nu),
             colour = "firebrick", shape = 17, size = 2.5) +
  scale_x_continuous(breaks = seq(0, 100, by = 5),
                     labels = c(as.character(seq(0, 95, by = 5)), "100+")) +
  scale_y_continuous(breaks = seq(0, 8.5, by = 0.5)) +
  labs(x = "Age (single years, 100+ open interval)", y = expression(nu[x])) +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave("fig10_nu_credible_bands.png", p10, width = 12, height = 6, dpi = 300)

# ---- Predictive interval widths: ZTCMP vs Poisson ----
# Plug-in predictive width at each fitted cell, using each model's own
# variance: Poisson Var = mu (half-width ~ 1.96*sqrt(mu)); ZTCMP
# Var ~ mu/nu (half-width ~ 1.96*sqrt(mu/nu)). Shows how the age-specific
# dispersion reshapes the uncertainty - wider where overdispersed (nu<1),
# narrower where underdispersed (nu>1). This is the fitted-cell width, not
# the posterior-predictive forecast interval (that is the forecast step).
mu_p      <- exp(pois$log_mu)
mu_z      <- exp(ztcmp$log_mu)
nu_by_obs <- nu_est[design$Age_idx]
width_pois  <- 2 * 1.96 * sqrt(mu_p)
width_ztcmp <- 2 * 1.96 * sqrt(mu_z / nu_by_obs)

relw_p <- tapply(width_pois  / mu_p, design$Age, mean)
relw_z <- tapply(width_ztcmp / mu_z, design$Age, mean)
dfw <- rbind(
  data.frame(Age = ages, relwidth = as.numeric(relw_p), Model = "Poisson-BGAM"),
  data.frame(Age = ages, relwidth = as.numeric(relw_z), Model = "ZTCMP-BGAM")
)
p9 <- ggplot(dfw, aes(Age, relwidth, colour = Model)) +
  geom_line() + geom_point(size = 1) +
  scale_x_continuous(breaks = seq(0, 100, by = 5),
                     labels = c(as.character(seq(0, 95, by = 5)), "100+")) +
  scale_y_continuous(breaks = seq(0, 2, by = 0.1)) +
  scale_colour_manual(values = c("Poisson-BGAM" = "grey45", "ZTCMP-BGAM" = "purple4")) +
  labs(x = "Age (single years, 100+ open interval)",
       y = "relative 95% interval width", colour = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")
ggsave("fig9_interval_widths.png", p9, width = 12, height = 6, dpi = 300)

# ---- Table 10: key parameter estimates with SEs ----
fe <- summary(sdr, "fixed")
keep_par <- rownames(fe) %in% c("alpha", "log_tau")
tab_fe <- data.frame(
  Parameter = rownames(fe)[keep_par],
  Estimate  = round(fe[keep_par, "Estimate"], 5),
  Std_Error = round(fe[keep_par, "Std. Error"], 3)
)
is_lt <- tab_fe$Parameter == "log_tau"
tab_fe$Parameter[is_lt] <- paste0("log_tau[", seq_len(sum(is_lt)), "] (",
                                  design$S_label, ")")
tab_fe$tau <- ifelse(grepl("log_tau", tab_fe$Parameter),
                     round(exp(tab_fe$Estimate), 3), NA)
cat("=============== TABLE 10: KEY PARAMETER ESTIMATES ===============\n")
print(tab_fe, row.names = FALSE)
cat("\n (full per-age nu with intervals in the appendix)\n\n")
write.xlsx(tab_fe, "table10_parameter_estimates.xlsx")

# ---- Appendix: nu with credible intervals, all ages ----
age_lab <- as.character(ages)
age_lab[length(age_lab)] <- "100+"          # terminal age is the open interval
appendix <- data.frame(
  Age       = age_lab,
  nu_x      = round(nu_est, 4),
  Std_Error = round(nu_se, 4),
  CI_lower  = round(nu_lo, 4),
  CI_upper  = round(nu_hi, 4),
  Note      = ifelse(fixed0, "fixed", ""),
  check.names = FALSE
)
write.xlsx(appendix, "appendix_nu_credible_intervals.xlsx")

cat("Saved: table10_parameter_estimates.xlsx,",
    "appendix_nu_credible_intervals.xlsx,",
    "fig9_interval_widths.png, fig10_nu_credible_bands.png\n")