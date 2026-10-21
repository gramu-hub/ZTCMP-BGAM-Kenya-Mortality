# =====================================================================
# STEP 9 (b): BACKTEST - CMP-BGAM
# Probabilistic Forecasting of Mortality in Kenya (ZTCMP-BGAM)
# ---------------------------------------------------------------------
# Run in a fresh session. Refits the CMP-BGAM to the fitting period
# 1950-1993 only, then forecasts the held-out period 1994-2023 and
# scores the forecast rate against the observed rate.
#
# Fit uses only the Year <= 1993 rows; the penalty and the age-specific
# dispersion nu_x are on the parameters, so dropping rows leaves them
# unchanged. The basis is already evaluated at the holdout years in
# design_matrices.rds, so the fitted predictor extrapolates directly:
#   forecast rate  m_hat = exp(eta_hat) = exp(X_int a + X_smooth b)
# (no future exposure needed). Observed rate = D / E in the holdout.
# Forecast relative error = mean(|m_obs - m_hat| / m_obs) over 1994-2023.
#
# nu at age 0 is fixed (index 1) as in the main fit; gamma_nu starts warm
# from the fitted Poisson training run. Reuses the cmp_bgam DLL from 5c.
#
# Input  : design_matrices.rds, poisson_bgam_fit.rds,
#          the cmp_bgam DLL (from 5c)
# Output : backtest_cmp_bgam_1994.rds
# =====================================================================

library(TMB)

design <- readRDS("design_matrices.rds")

CUT   <- 1993                       # last fitting year; holdout is 1994-2023
train <- design$Year <= CUT
hold  <- design$Year >  CUT

K         <- ncol(design$X_smooth)
n_pen     <- design$n_pen
A         <- design$n_ages
S_stacked <- do.call(rbind, design$S_list)

# log factorials for the CMP normaliser (up to max training Y + window).
Yt    <- as.integer(design$Y[train])
max_j <- max(Yt) + 120L
log_fac <- lgamma(seq(0, max_j) + 1)

# Fitting-period data only (free dispersion: add age index + window).
data_train <- list(
  Y            = Yt,
  X_int        = design$X_int[train, , drop = FALSE],
  X_smooth     = design$X_smooth[train, , drop = FALSE],
  log_exposure = design$logExp[train],
  Age_id       = as.integer(design$Age_idx[train]),
  A            = as.integer(A),
  window_lo    = 40L,
  window_hi    = 120L,
  log_fac      = log_fac,
  S_stacked    = S_stacked,
  S_rank       = as.integer(design$S_rank),
  n_pen        = as.integer(n_pen),
  K            = as.integer(K)
)

# Warm start: converged Poisson mean, nu_x = 1 (gamma_nu = 0).
if (!file.exists("poisson_bgam_fit.rds"))
  stop("Run 5a (Poisson) first - it supplies the starting values.")
pf <- readRDS("poisson_bgam_fit.rds")$par_full
par_free <- list(
  alpha       = unname(pf[names(pf) == "alpha"]),
  beta_smooth = unname(pf[names(pf) == "beta_smooth"]),
  log_tau     = unname(pf[names(pf) == "log_tau"]),
  gamma_nu    = rep(0, A)
)

# Fix nu at age 0 (index 1) via map, as in the main fit.
fix_idx <- 1
fix_val <- 0.15
par_free$gamma_nu[fix_idx] <- log(fix_val)
fac <- seq_along(par_free$gamma_nu)
fac[fix_idx] <- NA
map <- list(gamma_nu = factor(fac))

# ---- Load the CMP-BGAM DLL (compiled by 5c) ----
if (!file.exists(dynlib("cmp_bgam")))
  stop("cmp_bgam DLL not found - run 5c first to compile it.")
dyn.load(dynlib("cmp_bgam"))

# ---- Fit to 1950-1993 ----
obj <- MakeADFun(data = data_train, parameters = par_free, DLL = "cmp_bgam",
                 map = map, silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr,
              control = list(iter.max = 2000, eval.max = 2000))
cat("===== BACKTEST: CMP-BGAM =====\n")
cat(" fit convergence code:", opt$convergence, "\n")

# ---- Training-fit convergence guard ----
# nlminb's code 1 is the benign intercept-smooth ridge; the substantive
# checks are a small relative gradient and a positive-definite Hessian.
sdr   <- sdreport(obj)
rel_g <- max(abs(obj$gr(opt$par))) / abs(opt$objective)
cat(" relative gradient:", signif(rel_g, 3), " pdHess:", sdr$pdHess, "\n")
if (!sdr$pdHess || rel_g > 1e-3)
  cat(" *** WARNING: training fit may not have converged cleanly -",
      "treat this backtest with caution.\n")

# ---- Forecast the holdout period ----
# eta_hat = X_int a + X_smooth b (no exposure); rate = exp(eta_hat).
a <- opt$par[names(opt$par) == "alpha"]
b <- opt$par[names(opt$par) == "beta_smooth"]
eta_hat  <- as.numeric(design$X_int %*% a + design$X_smooth %*% b)
rate_hat <- exp(eta_hat)
rate_obs <- design$Y / exp(design$logExp)      # D / E

# ---- Forecast relative error over 1994-2023 (Eq. 3.65) ----
fre <- mean(abs(rate_obs[hold] - rate_hat[hold]) / rate_obs[hold])

cat(" fitting period : 1950 -", CUT, "(", sum(train), "cells )\n")
cat(" holdout period :", CUT + 1, "- 2023 (", sum(hold), "cells )\n")
cat(" forecast relative error :", round(fre, 4), "\n")

out <- list(
  model    = "CMP-BGAM",
  cut      = CUT,
  fre      = fre,
  Age      = design$Age[hold],
  Year     = design$Year[hold],
  rate_obs = rate_obs[hold],
  rate_hat = rate_hat[hold]
)
saveRDS(out, "backtest_cmp_bgam_1994.rds")

dyn.unload(dynlib("cmp_bgam"))
gc()
cat("\nSaved: backtest_cmp_bgam_1994.rds\n")