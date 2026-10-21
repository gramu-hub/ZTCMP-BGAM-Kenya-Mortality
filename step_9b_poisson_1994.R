# =====================================================================
# STEP 9 (b): BACKTEST - POISSON-BGAM
# Probabilistic Forecasting of Mortality in Kenya (ZTCMP-BGAM)
# ---------------------------------------------------------------------
# Run in a fresh session. Refits the Poisson-BGAM to the fitting period
# 1950-1993 only, then forecasts the held-out period 1994-2023 and
# scores the forecast rate against the observed rate.
#
# Fit uses only the Year <= 1993 rows; the penalty is on the coefficients,
# so dropping rows leaves it unchanged. The basis is already evaluated at
# the holdout years in design_matrices.rds, so the fitted predictor
# extrapolates to them directly:
#   forecast rate  m_hat = exp(eta_hat) = exp(X_int a + X_smooth b)
# (no future exposure needed). Observed rate = D / E in the holdout.
# Forecast relative error = mean(|m_obs - m_hat| / m_obs) over 1994-2023.
#
# Reuses the poisson_bgam DLL compiled in 5a (does not recompile).
#
# Input  : design_matrices.rds, the poisson_bgam DLL (from 5a)
# Output : backtest_poisson_bgam_1994.rds
# =====================================================================

library(TMB)

design <- readRDS("design_matrices.rds")

CUT   <- 1993                       # last fitting year; holdout is 1994-2023
train <- design$Year <= CUT
hold  <- design$Year >  CUT

K         <- ncol(design$X_smooth)
n_pen     <- design$n_pen
S_stacked <- do.call(rbind, design$S_list)

# Fitting-period data only.
data_train <- list(
  Y            = as.integer(design$Y[train]),
  X_int        = design$X_int[train, , drop = FALSE],
  X_smooth     = design$X_smooth[train, , drop = FALSE],
  log_exposure = design$logExp[train],
  S_stacked    = S_stacked,
  S_rank       = as.integer(design$S_rank),
  n_pen        = as.integer(n_pen),
  K            = as.integer(K)
)
par_fixed <- list(alpha = -2.6, beta_smooth = rep(0, K), log_tau = rep(0, n_pen))

# ---- Load the Poisson-BGAM DLL (compiled by 5a) ----
if (!file.exists(dynlib("poisson_bgam")))
  stop("poisson_bgam DLL not found - run 5a first to compile it.")
dyn.load(dynlib("poisson_bgam"))

# ---- Fit to 1950-1993 ----
obj <- MakeADFun(data = data_train, parameters = par_fixed,
                 DLL = "poisson_bgam", silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr,
              control = list(iter.max = 2000, eval.max = 2000))
cat("===== BACKTEST: POISSON-BGAM =====\n")
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
  model    = "Poisson-BGAM",
  cut      = CUT,
  fre      = fre,
  Age      = design$Age[hold],
  Year     = design$Year[hold],
  rate_obs = rate_obs[hold],
  rate_hat = rate_hat[hold]
)
saveRDS(out, "backtest_poisson_bgam_1994.rds")

dyn.unload(dynlib("poisson_bgam"))
gc()
cat("\nSaved: backtest_poisson_bgam_1994.rds\n")