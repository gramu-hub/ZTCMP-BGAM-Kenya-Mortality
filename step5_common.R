# =====================================================================
# STEP 5 - SHARED SETUP
# Probabilistic Forecasting of Mortality in Kenya (ZTCMP-BGAM)
# ---------------------------------------------------------------------
# Sourced by each of the four fitting scripts so the design, priors and
# diagnostics are identical across models. One model per session (TMB is
# sensitive to several compiled DLLs attached at once).
#
# beta is optimised jointly with the smoothing parameters tau under the
# P-spline penalty; joint estimation keeps the mean near the observed
# count, so the windowed CMP normalising constant stays valid.
# =====================================================================

library(TMB)

design <- readRDS("design_matrices.rds")
K     <- ncol(design$X_smooth)   # 119 smooth columns: f1 (19), f2 (19), f3 (81)
n_pen <- design$n_pen            # 4 smoothing parameters

# Penalty matrices stacked into one (K*n_pen) x K matrix; the .cpp pulls
# out block k with S_stacked.block(k*K, 0, K, K).
S_stacked <- do.call(rbind, design$S_list)
stopifnot(nrow(S_stacked) == K * n_pen, ncol(S_stacked) == K)

# Data for the fixed-dispersion models (Poisson, ZTP)
data_fixed <- list(
  Y            = as.integer(design$Y),
  X_int        = design$X_int,
  X_smooth     = design$X_smooth,
  log_exposure = design$logExp,
  S_stacked    = S_stacked,
  S_rank       = as.integer(design$S_rank),   # 18, 18, 72, 72
  n_pen        = as.integer(n_pen),
  K            = as.integer(K)
)

# Free-dispersion models (CMP, ZTCMP): add the age index for nu_x and
# the CMP normaliser window.
data_free <- c(data_fixed, list(
  Age_id    = as.integer(design$Age_idx),   # 1..101
  A         = as.integer(design$n_ages),    # 101
  window_lo = 40L,                          # CMP normaliser window below the count
  window_hi = 120L                          # CMP normaliser window above the count
))

# Starting values for the fixed-dispersion models
par_fixed <- list(alpha       = -2.6,
                  beta_smooth = rep(0, K),
                  log_tau     = rep(0, n_pen))

# ---------------------------------------------------------------------
# fit_model(): one routine for all four models, so the optimiser
# settings, diagnostics and reporting are identical between them.
# ---------------------------------------------------------------------
fit_model <- function(dll, params, data, label, save_as, keep_sdr = FALSE,
                      fix_nu_ages = NULL, fix_nu_value = 0.15) {
  cat("=====", label, "=====\n"); flush.console()
  t0 <- Sys.time()
  
  # Fix nu at the given ages by holding gamma_nu constant via TMB's map
  # argument (those ages leave the optimisation); the rest are estimated.
  # fix_nu_ages: positions into gamma_nu (age + 1). NULL estimates all.
  map <- NULL
  if (!is.null(fix_nu_ages)) {
    g0 <- params$gamma_nu
    g0[fix_nu_ages] <- log(fix_nu_value)     # set the fixed ages
    params$gamma_nu <- g0
    fac <- seq_along(g0)
    fac[fix_nu_ages] <- NA                    # NA => held constant by TMB
    map <- list(gamma_nu = factor(fac))
    cat("  fixing nu at ages (index):", paste(fix_nu_ages, collapse = ", "),
        "to", fix_nu_value, "\n")
  }
  
  obj <- MakeADFun(data = data, parameters = params, DLL = dll,
                   map = map, silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr,
                control = list(iter.max = 2000, eval.max = 2000))
  
  # ---- Diagnostics ----
  # TMB optimises rather than samples, so convergence is judged by the
  # gradient (relative to the objective, and by parameter block - gamma_nu
  # is the block that matters for the CMP models) and a positive-definite
  # Hessian. Two checks override the convergence code: loglik must be
  # negative, and the relative error must be plausible.
  g <- obj$gr(opt$par)
  mgc <- max(abs(g))
  rel_g <- mgc / abs(opt$objective)
  g_by_block <- tapply(abs(g), names(opt$par), max)
  sdr <- sdreport(obj)
  par_full <- obj$env$last.par.best
  rep <- obj$report(par_full)
  mu <- exp(rep$log_mu)
  rel_err <- mean(abs((data$Y - mu) / data$Y))
  corr <- cor(data$Y, mu)
  el <- round(difftime(Sys.time(), t0, units = "mins"), 2)
  
  cat("  convergence code  :", opt$convergence, " (0 = success)\n")
  cat("  message           :", opt$message, "\n")
  cat("  max |gradient|    :", format(mgc, digits = 3), "\n")
  cat("  relative gradient :", format(rel_g, digits = 3), " (want < 1e-5)\n")
  cat("  gradient by block :",
      paste(names(g_by_block), format(g_by_block, digits = 2),
            sep = " = ", collapse = ", "), "\n")
  cat("  pdHess            :", sdr$pdHess, " (want TRUE)\n")
  cat("  objective         :", format(opt$objective, digits = 8),
      " (-log posterior)\n")
  cat("  loglik            :", format(rep$loglik, digits = 8),
      " (likelihood only, for the information criteria)\n")
  cat("  parameters        :", length(opt$par), "\n")
  lt <- opt$par[names(opt$par) == "log_tau"]
  cat("  log_tau           :", paste(format(lt, digits = 3), collapse = ", "), "\n")
  cat("  tau = exp(log_tau):", paste(format(exp(lt), digits = 3), collapse = ", "), "\n")
  cat("  corr(obs, fitted) :", format(corr, digits = 4), "\n")
  cat("  relative error    :", format(rel_err, digits = 4), "\n")
  cat("  time              :", el, "min\n")
  
  # The two checks that override the convergence code.
  if (rep$loglik > 0)
    cat("  *** WARNING: loglik is positive - the likelihood has overflowed",
        "and this fit is not usable, whatever the convergence code says.\n")
  if (rel_err > 0.5)
    cat("  *** WARNING: relative error is implausible - check the fit.\n")
  cat("\n")
  
  # Everything the later steps need, so none of them has to refit.
  out <- list(
    opt        = opt,
    loglik     = rep$loglik,   # likelihood only - for the information criteria
    mgc        = mgc,
    rel_g      = rel_g,
    g_by_block = g_by_block,
    pdHess     = sdr$pdHess,
    npar       = length(opt$par),
    rel_err    = rel_err,
    corr       = corr,
    log_tau    = opt$par[names(opt$par) == "log_tau"],
    minutes    = as.numeric(el),
    par_full   = par_full,
    log_mu     = rep$log_mu,
    nu         = if (!is.null(rep$nu)) rep$nu else NULL,
    label      = label
  )
  if (keep_sdr) out$sdr <- sdr   # retained for the ZTCMP uncertainty step
  
  saveRDS(out, save_as)
  cat("Saved:", save_as, "\n")
  invisible(out)
}