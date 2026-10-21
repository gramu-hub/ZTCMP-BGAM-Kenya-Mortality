# =====================================================================
# STEP 5D: 4. ZTCMP-BGAM (proposed)
# Probabilistic Forecasting of Mortality in Kenya (ZTCMP-BGAM)
# ---------------------------------------------------------------------
# Run in a fresh R session (one model per session: TMB accumulates
# memory across MakeADFun objects and is sensitive to several compiled
# DLLs loaded at once). Order: 5a -> 5b -> 5c -> 5d -> 5e.
#
# The proposed model: zero-truncated with age-specific dispersion (the
# CMP plus truncation). sdreport is kept (keep_sdr = TRUE) for the later
# uncertainty step, which draws the credible intervals from it.
#
# Input  : design_matrices.rds, step5_common.R, poisson_bgam_fit.rds
# Output : ztcmp_bgam.cpp, ztcmp_bgam_fit.rds
# =====================================================================

source("step5_common.R")

# ---- Precompute log factorials ----
# log(j!) for every j the normaliser needs (up to max Y + window_hi).
# Depends only on the data.
max_j <- max(data_free$Y) + data_free$window_hi
data_free$log_fac <- lgamma(seq(0, max_j) + 1)
cat("log_fac precomputed up to j =", max_j, "\n")

# ---- Write the .cpp -------------------------------------------------
writeLines('
#include <TMB.hpp>

// ---------------------------------------------------------------------
// CMP log normalising constant, summed over a window around the observed
// count Y (the mass lies there; a full sum is infeasible as counts reach
// ~78,730). The window is set from Y, which is data, so the objective
// stays smooth in the parameters. The sum is anchored at j = Y for
// numerical stability.
// ---------------------------------------------------------------------
template<class Type>
Type cmp_logZ(int Y, Type log_lambda, Type nu, int window_lo, int window_hi,
              vector<Type> log_fac) {
  int j_start = Y - window_lo;
  if (j_start < 0) j_start = 0;
  int j_end = Y + window_hi;
  Type anchor = Type(Y) * log_lambda - nu * log_fac(Y);
  Type sum_rel = Type(0.0);
  for (int j = j_start; j <= j_end; j++) {
    Type lt = Type(j) * log_lambda - nu * log_fac(j);
    sum_rel += exp(lt - anchor);
  }
  return anchor + log(sum_rel);
}

template<class Type>
Type objective_function<Type>::operator() () {

  DATA_IVECTOR(Y);
  DATA_MATRIX(X_int);
  DATA_MATRIX(X_smooth);
  DATA_VECTOR(log_exposure);
  DATA_IVECTOR(Age_id);            // 1..101, indexes nu_x
  DATA_INTEGER(A);
  DATA_INTEGER(window_lo);
  DATA_INTEGER(window_hi);
  DATA_VECTOR(log_fac);            // log_fac(j) = log(j!)
  DATA_MATRIX(S_stacked);
  DATA_IVECTOR(S_rank);
  DATA_INTEGER(n_pen);
  DATA_INTEGER(K);

  PARAMETER_VECTOR(alpha);
  PARAMETER_VECTOR(beta_smooth);
  PARAMETER_VECTOR(log_tau);       // smoothing parameters
  PARAMETER_VECTOR(gamma_nu);      // nu_x = exp(gamma_nu)

  Type nll = Type(0.0);
  vector<Type> log_mu = X_int * alpha + X_smooth * beta_smooth + log_exposure;
  vector<Type> nu = exp(gamma_nu);

  // ---- Likelihood ----
  // loglik (likelihood only) is reported separately for the information
  // criteria; nll is the negative log posterior.
  Type loglik = Type(0.0);
  for (int i = 0; i < Y.size(); i++) {
    Type nu_i = nu(Age_id(i) - 1);   // R 1-indexed, C++ 0-indexed

    // Mean-to-lambda: the model specifies the mean, but the CMP is in
    // lambda with no closed-form inverse, so use the second-order
    // approximation lambda = (mu + (nu-1)/(2*nu))^nu. The correction
    // term (nu-1)/(2*nu) is retained (it is large when nu is far from 1).
    // The clamp guards base <= 0, which can occur transiently when nu and
    // mu are both small.
    Type base = exp(log_mu(i)) + (nu_i - Type(1.0)) / (Type(2.0) * nu_i);
    if (base < Type(1e-8)) base = Type(1e-8);
    Type log_lambda = nu_i * log(base);
    Type logZ = cmp_logZ(Y(i), log_lambda, nu_i, window_lo, window_hi, log_fac);
    Type cmp_lp = Type(Y(i)) * log_lambda - nu_i * log_fac(Y(i)) - logZ;

    // Zero truncation: for the CMP, P(Y=0) = 1/Z, so log P(Y=0) = -logZ;
    // divide the density by P(Y>0) = 1 - 1/Z. The guards prevent log(0)
    // if Z were ever near 1.
    Type log_P0 = -logZ;
    if (log_P0 > Type(-1e-5)) log_P0 = Type(-1e-5);
    Type log_p_gt0 = log(Type(1.0) - exp(log_P0));
    if (log_p_gt0 > Type(-1e-10)) log_p_gt0 = Type(-1e-10);
    loglik += cmp_lp - log_p_gt0;
  }
  nll -= loglik;
  REPORT(loglik);
  REPORT(log_mu);
  REPORT(nu);
  ADREPORT(nu);

  // ---- P-spline penalty ----
  // beta_k ~ N(0, tau_k^{-1} K_k^{-1}); each block contributes
  // 0.5*r_k*log(tau_k) - 0.5*tau_k*(beta^T K_k beta) to the objective.
  for (int k = 0; k < n_pen; k++) {
    matrix<Type> Sk = S_stacked.block(k * K, 0, K, K);
    Type quad = (beta_smooth.matrix().transpose() * Sk *
                 beta_smooth.matrix())(0, 0);
    nll -= Type(0.5) * Type(S_rank(k)) * log_tau(k)
         - Type(0.5) * exp(log_tau(k)) * quad;
  }

  // ---- Priors (identical across all four models) ----
  nll -= sum(dnorm(alpha, Type(-2.6), Type(2.0), true));   // b0
  nll -= sum(dnorm(log_tau, Type(0.0), Type(2.0), true));  // hyperprior on tau

  // nu_x = exp(gamma_nu): prior on the log scale, centred on nu = 1
  // (Poisson equivalence). N(0, 2) admits roughly nu in (0.02, 55), wide
  // enough for the dispersion the data exhibit while shrinking toward
  // Poisson without evidence.
  nll -= sum(dnorm(gamma_nu, Type(0.0), Type(2.0), true));

  return nll;
}
', "ztcmp_bgam.cpp")

compile("ztcmp_bgam.cpp")
dyn.load(dynlib("ztcmp_bgam"))

# ---- Warm start from the fitted Poisson ----
# Start at the converged Poisson solution with gamma_nu = 0 (nu_x = 1),
# so the optimiser only has to find the dispersion.
if (!file.exists("poisson_bgam_fit.rds"))
  stop("Run 5a (Poisson) first - it supplies the starting values.")
pf <- readRDS("poisson_bgam_fit.rds")$par_full
par_free <- list(
  alpha       = unname(pf[names(pf) == "alpha"]),
  beta_smooth = unname(pf[names(pf) == "beta_smooth"]),
  log_tau     = unname(pf[names(pf) == "log_tau"]),
  gamma_nu    = rep(0, data_free$A)
)

# ---- Fit ------------------------------------------------------------
fit <- fit_model(dll          = "ztcmp_bgam",
                 params       = par_free,
                 data         = data_free,
                 label        = "4. ZTCMP-BGAM (proposed)",
                 save_as      = "ztcmp_bgam_fit.rds",
                 keep_sdr     = TRUE,
                 fix_nu_ages  = 1,      # age 0 (index 1); estimate all others
                 fix_nu_value = 0.15)

# ---- Release --------------------------------------------------------
dyn.unload(dynlib("ztcmp_bgam"))
gc()
cat("\nDone. Next: restart R, then run step5e_tables_figures.R\n")