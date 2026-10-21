# =====================================================================
# STEP 5A: 1. POISSON-BGAM
# Probabilistic Forecasting of Mortality in Kenya (ZTCMP-BGAM)
# ---------------------------------------------------------------------
# Run in a fresh R session (one model per session: TMB accumulates
# memory across MakeADFun objects and is sensitive to several compiled
# DLLs loaded at once). Order: 5a -> 5b -> 5c -> 5d -> 5e; the CMP
# models warm-start from this fit, so 5a runs first.
#
# Baseline: no truncation, dispersion fixed at nu = 1.
#
# Input  : design_matrices.rds, step5_common.R
# Output : poisson_bgam.cpp, poisson_bgam_fit.rds
# =====================================================================

source("step5_common.R")

# ---- Write the .cpp -------------------------------------------------
writeLines('
#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() () {

  DATA_IVECTOR(Y);
  DATA_MATRIX(X_int);
  DATA_MATRIX(X_smooth);
  DATA_VECTOR(log_exposure);
  DATA_MATRIX(S_stacked);
  DATA_IVECTOR(S_rank);
  DATA_INTEGER(n_pen);
  DATA_INTEGER(K);

  PARAMETER_VECTOR(alpha);
  PARAMETER_VECTOR(beta_smooth);
  PARAMETER_VECTOR(log_tau);       // smoothing parameters

  Type nll = Type(0.0);
  vector<Type> log_mu = X_int * alpha + X_smooth * beta_smooth + log_exposure;

  // ---- Likelihood ----
  // Poisson log-density. loglik (likelihood only) is reported separately
  // for the information criteria; nll is the negative log posterior.
  Type loglik = Type(0.0);
  for (int i = 0; i < Y.size(); i++) {
    Type mu_i = exp(log_mu(i));
    loglik += Type(Y(i)) * log_mu(i) - mu_i - lgamma(Type(Y(i)) + Type(1.0));
  }
  nll -= loglik;
  REPORT(loglik);
  REPORT(log_mu);

  // ---- P-spline penalty ----
  // beta_k ~ N(0, tau_k^{-1} K_k^{-1}); each block contributes
  // 0.5*r_k*log(tau_k) - 0.5*tau_k*(quadratic form) to the objective.
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

  return nll;
}
', "poisson_bgam.cpp")

compile("poisson_bgam.cpp")
dyn.load(dynlib("poisson_bgam"))

# ---- Fit ------------------------------------------------------------
fit <- fit_model(dll     = "poisson_bgam",
                 params  = par_fixed,
                 data    = data_fixed,
                 label   = "1. POISSON-BGAM",
                 save_as = "poisson_bgam_fit.rds")

# ---- Release --------------------------------------------------------
dyn.unload(dynlib("poisson_bgam"))
gc()
cat("\nDone. Next: restart R, then run step5b_ztp_bgam.R\n")