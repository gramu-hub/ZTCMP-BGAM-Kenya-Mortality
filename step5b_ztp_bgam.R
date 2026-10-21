# =====================================================================
# STEP 5B: 2. ZTP-BGAM
# Probabilistic Forecasting of Mortality in Kenya (ZTCMP-BGAM)
# ---------------------------------------------------------------------
# Run in a fresh R session (one model per session: TMB accumulates
# memory across MakeADFun objects and is sensitive to several compiled
# DLLs loaded at once). Order: 5a -> 5b -> 5c -> 5d -> 5e.
#
# As the Poisson model but with zero truncation; dispersion fixed at
# nu = 1.
#
# Input  : design_matrices.rds, step5_common.R
# Output : ztp_bgam.cpp, ztp_bgam_fit.rds
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
  // Zero-truncated Poisson: Poisson log-density minus log P(Y > 0),
  // where P(Y > 0) = 1 - e^-mu. The guard caps log_p_gt0 below 0 to
  // avoid log(0) if mu were near zero. loglik (likelihood only) is
  // reported separately for the information criteria; nll is the
  // negative log posterior.
  Type loglik = Type(0.0);
  for (int i = 0; i < Y.size(); i++) {
    Type mu_i = exp(log_mu(i));
    Type pois = Type(Y(i)) * log_mu(i) - mu_i - lgamma(Type(Y(i)) + Type(1.0));
    Type log_p_gt0 = log(Type(1.0) - exp(-mu_i));
    if (log_p_gt0 > Type(-1e-10)) log_p_gt0 = Type(-1e-10);
    loglik += pois - log_p_gt0;
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
', "ztp_bgam.cpp")

compile("ztp_bgam.cpp")
dyn.load(dynlib("ztp_bgam"))

# ---- Fit ------------------------------------------------------------
fit <- fit_model(dll     = "ztp_bgam",
                 params  = par_fixed,
                 data    = data_fixed,
                 label   = "2. ZTP-BGAM",
                 save_as = "ztp_bgam_fit.rds")

# ---- Release --------------------------------------------------------
dyn.unload(dynlib("ztp_bgam"))
gc()
cat("\nDone. Next: restart R, then run step5c_cmp_bgam.R\n")