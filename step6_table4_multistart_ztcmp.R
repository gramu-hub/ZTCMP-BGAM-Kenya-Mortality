# =====================================================================
# STEP 6 (Table 4): MULTI-START ROBUSTNESS - ZTCMP
# Probabilistic Forecasting of Mortality in Kenya (ZTCMP-BGAM)
# ---------------------------------------------------------------------
# Run in a fresh session. Refits the ZTCMP from several starting values
# for the dispersion parameters gamma_nu, to confirm the optimiser
# returns to the same optimum regardless of where it begins (so the
# solution is not an artefact of the Poisson warm start).
#
# Each start uses the same converged Poisson mean (alpha, beta_smooth,
# log_tau) but a different constant nu across ages:
#   nu = 0.5 (overdispersed), 1.0 (equidispersion), 2.0 (underdispersed).
# Age 0 is fixed at 0.15 in every start, as in the main fit. Agreement in
# the objective (-log posterior) and nu profile across starts is the evidence.
#
# Input  : design_matrices.rds, poisson_bgam_fit.rds, step5_common.R,
#          the compiled ztcmp_bgam DLL (built by 5d)
# Output : table4_multistart_ztcmp.xlsx
# =====================================================================

source("step5_common.R")
library(openxlsx)

# log factorials for the normaliser
max_j <- max(data_free$Y) + data_free$window_hi
data_free$log_fac <- lgamma(seq(0, max_j) + 1)

# Load the compiled ztcmp_bgam DLL (built by 5d); this script does not
# recompile it.
if (!file.exists(dynlib("ztcmp_bgam"))) {
  stop("ztcmp_bgam DLL not found - run 5d first so the .cpp is compiled.")
}
dyn.load(dynlib("ztcmp_bgam"))

pf     <- readRDS("poisson_bgam_fit.rds")$par_full
alpha0 <- unname(pf[names(pf) == "alpha"])
beta0  <- unname(pf[names(pf) == "beta_smooth"])
ltau0  <- unname(pf[names(pf) == "log_tau"])

starts <- c(0.5, 1.0, 2.0)   # starting nu (constant across ages)
res <- data.frame()

for (s in starts) {
  cat("\n===== multi-start: nu0 =", s, "=====\n"); flush.console()
  par_free <- list(
    alpha       = alpha0,
    beta_smooth = beta0,
    log_tau     = ltau0,
    gamma_nu    = rep(log(s), data_free$A)
  )
  fit <- fit_model(dll          = "ztcmp_bgam",
                   params       = par_free,
                   data         = data_free,
                   label        = paste0("ZTCMP multi-start nu0=", s),
                   save_as      = paste0("ztcmp_multistart_nu", s, ".rds"),
                   fix_nu_ages  = 1,
                   fix_nu_value = 0.15)
  res <- rbind(res, data.frame(
    `Starting nu`   = s,
    Objective       = round(fit$opt$objective, 2),
    `nu_x (min)`    = round(min(fit$nu[-1]), 3),   # excl. fixed age 0
    `nu_x (max)`    = round(max(fit$nu[-1]), 3),
    `nu_x (median)` = round(median(fit$nu[-1]), 3),
    check.names = FALSE
  ))
}

dyn.unload(dynlib("ztcmp_bgam"))

cat("\n================= TABLE 4 (ZTCMP-BGAM): MULTI-START ROBUSTNESS =================\n")
print(res, row.names = FALSE)
cat("\n")

# Do all starts agree?
obj_spread <- max(res$Objective) - min(res$Objective)
cat("Objective (-log posterior) spread across starts:", round(obj_spread, 4), "\n")
if (obj_spread < 0.01) {
  cat("VERDICT: all starts converge to the same optimum (spread < 0.01",
      "objective units). The solution is robust to the starting value.\n")
} else {
  cat("VERDICT: starts disagree by", round(obj_spread, 4), "objective units -",
      "investigate before reporting.\n")
}

write.xlsx(res, "table4_multistart_ztcmp.xlsx")
cat("\nSaved: table4_multistart_ztcmp.xlsx\n")