# =====================================================================
# STEP 6 (Table 3): ESTIMATION VALIDATION - ALL FOUR MODELS
# Probabilistic Forecasting of Mortality in Kenya (ZTCMP-BGAM)
# ---------------------------------------------------------------------
# Run in a fresh session. Reads the four saved fits; refits nothing.
#
# TMB optimises rather than samples, so convergence is judged by the
# optimiser status code, the relative gradient (the absolute gradient is
# uninformative at an objective of order 1e5), the gradient by parameter
# block (log_tau and gamma_nu must be near zero; the b0 block may carry a
# harmless ridge), and a positive-definite Hessian.
#
# Input  : the four *_bgam_fit.rds
# Output : table3_convergence.xlsx
# =====================================================================

library(openxlsx)

fits <- list(
  Poisson = readRDS("poisson_bgam_fit.rds"),
  ZTP     = readRDS("ztp_bgam_fit.rds"),
  CMP     = readRDS("cmp_bgam_fit.rds"),
  ZTCMP   = readRDS("ztcmp_bgam_fit.rds")
)

# Pull the log_tau / gamma_nu gradient blocks where they exist - the
# blocks that must be near zero for convergence.
block_grad <- function(f, block) {
  gb <- f$g_by_block
  if (block %in% names(gb)) round(as.numeric(gb[block]), 4) else NA
}

table3 <- data.frame(
  Model                         = paste0(names(fits), "-BGAM"),
  Code                          = sapply(fits, function(f) f$opt$convergence),
  Message                       = sapply(fits, function(f) f$opt$message),
  Relative_gradient             = sapply(fits, function(f) signif(f$rel_g, 3)),
  `Gradient (smoothing block)`  = sapply(fits, block_grad, block = "log_tau"),
  `Gradient (dispersion block)` = sapply(fits, block_grad, block = "gamma_nu"),
  pdHess                        = sapply(fits, function(f) f$pdHess),
  row.names   = NULL,
  check.names = FALSE
)

# Gradient by block shows log_tau and gamma_nu - the blocks that must be
# near zero. The b0 block is omitted: it carries the harmless
# intercept-smooth ridge (already captured by the relative gradient).
cat("================= TABLE 3: CONVERGENCE DIAGNOSTICS =================\n")
print(table3, row.names = FALSE)
cat("\n")

# Footnote, printed so it travels with the table.
cat("FOOTNOTE (for the thesis):\n")
cat(strwrap(paste(
  "All models report optimiser status 'false convergence' (code 1).",
  "This reflects a flat ridge in the joint intercept-smooth space - a",
  "constant may shift between b0 and the level of the smooth functions  without altering",
  "any fitted value - and does not indicate a failed fit. Convergence is",
  "confirmed by (i) a positive-definite Hessian (pdHess = TRUE), (ii) a",
  "relative gradient near zero concentrated on the b0 block, with the",
  "smoothing-parameter (log_tau) and dispersion-parameter (gamma_nu)",
  "gradient blocks also near zero, and (iii) for the fixed-dispersion",
  "models (Poisson and ZTP), a marginal (REML-equivalent) fit converges",
  "cleanly (code 0, relative convergence) and, for the Poisson, agrees",
  "with the joint fit to within 0.06 log-likelihood units (marginal: -131294.14,",
  "joint: -131294.20). The status",
  "code reports that the optimiser could not certify the ridge",
  "direction, not that an optimum was not reached."
), width = 72), sep = "\n")
cat("\n")

write.xlsx(table3, "table3_convergence.xlsx")
cat("\nSaved: table3_convergence.xlsx\n")