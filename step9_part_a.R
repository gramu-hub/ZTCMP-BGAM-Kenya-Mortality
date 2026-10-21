# =====================================================================
# STEP 9 (a): IN-SAMPLE MODEL COMPARISON AND SELECTION (Objective 3)
# Probabilistic Forecasting of Mortality in Kenya (ZTCMP-BGAM)
# ---------------------------------------------------------------------
# Run in a fresh session. Reads the four saved fits; refits nothing.
#
# Compares the four fitted models and selects one. The 2x2 lets the
# comparison decompose any improvement into its truncation and
# dispersion components. Out-of-sample validation (Table 9, Fig 8) is
# handled separately in the step9b backtesting scripts.
#
# AIC and BIC use the log-likelihood, not the optimiser objective:
# opt$objective is a negative log posterior (likelihood + penalty +
# priors) and is not comparable across models. The likelihood was
# reported separately (REPORT(loglik) in the .cpp) for this step.
#   AIC = 2p - 2*logLik
#   BIC = p*log(N) - 2*logLik   (p = parameters, N = observations)
#
# Input  : the four *_bgam_fit.rds, design_matrices.rds
# Output : table7_master_comparison.xlsx,
#          table8_factorial_decomposition.xlsx
# =====================================================================

library(openxlsx)

design <- readRDS("design_matrices.rds")
Y <- design$Y
N <- length(Y)

pois  <- readRDS("poisson_bgam_fit.rds")
ztp   <- readRDS("ztp_bgam_fit.rds")
cmp   <- readRDS("cmp_bgam_fit.rds")
ztcmp <- readRDS("ztcmp_bgam_fit.rds")
fits  <- list(Poisson = pois, ZTP = ztp, CMP = cmp, ZTCMP = ztcmp)

# ---- helper: AIC / BIC / relative error / correlation ----
ic <- function(f) {
  p  <- f$npar
  ll <- f$loglik
  mu <- exp(f$log_mu)
  data.frame(
    logLik  = round(ll, 2),
    p       = p,
    AIC     = round(2*p - 2*ll, 2),
    BIC     = round(p*log(N) - 2*ll, 2),
    Rel_err = round(mean(abs((Y - mu)/Y)), 4),
    corr    = round(f$corr, 4)
  )
}

# ---- Table 7: master comparison ----
table7 <- do.call(rbind, lapply(fits, ic))
table7 <- cbind(Model = paste0(names(fits), "-BGAM"), table7)
rownames(table7) <- NULL

# Rank by BIC (harsher on the extra 100 dispersion parameters, so the
# more conservative choice).
table7$dBIC <- round(table7$BIC - min(table7$BIC), 2)
table7 <- table7[, c("Model", "logLik", "p", "AIC", "BIC", "dBIC", "Rel_err", "corr")]

cat("=================== TABLE 7: MASTER COMPARISON ===================\n")
print(table7, row.names = FALSE)
cat("\n N =", N, "observations\n")
cat(" dBIC = BIC - min(BIC); 0 marks the selected model.\n\n")
write.xlsx(table7, "table7_master_comparison.xlsx")

# ---- Table 8: factorial decomposition (log-likelihood) ----
# rows = truncation (no/yes), cols = dispersion (fixed/flexible)
LL <- matrix(c(pois$loglik, cmp$loglik,
               ztp$loglik,  ztcmp$loglik),
             nrow = 2, byrow = TRUE,
             dimnames = list(Truncation = c("Not truncated", "Zero-truncated"),
                             Dispersion = c("Fixed (Poisson-BGAM)",
                                            "Flexible (CMP-BGAM)")))
cat("========== TABLE 8: FACTORIAL DECOMPOSITION (log-likelihood) ==========\n")
print(round(LL, 2))
cat("\n")

# The two contrasts.
trunc_fixed <- ztp$loglik   - pois$loglik   # truncation effect, fixed dispersion
trunc_flex  <- ztcmp$loglik - cmp$loglik    # truncation effect, flexible dispersion
disp_notr   <- cmp$loglik   - pois$loglik   # dispersion effect, no truncation
disp_tr     <- ztcmp$loglik - ztp$loglik    # dispersion effect, zero-truncated

cat(" Truncation effect (reading down each column):\n")
cat("   fixed dispersion    : ZTP-BGAM minus Poisson-BGAM =", round(trunc_fixed, 2), "\n")
cat("   flexible dispersion : ZTCMP-BGAM minus CMP-BGAM   =", round(trunc_flex, 2), "\n")
cat(" Dispersion effect (reading across each row):\n")
cat("   no truncation  : CMP-BGAM minus Poisson-BGAM =", round(disp_notr, 2), "\n")
cat("   zero-truncated : ZTCMP-BGAM minus ZTP-BGAM   =", round(disp_tr, 2), "\n\n")

write.xlsx(as.data.frame.matrix(round(LL, 2)),
           "table8_factorial_decomposition.xlsx", rowNames = TRUE)

# ---- Selection ----
sel <- table7$Model[which.min(table7$BIC)]
cat("=================== SELECTION ===================\n")
cat(" Model selected by BIC:", sel, "\n")
cat(" -> carried forward to uncertainty quantification and the forecast.\n\n")

cat("Saved: table7_master_comparison.xlsx, table8_factorial_decomposition.xlsx\n")