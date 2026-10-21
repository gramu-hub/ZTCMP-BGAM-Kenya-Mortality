# =====================================================================
# SENSITIVITY (extract): all 101 nu values (ages 0-100+) from the
# FREE-age-0 fit, side by side with the FIXED-0.15 values.
# Reads saved .rds only - no refitting.
# Output: nu_free_vs_fixed15.xlsx  (terminal age shown as 100+)
# =====================================================================
library(openxlsx)

free    <- readRDS("ztcmp_bgam_freeage0_fit.rds")
nu_free <- free$nu
stopifnot(length(nu_free) == 101)

age_lab <- as.character(0:100)
age_lab[101] <- "100+"                 # terminal age is the open interval

out <- data.frame(
  Age      = age_lab,
  nu_free  = round(nu_free, 4),
  VMR_free = round(1 / nu_free, 4),
  check.names = FALSE
)

# Fixed-0.15 fit (the real model); fall back to the sensitivity file.
nu_fixed <- NULL
if (file.exists("ztcmp_bgam_fit.rds")) {
  fx <- readRDS("ztcmp_bgam_fit.rds")
  if (!is.null(fx$nu)) nu_fixed <- fx$nu
}
if (is.null(nu_fixed) && file.exists("ztcmp_age0_sens_0.15.rds")) {
  nu_fixed <- readRDS("ztcmp_age0_sens_0.15.rds")$nu
}
if (!is.null(nu_fixed)) {
  out$nu_fixed15 <- round(nu_fixed, 4)
  out$diff       <- round(nu_free - nu_fixed, 4)
  out$abs_diff   <- round(abs(nu_free - nu_fixed), 4)
}

write.xlsx(out, "nu_free_vs_fixed15.xlsx")
cat("Saved: nu_free_vs_fixed15.xlsx (all 101 ages, terminal = 100+)\n\n")

cat("Young ages 0-15 (where the age-0 choice has most effect):\n")
print(head(out, 16), row.names = FALSE)
if (!is.null(nu_fixed)) {
  cat(sprintf("\n  median |diff| : %.4f\n", median(out$abs_diff)))
  cat(sprintf("  max    |diff| : %.4f  at age %s\n",
              max(out$abs_diff), out$Age[which.max(out$abs_diff)]))
  cat(sprintf("  ages differing by > 0.1 : %d of 101\n", sum(out$abs_diff > 0.1)))
  cat(sprintf("  ages differing by > 0.5 : %d of 101\n", sum(out$abs_diff > 0.5)))
  cat(sprintf("  correlation of the two nu profiles : %.4f\n", cor(nu_free, nu_fixed)))
}