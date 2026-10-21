# =====================================================================
# SENSITIVITY (main-text table): free + 5 representative fixed values,
# with fit quality and how far the OTHER ages (1-100+) move relative to
# the chosen 0.15 model. Reads saved .rds - no refitting.
# Output: table_age0_sensitivity_main.xlsx
# =====================================================================
library(openxlsx)

files <- c(
  free   = "ztcmp_bgam_freeage0_fit.rds",
  "0.10" = "ztcmp_age0_sens_0.10.rds",
  "0.12" = "ztcmp_age0_fine_0.12.rds",
  "0.15" = "ztcmp_bgam_fit.rds",
  "0.17" = "ztcmp_age0_fine_0.17.rds",
  "0.19" = "ztcmp_age0_fine_0.19.rds"
)
ref <- readRDS(files[["0.15"]])$nu[2:101]   # reference profile, ages 1-100+

tab <- data.frame()
for (lab in names(files)) {
  f <- files[[lab]]
  if (!file.exists(f)) { cat("MISSING:", lab, "\n"); next }
  fit <- readRDS(f); nu <- fit$nu
  tab <- rbind(tab, data.frame(
    age0_fixed           = lab,
    corr                 = if (!is.null(fit$corr))    round(fit$corr, 4)    else NA,
    rel_err              = if (!is.null(fit$rel_err)) round(fit$rel_err, 5) else NA,
    regime               = ifelse(!is.null(fit$rel_err) && fit$rel_err > 0.06,
                                  "DEGRADED", "stable"),
    max_shift_ages_1_100 = round(max(abs(nu[2:101] - ref)), 4)
  ))
}
write.xlsx(tab, "table_age0_sensitivity_main.xlsx")
cat("Saved: table_age0_sensitivity_main.xlsx\n\n")
print(tab, row.names = FALSE)
cat("\n max_shift_ages_1_100 = largest change in any age-1..100+ nu relative\n")
cat(" to the 0.15 model; small values mean the age-0 choice leaves the rest\n")
cat(" of the profile essentially unchanged.\n")