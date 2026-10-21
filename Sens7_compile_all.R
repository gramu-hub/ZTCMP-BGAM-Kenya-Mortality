# =====================================================================
# SENSITIVITY (compile): all nu profiles into one Excel file.
# Settings: free + fixed 0.10-0.19. Reads saved .rds - no refitting.
# Classification uses the proposal's VMR thresholds (1.053 / 0.962).
# Terminal age shown as 100+.
# Output: nu_all_settings.xlsx (nu_values, classification, fit_metrics)
# =====================================================================
library(openxlsx)

files <- c(
  free   = "ztcmp_bgam_freeage0_fit.rds",
  "0.10" = "ztcmp_age0_sens_0.10.rds",
  "0.11" = "ztcmp_age0_fine_0.11.rds",
  "0.12" = "ztcmp_age0_fine_0.12.rds",
  "0.13" = "ztcmp_age0_fine_0.13.rds",
  "0.14" = "ztcmp_age0_fine_0.14.rds",
  "0.15" = "ztcmp_bgam_fit.rds",
  "0.16" = "ztcmp_age0_fine_0.16.rds",
  "0.17" = "ztcmp_age0_fine_0.17.rds",
  "0.18" = "ztcmp_age0_fine_0.18.rds",
  "0.19" = "ztcmp_age0_fine_0.19.rds"
)
cls <- function(nu) {
  vmr <- 1 / nu
  ifelse(vmr > 1.053, "over", ifelse(vmr < 0.962, "under", "equi"))
}
age_lab <- as.character(0:100); age_lab[101] <- "100+"

nu_tab  <- data.frame(Age = age_lab, check.names = FALSE)
cls_tab <- data.frame(Age = age_lab, check.names = FALSE)
meta    <- data.frame()

for (lab in names(files)) {
  f <- files[[lab]]
  if (!file.exists(f)) { cat("MISSING:", lab, "->", f, "\n"); next }
  fit <- readRDS(f); nu <- fit$nu
  re  <- if (!is.null(fit$rel_err)) fit$rel_err else NA
  co  <- if (!is.null(fit$corr))    fit$corr    else NA
  nu_tab[[lab]]  <- round(nu, 4)
  cls_tab[[lab]] <- cls(nu)
  meta <- rbind(meta, data.frame(
    setting = lab, corr = round(co, 4), rel_err = round(re, 5),
    regime  = ifelse(!is.na(re) && re > 0.06, "DEGRADED", "stable")
  ))
}

wb <- createWorkbook()
addWorksheet(wb, "nu_values");      writeData(wb, "nu_values", nu_tab)
addWorksheet(wb, "classification"); writeData(wb, "classification", cls_tab)
addWorksheet(wb, "fit_metrics");    writeData(wb, "fit_metrics", meta)
saveWorkbook(wb, "nu_all_settings.xlsx", overwrite = TRUE)
cat("Saved: nu_all_settings.xlsx (nu_values, classification, fit_metrics)\n\n")
print(meta, row.names = FALSE)