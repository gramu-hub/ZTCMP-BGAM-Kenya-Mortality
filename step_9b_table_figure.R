# =====================================================================
# STEP 9 (b): COLLATE BACKTESTS - Table 9 and Fig 8
# Probabilistic Forecasting of Mortality in Kenya (ZTCMP-BGAM)
# ---------------------------------------------------------------------
# Run after the four backtest scripts. Reads their saved outputs; refits
# nothing.
#
# Produces:
#   Table 9 - out-of-sample forecast relative error, four models
#   Fig 8   - forecast relative error by age, four models, held-out 2014-2023
# plus a diagnostic breakdown of the error by age and by year, to
# show where the differences between models come from.
#
# Input  : backtest_{poisson,ztp,cmp,ztcmp}_bgam.rds
# Output : table9_forecast_error.xlsx,
#          appendix_forecast_error_breakdown.xlsx,
#          fig8_forecast_vs_observed.png
# =====================================================================

library(openxlsx)
library(ggplot2)

files <- c(Poisson = "backtest_poisson_bgam.rds",
           ZTP     = "backtest_ztp_bgam.rds",
           CMP     = "backtest_cmp_bgam.rds",
           ZTCMP   = "backtest_ztcmp_bgam.rds")
bt <- lapply(files, readRDS)

# Per-cell relative error, long form.
long <- do.call(rbind, lapply(names(bt), function(nm) {
  b <- bt[[nm]]
  data.frame(Model    = paste0(nm, "-BGAM"),
             Age      = b$Age,
             Year     = b$Year,
             rate_obs = b$rate_obs,
             rate_hat = b$rate_hat,
             relerr   = abs(b$rate_obs - b$rate_hat) / b$rate_obs)
}))
model_order <- paste0(names(files), "-BGAM")

# ---- Table 9: overall forecast relative error ----
table9 <- data.frame(
  Model = model_order,
  Forecast_relative_error = round(sapply(bt, function(b) b$fre), 4),
  row.names = NULL
)
cat("=============== TABLE 9: OUT-OF-SAMPLE FORECAST ERROR ===============\n")
print(table9, row.names = FALSE)
cat("\n holdout: 2014-2023; lower is better.\n\n")
write.xlsx(table9, "table9_forecast_error.xlsx")

# ---- Diagnostic: error by single age ----
byage_single <- aggregate(relerr ~ Age + Model, data = long, FUN = mean)
byage_single$Model <- factor(byage_single$Model, levels = model_order)
byage_single$Age <- as.character(byage_single$Age)
byage_single$Age[byage_single$Age == "100"] <- "100+"
cat("--- forecast relative error by single age (first 10 rows) ---\n")
print(head(byage_single, 10))
cat("\n")

# ---- Diagnostic: error by holdout year ----
byyear <- round(tapply(long$relerr, list(long$Year, long$Model), mean), 4)
byyear <- byyear[, model_order]
cat("--- forecast relative error by holdout year ---\n")
print(byyear)
cat("\n")

wb <- createWorkbook()
addWorksheet(wb, "by_age")
writeData(wb, "by_age", byage_single)
addWorksheet(wb, "by_year")
writeData(wb, "by_year", as.data.frame.matrix(byyear), rowNames = TRUE)
saveWorkbook(wb, "appendix_forecast_error_breakdown.xlsx", overwrite = TRUE)

# ---- Figure 8: forecast relative error by single year of age ----
# Average the per-cell relative error over the holdout years 2014-2023
# for each age and model. One panel; four lines.

byage_single <- aggregate(relerr ~ Age + Model, data = long, FUN = mean)
byage_single$Model <- factor(byage_single$Model, levels = model_order)

p8 <- ggplot(byage_single, aes(x = Age, y = relerr, colour = Model)) +
  geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c("Poisson-BGAM" = "darkorange",
                                 "ZTP-BGAM"     = "firebrick",
                                 "CMP-BGAM"     = "steelblue",
                                 "ZTCMP-BGAM"   = "purple4"),
                      breaks = model_order, name = NULL) +
  scale_x_continuous(breaks = seq(0, 100, by = 5),
                     labels = c(as.character(seq(0, 95, by = 5)), "100+")) +
                     scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                        breaks = seq(0, 1.5, by = 0.1)) +
  labs(x        = "Age (single years, 100+ open interval)",
       y        = "Mean forecast relative error") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor  = element_blank(),
        plot.title        = element_text(face = "bold"),
        legend.position   = "top")

ggsave("fig8_forecast_vs_observed.png", p8, width = 12, height = 6, dpi = 300)

# ---- Figure 8b: forecast relative error by holdout year ----

byyear_long <- data.frame(
  Year  = as.integer(rownames(byyear)),
  as.data.frame(byyear),
  check.names = FALSE
)

byyear_long <- tidyr::pivot_longer(
  byyear_long,
  cols      = -Year,
  names_to  = "Model",
  values_to = "relerr"
)

byyear_long$Model <- factor(byyear_long$Model, levels = model_order)

p8b <- ggplot(byyear_long, aes(x = Year, y = relerr, colour = Model)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  scale_x_continuous(breaks = seq(2014, 2023, by = 1)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = seq(0, 0.8, by = 0.05)) +
  scale_colour_manual(values = c(
    "Poisson-BGAM"  = "#E69F00",
    "ZTP-BGAM"      = "#CC0000",
    "CMP-BGAM"      = "#56B4E9",
    "ZTCMP-BGAM"    = "#4B0082"
  )) +
  labs(x = "Holdout year",
       y = "Mean forecast relative error",
       colour = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top")

ggsave("fig8b_forecast_error_by_year.png", p8b, width = 10, height = 6, dpi = 300)

cat("Saved: table9_forecast_error.xlsx,",
    "appendix_forecast_error_breakdown.xlsx, fig8_forecast_vs_observed.png\n, fig8b_forecast_error_by_year.png\n")