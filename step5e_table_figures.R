# =====================================================================
# STEP 5E: TABLES AND FIGURES FROM THE FITTED MODELS
# Probabilistic Forecasting of Mortality in Kenya (ZTCMP-BGAM)
# ---------------------------------------------------------------------
# Run after 5a-5d in a fresh session. Reads the four saved fits; refits
# nothing.
#
# Produces:
#   Table 2 - estimated smoothing parameters tau_k
#   Fig 3   - fitted mortality surface (read against Fig 1, observed)
#   Fig 4   - the three estimated smooth terms f1(x), f2(t), f3(x,t)
#
# Input  : design_matrices.rds, the four *_bgam_fit.rds
# Output : table2_smoothing_parameters.xlsx,
#          fig3_fitted_surface.png, fig4_smooth_terms.png
# =====================================================================

library(openxlsx)
library(ggplot2)

design <- readRDS("design_matrices.rds")
pois  <- readRDS("poisson_bgam_fit.rds")
ztp   <- readRDS("ztp_bgam_fit.rds")
cmp   <- readRDS("cmp_bgam_fit.rds")
ztcmp <- readRDS("ztcmp_bgam_fit.rds")

# ---- Table 2: estimated smoothing parameters ----
# tau_k = exp(log_tau_k), reported across the four models.
smooth_labels <- design$S_label   # e.g. s(Age)_pen1, s(Year)_pen1, ...
tau_tab <- data.frame(
  Smoothing_parameter = smooth_labels,
  `Poisson-BGAM` = round(exp(pois$log_tau), 3),
  `ZTP-BGAM`     = round(exp(ztp$log_tau), 3),
  `CMP-BGAM`     = round(exp(cmp$log_tau), 3),
  `ZTCMP-BGAM`   = round(exp(ztcmp$log_tau), 3),
  check.names = FALSE
)
cat("\n========= TABLE 2: SMOOTHING PARAMETERS (tau_k) =========\n")
print(tau_tab, row.names = FALSE)
write.xlsx(tau_tab, "table2_smoothing_parameters.xlsx")

# ---- Figure 3: fitted mortality surface ----
# Fitted surface from the ZTCMP model; read against the observed surface
# (Figure 1). The four models share the mean structure and give near-
# identical surfaces, so only the proposed model is shown.
dat_fit <- data.frame(
  Age = design$Age,
  Year = design$Year,
  log_rate_fitted = as.vector(ztcmp$log_mu) - design$logExp
)

p3 <- ggplot(dat_fit, aes(x = Year, y = Age, fill = log_rate_fitted)) +
  geom_raster() +
  scale_fill_viridis_c(option = "magma") +
  scale_y_continuous(breaks = seq(0, 100, by = 5),
                     labels = c(as.character(seq(0, 95, by = 5)), "100+")) +
  scale_x_continuous(breaks = seq(1950, 2020, by = 5)) +
  labs(x = "Year", y = "Age (single years, 100+ open interval)",
       fill = "log(rate)") +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"))

ggsave("fig3_fitted_surface.png", p3, width = 12, height = 7, dpi = 300)

# ---- Figure 4: the estimated smooth terms ----
# The three estimated smooths shown individually (the standard GAM
# presentation): f1(x) age effect, f2(t) time effect, f3(x,t) age-time
# interaction. Each is reconstructed from its own block of basis columns
# and coefficients; the column boundaries were saved with the design.
beta_hat <- ztcmp$par_full[names(ztcmp$par_full) == "beta_smooth"]
si <- design$smooth_info

cols_f1 <- si[[1]]$cols
f1_vals <- as.vector(design$X_smooth[, cols_f1, drop = FALSE] %*% beta_hat[cols_f1])
df_f1 <- unique(data.frame(Age = design$Age, f1 = f1_vals))
p_f1 <- ggplot(df_f1, aes(Age, f1)) +
  geom_line(linewidth = 0.8, colour = "steelblue") +
  scale_x_continuous(breaks = seq(0, 100, by = 5),
                     labels = c(as.character(seq(0, 95, by = 5)), "100+")) +
  scale_y_continuous(breaks = seq(-2.5, 3, by = 0.5)) +
  labs(title = expression(paste("f"[1], "(x): age effect")),
       x = "Age (single years, 100+ open interval)", y = expression(paste("f"[1], "(x)"))) +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())

cols_f2 <- si[[2]]$cols
f2_vals <- as.vector(design$X_smooth[, cols_f2, drop = FALSE] %*% beta_hat[cols_f2])
df_f2 <- unique(data.frame(Year = design$Year, f2 = f2_vals))
p_f2 <- ggplot(df_f2, aes(Year, f2)) +
  geom_line(linewidth = 0.8, colour = "steelblue") +
  scale_x_continuous(breaks = seq(1950, 2020, by = 5)) +
  scale_y_continuous(breaks = seq(-0.3, 0.5, by = 0.1),
                     labels = function(x) round(x, 1)) +
  labs(title = expression(paste("f"[2], "(t): time effect")),
       x = "Year", y = expression(paste("f"[2], "(t)"))) +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(),
                                        axis.text.x = element_text(angle = 45, hjust = 1))

cols_f3 <- si[[3]]$cols
f3_vals <- as.vector(design$X_smooth[, cols_f3, drop = FALSE] %*% beta_hat[cols_f3])
df_f3 <- data.frame(Age = design$Age, Year = design$Year, f3 = f3_vals)
p_f3 <- ggplot(df_f3, aes(Year, Age, fill = f3)) +
  geom_raster() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_y_continuous(breaks = seq(0, 100, by = 5),
                     labels = c(as.character(seq(0, 95, by = 5)), "100+")) +
  scale_x_continuous(breaks = seq(1950, 2020, by = 5)) +
  labs(title = expression(paste("f"[3], "(x,t): age-time interaction")),
       x = "Year", y = "Age (single years, 100+ open interval)", fill = NULL) +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p4 <- (p_f1 | p_f2) / p_f3 +
    plot_annotation(theme = theme(plot.title = element_text(face = "bold",
                                                            size = 13)))
  ggsave("fig4_smooth_terms.png", p4, width = 12, height = 9, dpi = 300)
  cat("\nSaved: fig4_smooth_terms.png\n")
} else {
  message("patchwork not installed - saving the three panels separately.")
  ggsave("fig4a_f1_age.png", p_f1, width = 8, height = 5, dpi = 300)
  ggsave("fig4b_f2_time.png", p_f2, width = 8, height = 5, dpi = 300)
  ggsave("fig4c_f3_interaction.png", p_f3, width = 12, height = 7, dpi = 300)
}

cat("\nSaved: table2_smoothing_parameters.xlsx, fig3_fitted_surface.png\n")