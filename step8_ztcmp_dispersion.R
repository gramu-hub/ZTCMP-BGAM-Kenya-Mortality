# =====================================================================
# STEP 8: DISPERSION ASSESSMENT FROM THE FITTED ZTCMP (ZTCMP ONLY)
# Probabilistic Forecasting of Mortality in Kenya (ZTCMP-BGAM)
# ---------------------------------------------------------------------
# Run in a fresh session. Reads the saved ZTCMP and Poisson fits; refits
# nothing.
#
# Where the Poisson step diagnosed dispersion from residuals, this step
# reports it as estimated by the proposed model through nu_x, and
# compares the two on a common variance-to-mean scale. The CMP nu is
# identical (truncation inert), so it is not duplicated. nu at age 0 is
# fixed and shown as a boundary case.
#
# Classification threshold (same as the Poisson step, on the
# variance-to-mean scale):
#   over  if 1/nu > 1.053  (nu < 0.95)
#   under if 1/nu < 0.962  (nu > 1.04)
#   equi  otherwise
#
# Input  : ztcmp_bgam_fit.rds, poisson_bgam_fit.rds, design_matrices.rds
# Output : fig6_nu_by_age.png, fig7_overlay_poisson_vs_ztcmp.png,
#          table6_dispersion_agreement.xlsx,
#          appendix_ztcmp_nu_by_age.xlsx
# =====================================================================

library(openxlsx)
library(ggplot2)

design <- readRDS("design_matrices.rds")
ztcmp  <- readRDS("ztcmp_bgam_fit.rds")
pois   <- readRDS("poisson_bgam_fit.rds")

ages    <- 0:100
nu      <- ztcmp$nu
vmr_cmp <- 1 / nu          # CMP implied variance-to-mean ratio

OVER_T  <- 1.053
UNDER_T <- 0.962
class_vmr <- function(v) ifelse(v > OVER_T, "over",
                                ifelse(v < UNDER_T, "under", "equi"))
verdict_cmp <- class_vmr(vmr_cmp)
verdict_cmp[ages == 0] <- "over (fixed)"   # age 0 is imposed, flag it

# ---- Figure 6: nu by age ----
# Plotted on nu itself (the model parameter), with the nu = 1 line and
# the over/under-dispersed regions shaded. Age 0 marked as fixed.
dfnu <- data.frame(Age = ages, nu = nu,
                   region = ifelse(nu < 0.95, "over",
                                   ifelse(nu > 1.04, "under", "equi")))
dfnu$fixed <- ages == 0
p6 <- ggplot(dfnu, aes(Age, nu)) +
  annotate("rect", xmin = -2, xmax = 102, ymin = 1.04, ymax = Inf,
           alpha = 0.06, fill = "steelblue") +
  annotate("rect", xmin = -2, xmax = 102, ymin = -Inf, ymax = 0.95,
           alpha = 0.06, fill = "firebrick") +
  geom_hline(yintercept = 1, colour = "grey30", linetype = "dashed") +
  geom_point(aes(colour = region, shape = fixed), size = 1.8) +
  scale_colour_manual(values = c(over = "firebrick", equi = "grey50",
                                 under = "steelblue"), guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     labels = c(`FALSE` = "estimated", `TRUE` = "fixed (age 0)"),
                     name = NULL) +
  scale_x_continuous(breaks = seq(0, 100, by = 5),
                     labels = c(as.character(seq(0, 95, by = 5)), "100+")) +
  scale_y_continuous(breaks = seq(0, 7, by = 0.5)) +
  labs(x = "Age (single years, 100+ open interval)", y = expression(nu[x])) +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank())
ggsave("fig6_nu_by_age.png", p6, width = 12, height = 6, dpi = 300)

# ---- Poisson diagnostic ratio (the Poisson-step residual VMR) ----
mu_pois   <- exp(pois$log_mu)
pearson2  <- (design$Y - mu_pois)^2 / mu_pois
disp_pois <- tapply(pearson2, design$Age, mean)   # empirical VMR
verdict_pois <- class_vmr(as.numeric(disp_pois))

# ---- Figure 7: overlay on the common VMR scale ----
# Both methods on the variance-to-mean scale: Poisson residual ratio
# (diagnosed) against 1/nu_x (modelled).
dfov <- rbind(
  data.frame(Age = ages, vmr = as.numeric(disp_pois), method = "Poisson-BGAM diagnosed"),
  data.frame(Age = ages, vmr = vmr_cmp, method = "ZTCMP-BGAM modelled")
)
p7 <- ggplot(dfov, aes(Age, vmr, colour = method)) +
  geom_hline(yintercept = 1, colour = "grey30", linetype = "dashed") +
  geom_line(alpha = 0.5) + geom_point(size = 1.2) +
  scale_y_log10(breaks = 10^seq(-1, 3, by = 0.5),
                labels = as.character(seq(-1, 3, by = 0.5))) +
  scale_x_continuous(breaks = seq(0, 100, by = 5),
                     labels = c(as.character(seq(0, 95, by = 5)), "100+")) +
  scale_colour_manual(values = c("Poisson-BGAM diagnosed" = "darkorange",
                                 "ZTCMP-BGAM modelled" = "purple4"), name = NULL) +
  labs(x = "Age (single years, 100+ open interval)",
       y = "variance-to-mean ratio (log10 scale)")
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")
ggsave("fig7_overlay_poisson_vs_ztcmp.png", p7, width = 12, height = 6, dpi = 300)

# ---- Table 6: agreement ----
# Counts from each method and the per-age cross-classification. 
v_cmp_plain <- class_vmr(vmr_cmp)   # without the "(fixed)" label
keep <- ages >= 0
counts <- data.frame(
  Method = c("Poisson-BGAM diagnosed", "ZTCMP-BGAM modelled"),
  Over  = c(sum(verdict_pois == "over"),  sum(v_cmp_plain == "over")),
  Equi  = c(sum(verdict_pois == "equi"),  sum(v_cmp_plain == "equi")),
  Under = c(sum(verdict_pois == "under"), sum(v_cmp_plain == "under"))
)
cat("=============== TABLE 6a: DISPERSION COUNTS (of 101) ===============\n")
print(counts, row.names = FALSE)
cat("\n--- per-age cross-classification (all ages 0 to 100+) ---\n")
xtab <- table(`Poisson-BGAM` = verdict_pois[keep], `ZTCMP-BGAM` = v_cmp_plain[keep])
print(xtab)
agree <- sum(verdict_pois[keep] == v_cmp_plain[keep])
cat("\n ages classified identically:", agree, "of", sum(keep), "\n")
logcor <- cor(log(as.numeric(disp_pois)[keep]), log(vmr_cmp[keep]))
cat(" log-correlation (Poisson ratio vs 1/nu_x):", round(logcor, 3), "\n")
cat(" (positive = the two methods agree on the pattern)\n")

wb <- createWorkbook()
addWorksheet(wb, "counts"); writeData(wb, "counts", counts)
addWorksheet(wb, "crosstab"); writeData(wb, "crosstab", as.data.frame.matrix(xtab),
                                        rowNames = TRUE)
saveWorkbook(wb, "table6_dispersion_agreement.xlsx", overwrite = TRUE)

# ---- Appendix: full nu by age ----
# The terminal age (100) is the open interval 100+.
age_label <- as.character(ages)
age_label[ages == max(ages)] <- paste0(max(ages), "+")
appendix <- data.frame(
  Age = age_label,
  `Dispersion parameter (nu_x)` = round(nu, 4),
  implied_VMR = round(vmr_cmp, 4),
  Classification = verdict_cmp,
  check.names = FALSE
)
write.xlsx(appendix, "appendix_ztcmp_nu_by_age.xlsx")

cat("\nSaved: fig6_nu_by_age.png, fig7_overlay_poisson_vs_ztcmp.png,",
    "table6_dispersion_agreement.xlsx, appendix_ztcmp_nu_by_age.xlsx\n")