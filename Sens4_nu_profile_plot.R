# =====================================================================
# SENSITIVITY (plot): nu_x by age, FREE age-0 vs FIXED-at-0.15.
# Reads nu_free_vs_fixed15.xlsx - no refitting. Terminal age shown 100+.
# Classification flips computed with the proposal's VMR thresholds
# (over if VMR > 1.053, under if VMR < 0.962).
# =====================================================================
library(openxlsx)
library(ggplot2)
library(tidyr)

d <- read.xlsx("nu_free_vs_fixed15.xlsx")
d$AgeNum <- as.numeric(sub("\\+", "", d$Age))     # "100+" -> 100 for plotting

cls <- function(nu) {
  vmr <- 1 / nu
  ifelse(vmr > 1.053, "over", ifelse(vmr < 0.962, "under", "equi"))
}
d$cls_free  <- cls(d$nu_free)
d$cls_fixed <- cls(d$nu_fixed15)
flip_ages   <- d$AgeNum[d$cls_free != d$cls_fixed]
cat("Ages that flip classification (free vs fixed):",
    if (length(flip_ages)) paste(flip_ages, collapse = ", ") else "none", "\n")

prof_cor <- round(cor(d$nu_free, d$nu_fixed15), 4)
long <- pivot_longer(d[, c("AgeNum", "nu_free", "nu_fixed15")],
                     cols = c("nu_free", "nu_fixed15"),
                     names_to = "version", values_to = "nu")
long$version <- factor(long$version, levels = c("nu_fixed15", "nu_free"),
                       labels = c("Age 0 fixed at 0.15 (model used)",
                                  "Age 0 free (diagnostic)"))

xb <- seq(0, 100, by = 5)
xl <- c(as.character(seq(0, 95, by = 5)), "100+")

p <- ggplot(long, aes(AgeNum, nu, colour = version)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  annotate("text", x = 86, y = 1.2, label = "equidispersion",
           size = 3, colour = "grey40") +
  geom_line(linewidth = 0.7) + geom_point(size = 0.9) +
  scale_x_continuous(breaks = xb, labels = xl) +
  scale_y_continuous(breaks = seq(0, 7, by = 0.5)) +
  scale_colour_manual(values = c("Age 0 fixed at 0.15 (model used)" = "#065A82",
                                 "Age 0 free (diagnostic)"          = "#B0413E")) +
  labs(x = "Age (single years, 100+ open interval)",
       y = expression(nu[x]), colour = NULL) +
  coord_cartesian(ylim = c(0, 7)) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "top", plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank())

if (length(flip_ages) > 0) {
  fl <- d[d$AgeNum %in% flip_ages, ]
  p <- p +
    geom_point(data = fl, aes(AgeNum, nu_free),   inherit.aes = FALSE, shape = 1, size = 3.5) +
    geom_point(data = fl, aes(AgeNum, nu_fixed15), inherit.aes = FALSE, shape = 1, size = 3.5) +
    annotate("text", x = fl$AgeNum, y = pmax(fl$nu_free, fl$nu_fixed15) + 0.4,
             label = paste("age", fl$AgeNum), size = 3)
}
ggsave("nu_free_vs_fixed_profile.png", p, width = 12, height = 7, dpi = 300)
cat("Saved: nu_free_vs_fixed_profile.png\n")

p2 <- p + coord_cartesian(xlim = c(0, 25), ylim = c(0, 3)) +
  scale_y_continuous(breaks = seq(0, 3, by = 0.5)) +
  labs(title = NULL, subtitle = NULL,
       x = "Age (single years)")
ggsave("nu_free_vs_fixed_young.png", p2, width = 8, height = 6, dpi = 300)
cat("Saved: nu_free_vs_fixed_young.png\n")