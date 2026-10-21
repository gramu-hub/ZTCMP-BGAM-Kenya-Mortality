# =====================================================================
# STEP 3: MODEL SPECIFICATION AND BASIS CONSTRUCTION
# Probabilistic Forecasting of Mortality in Kenya (ZTCMP-BGAM)
# ---------------------------------------------------------------------
# Builds the design matrix and penalty matrices with mgcv (fit = FALSE).
# Nothing is estimated here; the models are fitted later in TMB.
#
# Model: D_{x,t} ~ ZTCMP(mu_{x,t}, nu_x)
#   log(mu_{x,t}) = log(E_{x,t}) + b0 + f1(x) + f2(t) + f3(x,t)
#   f1 = s(Age), f2 = s(Year) on P-spline bases; f3 = ti(Age, Year)
#   with cubic-regression margins. Age-specific dispersion nu_x, with
#   nu at age 0 fixed and ages 1 to 100+ estimated.
#
# Input  : mortality_long.rds
# Output : design_matrices.rds
# =====================================================================

library(mgcv)

dat <- readRDS("mortality_long.rds")

# ---- 1. Model variables ----
dat$logExp  <- log(dat$Exposure)            # exposure offset
dat$Age_idx <- as.integer(factor(dat$Age))  # 1..101, indexes nu_x
n_ages <- length(unique(dat$Age))

# ---- 2. Basis construction (fit = FALSE builds the basis only) ----
basis_spec <- gam(
  Deaths ~ s(Age,  k = 20, bs = "ps") +                # f1(x)
    s(Year, k = 20, bs = "ps") +                # f2(t)
    ti(Age, Year, k = c(10, 10), bs = "cr"),    # f3(x,t)
  data = dat, fit = FALSE
)

# Split the intercept column (nsdf = 1) from the smooth columns.
np <- basis_spec$nsdf
X_int    <- basis_spec$X[, 1:np, drop = FALSE]                        # b0
X_smooth <- basis_spec$X[, (np + 1):ncol(basis_spec$X), drop = FALSE] # f1,f2,f3
K <- ncol(X_smooth)

cat("================ BASIS DIMENSIONS ================\n")
cat("Observations (age-year cells):", nrow(basis_spec$X), "\n")
cat("Parametric columns (X_int) :", ncol(X_int), " [intercept b0]\n")
cat("Smooth columns (X_smooth)  :", K, "\n\n")

for (s in basis_spec$smooth) {
  cat("  ", s$label, ": columns ", s$first.para - np, "-", s$last.para - np,
      " (", s$last.para - s$first.para + 1, " cols), ",
      length(s$S), " penalty matrix/matrices, rank ",
      paste(s$rank, collapse = "/"), "\n", sep = "")
}
cat("==================================================\n\n")

# ---- 3. Penalty matrices ----
# Each smooth's penalty is padded to full K x K (zeros outside its own
# block). ti() contributes two penalties, giving four in total.
S_list  <- list()
S_rank  <- integer(0)
S_label <- character(0)

for (s in basis_spec$smooth) {
  cols <- (s$first.para - np):(s$last.para - np)
  for (j in seq_along(s$S)) {
    S_pad <- matrix(0, K, K)
    S_pad[cols, cols] <- s$S[[j]]
    S_list[[length(S_list) + 1]] <- S_pad
    S_rank  <- c(S_rank, s$rank[j])
    S_label <- c(S_label, paste0(s$label, "_pen", j))
  }
}
n_pen <- length(S_list)

# Check the padded quadratic form matches the block version.
set.seed(1)
b_test <- rnorm(K)
ok  <- TRUE
idx <- 0
for (s in basis_spec$smooth) {
  cols <- (s$first.para - np):(s$last.para - np)
  for (j in seq_along(s$S)) {
    idx <- idx + 1
    q_block <- as.numeric(t(b_test[cols]) %*% s$S[[j]] %*% b_test[cols])
    q_pad   <- as.numeric(t(b_test) %*% S_list[[idx]] %*% b_test)
    if (abs(q_block - q_pad) > 1e-8) ok <- FALSE
  }
}

cat("================ PENALTY MATRICES ================\n")
cat("Number of smoothing parameters (tau_k):", n_pen, "\n")
for (i in seq_len(n_pen)) {
  cat("  ", S_label[i], ": ", K, "x", K, " padded, rank ", S_rank[i], "\n", sep = "")
}
cat("\nPadding verification (block vs padded quadratic form):",
    if (ok) "PASSED" else "*** FAILED ***", "\n")
cat("Age-specific dispersion parameters (nu_x):", n_ages, "\n")
cat("==================================================\n")

if (!ok) stop("Penalty padding is misaligned - do not proceed.")

# ---- 4. Save the design ----
design <- list(
  X_int    = X_int,
  X_smooth = X_smooth,
  S_list   = S_list,
  S_rank   = S_rank,
  S_label  = S_label,
  n_pen    = n_pen,
  logExp   = dat$logExp,
  Y        = dat$Deaths,
  Age_idx  = dat$Age_idx,
  n_ages   = n_ages,
  Age      = dat$Age,
  Year     = dat$Year,
  smooth_info = lapply(basis_spec$smooth,
                       function(s) list(label = s$label,
                                        cols = (s$first.para - np):(s$last.para - np)))
)
saveRDS(design, "design_matrices.rds")

cat("\nSaved: design_matrices.rds\n")