# =====================================================================
# STEP 1: DATA IMPORT & PREPARATION
# Probabilistic Forecasting of Mortality in Kenya (ZTCMP-BGAM)
# -----
# Source : UN World Population Prospects 2024, Kenya
#          Graduated (model-based) estimates, not raw vital registration.
# Output : mortality_long.rds / .csv (analysis dataset)
#          table1_data_summary.xlsx (Table 1)
# =====================================================================

# ---- 1. Packages ----
# install.packages("readxl")   # run once if not already installed
library(readxl)
library(openxlsx)   # for writing .xlsx tables

# ---- 2. File path ----
file_path <- "Data.xlsx"

# ---- 3. Read raw data ----
# The sheet has TWO header rows:
#   Row 1 = block titles ("Total deaths..." / "Total population...")
#   Row 2 = age labels (0, 1, ..., 100+ repeated twice)
# Skip BOTH (skip = 2); read only the numeric data rows (1950-2023).
raw <- read_excel(file_path, sheet = "Sheet1", col_names = FALSE, skip = 2)
raw <- as.data.frame(raw)
cat("Raw dimensions:", dim(raw)[1], "rows x", dim(raw)[2], "cols\n")

# ---- 4. Define structure ----
#   Col 1       = Year
#   Cols 2:102  = Deaths for ages 0, 1, ..., 99, 100+ (101 ages)
#   Cols 103:203 = Population for ages 0, 1, ..., 99, 100+ (101 ages)
# NOTE: "100+" is an OPEN age interval, coded here as 100.
ages   <- 0:100
years  <- raw[[1]]
n_ages <- length(ages)    # 101
n_yrs  <- length(years)   # 74

deaths_mat <- as.matrix(raw[, 2:102])     # 74 x 101
pop_mat    <- as.matrix(raw[, 103:203])   # 74 x 101

# ---- 5. Build long-format analysis data frame ----
# as.vector() unrolls a matrix column-by-column:
#   col 1 = all years for age 0, col 2 = all years for age 1, ...
# hence Year repeats (times = n_ages) and Age repeats (each = n_yrs).
dat <- data.frame(
  Year     = rep(years, times = n_ages),
  Age      = rep(ages,  each  = n_yrs),
  Deaths   = as.vector(deaths_mat),
  Exposure = as.vector(pop_mat)
)

# The two columns are non-integer for different reasons. Death counts are
# whole numbers stored with tiny floating-point error (e.g. 1006.9999999999999
# is really 1007), so they are ROUNDED to integers - not truncated: as.integer()
# alone truncates toward zero and would return 1006, an off-by-one undercount.
# Exposure (person-years) carries genuine fractional values (many mid-year
# estimates end in .5), so it is retained as numeric rather than rounded.
dat$Deaths <- as.integer(round(dat$Deaths))

# Derived quantities used downstream
dat$Rate   <- dat$Deaths / dat$Exposure   # crude central death rate
dat$logExp <- log(dat$Exposure)           # offset term for the GAM

# Order tidily by Year then Age
dat <- dat[order(dat$Year, dat$Age), ]
rownames(dat) <- NULL

# ---- 6. Data quality checks (validation: all should be zero) ----
cat("\n================ DATA QUALITY CHECKS ================\n")
checks <- data.frame(
  Check = c("Rows equal 74 x 101 = 7474",
            "Missing (NA) in Deaths",
            "Missing (NA) in Exposure",
            "Negative death counts",
            "Negative exposure",
            "Cells where Deaths > Exposure"),
  Result = c(nrow(dat) == 7474,
             sum(is.na(dat$Deaths)),
             sum(is.na(dat$Exposure)),
             sum(dat$Deaths < 0),
             sum(dat$Exposure < 0),
             sum(dat$Deaths > dat$Exposure))
)
print(checks, row.names = FALSE)
cat("=====================================================\n")

# ---- 7. Table 1: Data summary (description) ----
table1 <- data.frame(
  Characteristic = c("Number of age-year cells",
                     "Year range",
                     "Age range",
                     "Deaths: minimum",
                     "Deaths: maximum",
                     "Deaths: mean",
                     "Deaths: median",
                     "Exposure: minimum",
                     "Exposure: maximum",
                     "Cells with zero deaths",
                     "Missing values",
                     "Count dynamic range (max/min)"),
  Value = c(format(nrow(dat), big.mark = ","),
            paste(min(dat$Year), "-", max(dat$Year)),
            paste0(min(dat$Age), " - ", max(dat$Age), "+"),  # 100 is the open interval 100+
            format(min(dat$Deaths), big.mark = ","),
            format(max(dat$Deaths), big.mark = ","),
            format(round(mean(dat$Deaths)),   big.mark = ","),
            format(round(median(dat$Deaths)), big.mark = ","),
            format(round(min(dat$Exposure)),  big.mark = ","),
            format(round(max(dat$Exposure)),  big.mark = ","),
            sum(dat$Deaths == 0),
            sum(is.na(dat)),
            paste0(format(round(max(dat$Deaths) / min(dat$Deaths)),
                          big.mark = ","), "x"))
)
cat("\n==================== TABLE 1: DATA SUMMARY ====================\n")
print(table1, row.names = FALSE)
cat("==============================================================\n")

# Two rows carry the key methodological arguments:
#   "Cells with zero deaths = 0"  -> zeros are structurally absent,
#      justifying the zero-truncated formulation.
#   "Count dynamic range"         -> the span of counts underlies the
#      windowed normalising constant of the CMP likelihood.

# ---- 8. Preview ----
cat("\nFirst 6 rows:\n"); print(head(dat))
cat("\nLast 6 rows:\n");  print(tail(dat))
cat("\nStructure:\n");    str(dat)

# ---- 9. Save ----
saveRDS(dat, "mortality_long.rds")
write.csv(dat, "mortality_long.csv", row.names = FALSE)
write.xlsx(table1, "table1_data_summary.xlsx")
cat("\nSaved: mortality_long.rds, mortality_long.csv, table1_data_summary.xlsx\n")