README
======
Project : Probabilistic Forecasting of Mortality in Kenya Using a
          Zero-Truncated Conway-Maxwell-Poisson Bayesian Generalized
          Additive Model (ZTCMP-BGAM)
Author  : Grace Mukami Mwangi
          Department of Mathematics
          Multimedia University of Kenya
Contact : graciousmwas3@gmail.com

Software
--------
R version        : 4.6.0
RStudio version  : RStudio-2026.05.1-225
RTools version   : rtools45-6768-6492
Key packages     : TMB, ggplot2, dplyr, openxlsx, mgcv, tidyr, MASS, scales

Data
----
Two files were downloaded from the UN World Population Prospects 2024:

Deaths:
  WPP2024_MORT_F01_1_DEATHS_SINGLE_AGE_BOTH_SEXES.xlsx
  https://population.un.org/wpp/downloads?folder=Standard%20Projections&group=Mortality

Population:
  WPP2024_POP_F01_1_POPULATION_SINGLE_AGE_BOTH_SEXES.xlsx
  https://population.un.org/wpp/downloads?folder=Standard%20Projections&group=Population

Both files were combined in Excel and values multiplied by 1,000 to
convert from thousands to persons, producing Data.xlsx.
step1_data_preparation.R reads Data.xlsx and produces mortality_long.csv and
mortality_long.rds (7,474 age-year cells, Kenya 1950-2023).

How to Run
----------
Run scripts in order: step1 through step11, then sensitivity scripts.
NOTE: Do not run the common setup code block (step5_common)
