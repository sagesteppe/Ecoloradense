# Smoke test of Phase 1 (census_uncertainty_roadmap.md) - now the full
# 4-family x 2 CV-mode brms sweep (Poisson/NegBinomial/HurdlePoisson/
# HurdleNegBinomial), letting elpd actually pick a winner instead of assuming
# the roadmap's hurdle-negbinomial default. brms_cv_compare()'s per-candidate
# metrics are now scored against the real held-out `test` set (previously
# in-sample train-on-train, which inflated every brms row relative to the ML
# candidates' genuine held-out MAE/MSE/RMSE - see density_bayes.R). The final
# promote-and-refit is trimmed from the full 4-chain/2000-iter precision to
# 2/1000/500 for speed. Still fits the 5 existing ML candidates (poisson/
# tweedie/lgbm, spatial + non-spatial) as densityModeller() normally does,
# since that comparison is the point.

suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(terra); library(brms)
  library(caret); library(future); library(CAST); library(gstat); library(twinning)
  library(recipes); library(rsample); library(parsnip); library(tune)
  library(finetune); library(dials); library(yardstick); library(workflows)
  library(bonsai); library(lightgbm); library(xgboost)
})

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')
source('functions.R')
source('Census_uncertainty/density_bayes.R')
source('Census_uncertainty/density_conformal.R')

x <- '1-3arc-Iteration1-PA1:1DO:0-Seed1021-Pr.tif'
bn <- gsub('DO.*$', '', x)

df <- wrapper(x, return_early = TRUE)
coords <- sf::st_coordinates(df)
df$Longitude <- coords[, 1]
df$Latitude  <- coords[, 2]

result <- densityModeller(
  df, bn = bn, fp = file.path('..', 'results', 'count_models'),
  brms_refit_args = list(full_chains = 2, full_iter = 1000, full_warmup = 500)
)

message('\n=== comparison table ===')
print(result$metrics)

message('\nDone. Promoted brms model: ', result$brms_promoted$model_path)
