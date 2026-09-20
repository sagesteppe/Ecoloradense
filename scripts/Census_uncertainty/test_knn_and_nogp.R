## Quick verification of the two things added after the split-strategy
## comparison + real 2026 field validation (finding_density_model_extrapolation_limits.md):
##   1. k-NN density null models (functions.R:knn_density_null(), roadmap Phase 1b) -
##      within-population and across-population variants.
##   2. brms NegBinomial/HurdleNegBinomial WITHOUT the gp() spatial smoother
##      (density_bayes.R:brms_density_formula(use_gp=), brms_cv_compare(no_gp=)) -
##      a fixed-effects-only ablation, to see whether dropping the (fragile)
##      smoother is viable or whether spatial autocorrelation leaks into the
##      residuals as expected given this project's weak covariates.
## Deliberately skips re-running the ML candidates and the existing gp()
## brms families - those are already comprehensively covered by
## compare_split_strategies.R this session. This is scoped to just the two
## new things, on the canonical Iteration1 twinning-replicate-1 split
## (reuses the cached split if present).

suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(terra); library(brms); library(CAST); library(gstat)
})

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')
source('functions.R')
source('Census_uncertainty/density_bayes.R')

x <- '1-3arc-Iteration1-PA1:1DO:0-Seed1021-Pr.tif'
bn <- gsub('DO.*$', '', x)
fp <- file.path('..', 'results', 'count_models')

df <- wrapper(x, return_early = TRUE)
coords <- sf::st_coordinates(df)
df$Longitude <- coords[, 1]
df$Latitude  <- coords[, 2]

dsplit <- splitData(df, fp = fp, bn = bn, method = 'twinning', replicate = 1)
train <- dsplit$train
test  <- dsplit$test
train.sf <- dsplit$train.sf

## 1. k-NN null models -------------------------------------------------------

for(k in c(3, 5, 10)){
  knn <- knn_density_null(train.sf, test, k = k)
  message(sprintf('k=%d  within-pop MAE=%.3f   across-pop MAE=%.3f',
                   k, Metrics::mae(knn$within_population$Observed, knn$within_population$Predicted),
                   Metrics::mae(knn$across_population$Observed, knn$across_population$Predicted)))
}

## 2. brms NegBinomial / HurdleNegBinomial WITHOUT gp() ----------------------

## Lctn_bb is the population grouping id (needed above for knn_density_null()),
## not a covariate - brms_cv_compare()/standardize_covariates() would otherwise
## try to treat it as one. densityModeller() drops it the same way before its
## own brms_cv_compare() call.
train_nogeom <- sf::st_drop_geometry(train) |> dplyr::select(-Lctn_bb)
test_nogeom  <- sf::st_drop_geometry(test) |> dplyr::select(-Lctn_bb)

nogp_cmp <- brms_cv_compare(
  train_nogeom, test_nogeom,
  families = list(
    NegBinomial_NoGP       = brms::negbinomial(),
    HurdleNegBinomial_NoGP = brms::hurdle_negbinomial()
  ),
  no_gp = c('NegBinomial_NoGP', 'HurdleNegBinomial_NoGP')
)

message('\n=== brms, no gp() smoother ===')
print(nogp_cmp$table)

message('\nDone.')
