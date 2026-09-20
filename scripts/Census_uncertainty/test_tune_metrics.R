## Quick correctness check of the new tune_metric plumbing (functions.R::
## resolve_tune_metric(), threaded through poiss()/tweed()/gbs()/
## brms_cv_compare()) before committing to a full split-strategy sweep under
## 'huber' and 'poisson' - one replicate each, twinning split, same pattern
## as every other new-capability smoke test this session
## (test_knn_and_nogp.R, test_phase1_density_scoped.R).

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
bn_base <- gsub('DO.*$', '', x)
fp <- file.path('..', 'results', 'count_models')

df <- wrapper(x, return_early = TRUE)
coords <- sf::st_coordinates(df)
df$Longitude <- coords[, 1]
df$Latitude  <- coords[, 2]

for(tm in c('huber', 'poisson')){
  bn_i <- paste0(bn_base, '-twinningRep1-', tm, '-metrictest')
  message('\n===== tune_metric = ', tm, ' =====')
  t0 <- Sys.time()
  res <- densityModeller(
    df, bn = bn_i, fp = fp,
    brms_refit_args = list(full_chains = 2, full_iter = 1000, full_warmup = 500),
    split_method = 'twinning', split_replicate = 1, tune_metric = tm
  )
  message(sprintf('[timing] tune_metric=%s total: %.1fs', tm, as.numeric(difftime(Sys.time(), t0, units = 'secs'))))
  print(res$metrics)
}

message('\nDone.')
