# Diagnostic-only: reproduce the exact train_brms densityModeller() builds,
# then call brms::brm() directly (no silent=2, no cheap-search wrapper) to see
# the real per-chain error test_phase1_density_scoped.R's silent=2 hid.

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

x <- '1-3arc-Iteration1-PA1:1DO:0-Seed1021-Pr.tif'
bn <- gsub('DO.*$', '', x)
fp <- file.path('..', 'results', 'count_models')

df <- wrapper(x, return_early = TRUE)
coords <- sf::st_coordinates(df)
df$Longitude <- coords[, 1]
df$Latitude  <- coords[, 2]

dsplit <- splitData(df, fp = fp, bn = bn)
train <- dsplit$train

train <- sf::st_drop_geometry(train) |> select(-Lctn_bb)
train_coords <- train[, c('Longitude', 'Latitude')]

rfProfile <- readRDS(file.path(fp, 'modelsTune', paste0(bn, '.rds')))
train <- dplyr::select(train, all_of(c('Prsnc_All', predictors(rfProfile))))

train_brms <- train
train_brms$Longitude <- train_coords[, 1]
train_brms$Latitude  <- train_coords[, 2]

cat('train_brms dim:', dim(train_brms), '\n')
cat('any NA:', anyNA(train_brms), '\n')
cat('Prsnc_All class/range:', class(train_brms$Prsnc_All), range(train_brms$Prsnc_All), '\n')
str(train_brms)

form <- brms_density_formula(train_brms)
cat('\nformula:\n')
print(form)

cat('\n=== fitting verbosely, 1 chain, short ===\n')
fit <- brms::brm(
  form, data = train_brms, family = brms::hurdle_negbinomial(), backend = 'cmdstanr',
  chains = 1, iter = 200, warmup = 100, cores = 1, seed = 1,
  silent = 0, refresh = 10
)
print(fit)
