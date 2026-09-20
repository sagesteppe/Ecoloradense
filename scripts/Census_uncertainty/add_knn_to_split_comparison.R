## Adds the k-NN density null models (functions.R:knn_density_null(), roadmap
## Phase 1b) to the already-completed split-strategy comparison results
## (Census_uncertainty/compare_split_strategies.R). k-NN is pure computation
## (no fitting), and splitData() caches its train/test partition to disk per
## (method, replicate) - so this reuses the *exact same* held-out sets the
## original ML/Bayesian comparison scored against, without repeating any of
## the expensive RFE/XGBoost/LightGBM/brms fitting.
##
## k fixed at 10 (within-population) / 2 (across-population) - both the
## twinning-replicate-1 sweep's best point (test_knn_and_nogp.R); not
## re-swept per strategy here, same spirit as Kriging using gstat's defaults
## rather than being re-tuned per split.

suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(terra); library(CAST); library(gstat)
})

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')
source('functions.R')

x <- '1-3arc-Iteration1-PA1:1DO:0-Seed1021-Pr.tif'
bn_base <- gsub('DO.*$', '', x)
fp <- file.path('..', 'results', 'count_models')

df <- wrapper(x, return_early = TRUE)
coords <- sf::st_coordinates(df)
df$Longitude <- coords[, 1]
df$Latitude  <- coords[, 2]

score_knn <- function(dsplit, replicate){
  knn_within <- knn_density_null(dsplit$train.sf, dsplit$test, k = 10)$within_population
  knn_across <- knn_density_null(dsplit$train.sf, dsplit$test, k = 2)$across_population
  dplyr::bind_rows(
    mets(knn_within) |> dplyr::mutate(Model = 'kNN (within pop)', .before = 1),
    mets(knn_across) |> dplyr::mutate(Model = 'kNN (across pop)', .before = 1)
  ) |> dplyr::mutate(replicate = replicate, n_test = nrow(dsplit$test))
}

for(method in c('twinning', 'classic', 'spatial_knn', 'population')){
  out_csv <- file.path('..', 'results', 'tables', paste0(bn_base, '-splitcompare-', method, '.csv'))
  if(!file.exists(out_csv)){ message('skip ', method, ' - no existing results'); next }

  existing <- readr::read_csv(out_csv, show_col_types = FALSE)
  reps <- unique(existing$replicate)

  knn_rows <- purrr::map_dfr(reps, function(r){
    bn_i <- paste0(bn_base, '-', method, 'Rep', r)
    dsplit <- splitData(df, fp = fp, bn = bn_i, method = method, replicate = r, k = 10)
    score_knn(dsplit, r)
  })

  combined <- dplyr::bind_rows(
    dplyr::filter(existing, !Model %in% c('kNN (within pop)', 'kNN (across pop)')),
    knn_rows
  )
  write.csv(combined, out_csv, row.names = FALSE)

  summary_tab <- knn_rows |>
    dplyr::filter(Metric == 'MAE') |>
    dplyr::group_by(Model) |>
    dplyr::summarize(mean_mae = mean(Value), sd_mae = stats::sd(Value), n = dplyr::n(), .groups = 'drop')
  cat('\n===', method, '===\n')
  print(as.data.frame(summary_tab))
}

message('\nDone.')
