## Scores the canonical Iteration1 density candidates (bn = 1-3arc-Iteration1-PA1:1)
## against the 2026 field census (data/GroundTruthing/2026_groundtruth_points.gpkg) -
## genuinely new plots at new sites, not a resample of the Iteration1 training data,
## so this is real held-out MAE rather than split-strategy noise.
## Extraction pattern (raster -> points, CRS) copied from
## Second_groundtruth/digitize_evaluate_2026_groundtruth.R.

suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(terra); library(brms); library(Metrics)
  library(parsnip); library(workflows); library(xgboost); library(bonsai); library(lightgbm)
})

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')
source('functions.R')
source('Census_uncertainty/density_bayes.R')

bn <- '1-3arc-Iteration1-PA1:1'
fp <- file.path('..', 'results', 'count_models')
p2proc <- file.path(PROJ_ROOT, 'data', 'spatial', 'processed')

## 1. Assemble the held-out covariate table -------------------------------

gt <- sf::st_read(file.path(PROJ_ROOT, 'data', 'GroundTruthing', '2026_groundtruth_points.gpkg'),
                   quiet = TRUE) |>
  dplyr::filter(!is.na(Plants))

rast_dat <- rastReader('dem_1-3arc', p2proc)
suit_r <- terra::rast(file.path(PROJ_ROOT, 'results', 'suitability_maps',
                                 paste0(bn, 'DO:0-Seed1021-Pr.tif')))

coords <- sf::st_coordinates(gt)
df <- dplyr::bind_cols(
  Observed = gt$Plants,
  Site = gt$Site,
  dplyr::select(terra::extract(rast_dat, terra::vect(gt)), -ID),
  Pr.SuitHab = terra::extract(suit_r, terra::vect(gt))[, 2]
)
# rast_dat already carries its own Longitude/Latitude raster layers (per-cell
# coordinate surfaces); overwrite them with the point's exact coordinates,
# matching wrapper()'s convention (functions.R) rather than the raster's
# cell-binned approximation.
df$Longitude <- coords[, 1]
df$Latitude <- coords[, 2]

n_before <- nrow(df)
df <- tidyr::drop_na(df)
df$point_id <- seq_len(nrow(df))
message(sprintf('%d of %d field plots retained after covariate extraction (%d dropped - outside raster extent or missing covariate).',
                 nrow(df), n_before, n_before - nrow(df)))

## 2. Score each canonical ML candidate ------------------------------------

score_ml <- function(model_suffix, label){
  f <- file.path(fp, 'models', paste0(bn, '-', model_suffix, '.rds'))
  if(!file.exists(f)) return(NULL)
  fit <- readRDS(f)$Model
  preds <- as.numeric(stats::predict(fit, new_data = df)[[1]])
  data.frame(Model = label, point_id = df$point_id, Site = df$Site, Observed = df$Observed, Predicted = preds)
}

ml_results <- dplyr::bind_rows(
  score_ml('poisson_spat', 'XGB Poisson Spat.'),
  score_ml('poisson', 'XGB Poisson'),
  score_ml('tweedie_spat', 'XGB Tweedie Spat.'),
  score_ml('tweedie', 'XGB Tweedie'),
  score_ml('lgbm-poisson_spat', 'LGBM Poisson Spat.')
)

## 2a. Score conformal interval coverage for the ML candidates -------------
## (census_uncertainty_roadmap.md Phase 1c) - the internal test/CV rows used
## to calibrate these (functions.R's densityModeller()) can't judge real
## coverage (circular), so this is where it's actually checked: genuinely new
## field plots, the same place the brms gp() term's extrapolation failure
## first surfaced.

conformal_level <- 0.9

score_conformal <- function(model_suffix, label){
  f <- file.path(fp, 'models', paste0(bn, '-', model_suffix, '-conformal.rds'))
  if(!file.exists(f)) return(NULL)
  cf <- readRDS(f)
  alpha <- 1 - conformal_level

  split_pred <- as.numeric(stats::predict(cf$split$model, new_data = df)[[1]])
  split_q <- stats::quantile(cf$split$residuals, probs = c(alpha / 2, 1 - alpha / 2),
                              na.rm = TRUE, names = FALSE)

  # Full Barber/Candes/Ramdas/Tibshirani CV+ interval: every pooled
  # out-of-fold (fold model, residual) pair contributes its own order
  # statistic (mu_{-k(i)}(x0) +- R_i), not just one randomly-chosen fold -
  # unlike combine_census.R's predict_density_draw_conformal(), which only
  # needs one fold per pseudo-draw and builds the same pooled quantile via
  # many Monte Carlo draws instead. Cheap here since `df` (field plots) and
  # the fold count are both small.
  fold_preds <- sapply(cf$cv_plus$fold_models, function(m) as.numeric(stats::predict(m, new_data = df)[[1]]))
  mu_minus_k <- fold_preds[, cf$cv_plus$fold_id, drop = FALSE]  # n_plots x n_pooled_residuals
  lower_mat <- sweep(mu_minus_k, 2, cf$cv_plus$residuals, '-')
  upper_mat <- sweep(mu_minus_k, 2, cf$cv_plus$residuals, '+')
  cvplus_lower <- apply(lower_mat, 1, stats::quantile, probs = alpha / 2, na.rm = TRUE)
  cvplus_upper <- apply(upper_mat, 1, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)

  dplyr::bind_rows(
    data.frame(Model = label, Method = 'split', point_id = df$point_id, Site = df$Site,
               Observed = df$Observed,
               lower = pmax(split_pred + split_q[1], 0), upper = split_pred + split_q[2]),
    data.frame(Model = label, Method = 'cv_plus', point_id = df$point_id, Site = df$Site,
               Observed = df$Observed,
               lower = pmax(cvplus_lower, 0), upper = cvplus_upper)
  )
}

conformal_results <- dplyr::bind_rows(
  score_conformal('poisson_spat', 'XGB Poisson Spat.'),
  score_conformal('poisson', 'XGB Poisson'),
  score_conformal('tweedie_spat', 'XGB Tweedie Spat.'),
  score_conformal('tweedie', 'XGB Tweedie'),
  score_conformal('lgbm-poisson_spat', 'LGBM Poisson Spat.')
)

if(!is.null(conformal_results) && nrow(conformal_results) > 0){
  conformal_summary <- conformal_results |>
    dplyr::group_by(Model, Method) |>
    dplyr::summarize(
      nominal_level = conformal_level,
      coverage_pct = mean(Observed >= lower & Observed <= upper) * 100,
      mean_width = mean(upper - lower),
      n = dplyr::n(),
      .groups = 'drop'
    ) |>
    dplyr::arrange(Model, Method)

  cat('\n=== conformal interval coverage/width vs. 2026 field groundtruth (nominal ',
      conformal_level * 100, '%) ===\n', sep = '')
  print(as.data.frame(conformal_summary))
  write.csv(conformal_summary,
            file.path('..', 'results', 'tables', paste0(bn, '-2026groundtruth-conformal.csv')),
            row.names = FALSE)
}

## 2b. Score the k-NN null models ---------------------------------------

# Every 2026 field site is new (not a training population), so
# within-population k-NN has zero same-Lctn_bb training rows for every test
# plot and always falls back to across-population by construction - included
# anyway for a direct comparability check (the two columns should end up
# identical here), with across-population as the meaningful score: unlike
# the gp() term, it structurally cannot blow up when extrapolating to an
# unmodelled population (finding_density_model_extrapolation_limits.md).
train_full_sf <- wrapper('1-3arc-Iteration1-PA1:1DO:0-Seed1021-Pr.tif', return_early = TRUE)
knn_test_sf <- sf::st_as_sf(
  data.frame(Prsnc_All = df$Observed, Lctn_bb = df$Site, Pr.SuitHab = df$Pr.SuitHab,
             Longitude = df$Longitude, Latitude = df$Latitude),
  coords = c('Longitude', 'Latitude'), crs = 32613, remove = FALSE
)

knn_scored <- knn_density_null(train_full_sf, knn_test_sf, k = 5)
knn_results <- dplyr::bind_rows(
  data.frame(Model = 'kNN (within pop)', point_id = df$point_id, Site = df$Site,
             Observed = knn_scored$within_population$Observed, Predicted = knn_scored$within_population$Predicted),
  data.frame(Model = 'kNN (across pop)', point_id = df$point_id, Site = df$Site,
             Observed = knn_scored$across_population$Observed, Predicted = knn_scored$across_population$Predicted)
)

## 3. Score the promoted brms candidate -------------------------------------

brms_family_label <- 'poisson'  # promoted family for this bn per Phase 1's MAE-based selection
brms_fit <- readRDS(file.path(fp, 'models', paste0(bn, '-brms-', brms_family_label, '.rds')))
scaling  <- readRDS(file.path(fp, 'models', paste0(bn, '-brms-', brms_family_label, '-scaling.rds')))

std_newdata <- apply_standardization(dplyr::select(df, -Site, -point_id), scaling$center, scaling$scale)
brms_draws <- brms::posterior_epred(brms_fit, newdata = std_newdata, allow_new_levels = TRUE)
# The approximate Hilbert-space gp() term (density_bayes.R) extrapolates via
# sinusoidal basis functions that are only well-behaved inside the training
# spatial domain's boundary factor. A handful of 2026 field plots sit clearly
# outside that domain (new survey sites), and for those the basis blows up -
# individual posterior draws hit ~1e308, so colMeans() on that column is Inf.
# This isn't a coding bug, it's a genuine extrapolation failure of the GP
# term - flagged explicitly below rather than silently producing Inf/NaN MAE.
brms_preds <- colMeans(brms_draws)

brms_results <- data.frame(Model = paste0('brms ', tools::toTitleCase(brms_family_label)),
                            point_id = df$point_id, Site = df$Site,
                            Observed = df$Observed, Predicted = brms_preds)

n_blown_up <- sum(!is.finite(brms_results$Predicted))
if(n_blown_up > 0){
  message(sprintf(
    'brms %s: GP term produced non-finite predictions at %d/%d field plots (site(s): %s) - these plots sit outside the training spatial domain the approximate gp() term was fit on.',
    tools::toTitleCase(brms_family_label), n_blown_up, nrow(brms_results),
    paste(unique(brms_results$Site[!is.finite(brms_results$Predicted)]), collapse = ', ')
  ))
}

## 4. Combine and report ----------------------------------------------------

all_preds <- dplyr::bind_rows(ml_results, knn_results, brms_results)

summarize_preds <- function(preds){
  preds |>
    dplyr::group_by(Model) |>
    dplyr::summarize(
      MAE = Metrics::mae(Observed, Predicted),
      MSE = Metrics::mse(Observed, Predicted),
      RMSE = Metrics::rmse(Observed, Predicted),
      .groups = 'drop'
    ) |>
    dplyr::arrange(MAE)
}

table_raw <- summarize_preds(all_preds)
cat('=== full 2026 field validation (', dplyr::n_distinct(all_preds$point_id), ' plots; brms Inf reflects GP extrapolation failure, see message above) ===\n', sep = '')
print(as.data.frame(table_raw))

## Cochetopa Dome alone carries far more individuals than every other known
## population combined (see finding_density_model_extrapolation_limits.md) -
## a bare all-site MAE is dominated by it. Report the same table excluding
## Coche/Coche Highway plots so "normal"-population performance isn't masked,
## plus a median-based view: the non-Coche Observed distribution is itself
## heavily right-skewed (median 1 plant, mean ~10.5), so MAE alone there is
## still driven by a handful of higher-count plots, not representative of the
## typical plot - MdAE and hit-rate-within-N give a less mean-dominated read.
non_coche <- dplyr::filter(all_preds, !Site %in% c('Coche', 'Coche Highway'))
table_ex_coche <- non_coche |>
  dplyr::group_by(Model) |>
  dplyr::summarize(
    MAE = Metrics::mae(Observed, Predicted),
    MdAE = stats::median(abs(Observed - Predicted)[is.finite(Predicted)]),
    pct_within_2 = mean(abs(Observed - Predicted)[is.finite(Predicted)] <= 2) * 100,
    pct_within_5 = mean(abs(Observed - Predicted)[is.finite(Predicted)] <= 5) * 100,
    n = dplyr::n(),
    n_finite = sum(is.finite(Predicted)),
    .groups = 'drop'
  ) |>
  dplyr::arrange(MAE)
cat('\n=== ex-Cochetopa comparison (', dplyr::n_distinct(non_coche$point_id),
    ' plots; non-Coche Observed mean ', round(mean(dplyr::distinct(non_coche, point_id, .keep_all = TRUE)$Observed), 2),
    ', median ', stats::median(dplyr::distinct(non_coche, point_id, .keep_all = TRUE)$Observed),
    ' - MAE alone is still mean-dominated by a right-skewed tail, MdAE/hit-rate give the typical-plot picture) ===\n', sep = '')
print(as.data.frame(table_ex_coche))
write.csv(table_ex_coche, file.path('..', 'results', 'tables', paste0(bn, '-2026groundtruth-validation-excoche.csv')), row.names = FALSE)

bad_ids <- unique(all_preds$point_id[!is.finite(all_preds$Predicted)])
table <- table_raw
if(length(bad_ids) > 0){
  trimmed_preds <- dplyr::filter(all_preds, !point_id %in% bad_ids)
  table <- summarize_preds(trimmed_preds)
  cat('\n=== trimmed comparison, excluding the ', length(bad_ids),
      ' plot(s) where any model produced a non-finite prediction (', dplyr::n_distinct(trimmed_preds$point_id), ' plots) ===\n', sep = '')
  print(as.data.frame(table))
}

dir.create(file.path('..', 'results', 'tables'), showWarnings = FALSE)
out_csv <- file.path('..', 'results', 'tables', paste0(bn, '-2026groundtruth-validation.csv'))
write.csv(table_raw, file.path('..', 'results', 'tables', paste0(bn, '-2026groundtruth-validation-raw.csv')), row.names = FALSE)
write.csv(table, out_csv, row.names = FALSE)
write.csv(all_preds, file.path('..', 'results', 'tables', paste0(bn, '-2026groundtruth-predictions.csv')), row.names = FALSE)
message('Written to: ', out_csv)
