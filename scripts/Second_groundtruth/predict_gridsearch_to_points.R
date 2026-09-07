## Point-prediction comparison of the PA-ratio grid search models (every seed and
## resolution fit during the search, not just the handful promoted to full suitability
## rasters) against the 2026 ground-truth points. Predicting at point locations rather
## than building a full raster surface per model is what makes it feasible to run this
## across the whole grid search (500+ fitted models).
##
## CV structure is recorded per model, inferred from which results tree it was fit under:
##   results/              -> 'spatial knndm'  (CAST::knndm() class-stratified split)
##   results_classicsplit/ -> 'classic'        (caret::createDataPartition() random split)
##   results_spatialblock/ -> 'spatial block'  (not fit yet - see note below)
## A future spatial-block CV needs to block/permute on POPULATION, not individual site -
## sites cluster tightly within a population, so a site-level block would leak population
## identity across folds. That grouping isn't implemented yet, so no 'spatial block' rows
## exist until scripts/First_modelling fits models under a results_spatialblock/ tree.

library(sf)
library(terra)
library(tidyverse)
library(yardstick)
library(ranger)

p2proj <- '/media/steppe/hdd/EriogonumColoradenseTaxonomy'
p2proc <- file.path(p2proj, 'data', 'spatial', 'processed')
source(file.path(p2proj, 'scripts', 'functions.R'))

## 1. Ground-truth points (same construction as digitize_evaluate_2026_groundtruth.R) ----

gt_raw <- read.csv(file.path(p2proj, 'data', 'GroundTruthing', '2026_groundtruth_points.csv')) |>
  rename(Plot = `Plot..if.tissue.`) |>
  mutate(Date = as.Date(Date, format = '%m/%d/%y'),
         Presence = as.integer(Plants > 0))

gt <- gt_raw |>
  st_as_sf(coords = c('Long', 'Lat'), crs = 4326) |>
  st_transform(32613)

visited_gt <- st_read(file.path(p2proj, 'data', 'GroundTruthPts', 'visited_subset.gpkg'),
                       layer = 'visited_points', quiet = TRUE) |>
  st_transform(st_crs(gt)) |>
  mutate(Presence = as.integer(if_any(c(Prsnc_M, Prsnc_J, Prsnc_S), ~ replace_na(. > 0, FALSE))))
visited_gt <- st_sf(Presence = visited_gt$Presence, geometry = st_geometry(visited_gt))

gt <- bind_rows(gt |> select(Presence), visited_gt)

known_occ <- st_read(file.path(p2proj, 'data', 'collections', 'occurrences_coloradense', 'occurrences.shp'),
                      quiet = TRUE) |>
  st_transform(32613)

gt <- gt |>
  mutate(dist_known_occ = as.numeric(st_distance(gt, known_occ)[cbind(seq_len(nrow(gt)), st_nearest_feature(gt, known_occ))]))

## 2. Locate every fitted grid-search model, across both CV trees -----------------------

cv_roots <- c(
  'results'              = 'spatial knndm',
  'results_classicsplit' = 'classic'
  # 'results_spatialblock' = 'spatial block'  # add once population-blocked CV is fit
)

parse_model_id <- function(f) {
  id <- sub('\\.rds$', '', basename(f))
  tibble(
    path = f,
    model = id,
    resolution = sub('-Iteration.*$', '', id),
    iteration = as.integer(sub('.*-Iteration([0-9]+)-PA.*', '\\1', id)),
    pa_ratio = sub('.*-PA(.+)DO:.*', '\\1', id),
    dist_order = sub('.*DO:([^-]+)-Seed.*', '\\1', id),
    seed = as.integer(sub('.*-Seed([0-9]+)$', '\\1', id))
  )
}

model_files <- map_dfr(names(cv_roots), function(root) {
  f <- list.files(file.path(p2proj, root, 'models'), pattern = '\\.rds$', full.names = TRUE)
  bind_cols(map_dfr(f, parse_model_id), root = root, cv_structure = unname(cv_roots[root]))
})

## modeller()'s own holdout evaluation (results_root/tables/<model>.csv) - fit-time
## pr_auc/roc_auc on that fit's internal train/test split, for comparing against how the
## same model does on the independent 2026 ground truth below.
read_orig_eval <- function(root, model) {
  f <- file.path(p2proj, root, 'tables', paste0(model, '.csv'))
  if (!file.exists(f)) return(tibble(orig_pr_auc = NA_real_, orig_roc_auc = NA_real_))
  oe <- read.csv(f)
  tibble(orig_pr_auc = oe$estimate[oe$metric == 'pr_auc'][1],
         orig_roc_auc = oe$estimate[oe$metric == 'roc_auc'][1])
}
model_files <- bind_cols(model_files, pmap_dfr(model_files, function(root, model, ...) read_orig_eval(root, model)))

## 3. Extract point-level predictor values once per resolution, reused by every model ----
## fit at that resolution - the predictor stack doesn't change with split/PAratio/seed,
## so this is the only per-resolution raster work needed (no per-model raster predict()).

res_labels <- c('3arc' = '3 arc-second', '1arc' = '1 arc-second', '1-3arc' = '1/3 arc-second')

point_covs_by_res <- map(unique(model_files$resolution), function(res) {
  rast_dat <- rastReader(paste0('dem_', res), p2proc)
  terra::extract(rast_dat, vect(gt))
}) |> set_names(unique(model_files$resolution))

## 4. Predict each model at the ground-truth points, evaluate PR-AUC --------------------

evaluate_model <- function(path, model, resolution, iteration, pa_ratio, dist_order, seed,
                            cv_structure, orig_pr_auc, orig_roc_auc, within_m = NULL, ...) {
  covs <- point_covs_by_res[[resolution]]
  keep_dist <- if (is.null(within_m)) rep(TRUE, nrow(gt)) else gt$dist_known_occ <= within_m
  rf_model <- readRDS(path)
  pred <- predict(rf_model, data = covs, type = 'response')$predictions[, '1']
  keep <- keep_dist & !is.na(pred)

  base <- tibble(model = model, resolution = unname(res_labels[resolution]), iteration = iteration,
                 pa_ratio = pa_ratio, dist_order = dist_order, seed = seed, cv_structure = cv_structure,
                 orig_pr_auc = orig_pr_auc, orig_roc_auc = orig_roc_auc)

  if (sum(keep) == 0 || length(unique(gt$Presence[keep])) < 2) {
    return(bind_cols(base, tibble(n = sum(keep), pr_auc = NA_real_)))
  }
  df <- tibble(truth = factor(gt$Presence[keep], levels = c(0, 1)), pred = pred[keep])
  bind_cols(base, tibble(n = nrow(df),
    pr_auc = yardstick::pr_auc(df, truth, pred, event_level = 'second')$.estimate))
}

eval_all <- pmap_dfr(model_files, evaluate_model)
eval_270 <- pmap_dfr(model_files, evaluate_model, within_m = 270)

## 5. Aggregate over seed --------------------------------------------------------------
## A single (resolution, cv_structure, iteration, pa_ratio) row is a seed replicate, not
## a distinct model - seed-to-seed noise is often as large as the PA-ratio effect itself
## (see adaptive_PAratio_search() Details), so comparing individual rows overstates how
## much any one of them can be trusted. Summarise across seeds before comparing ratios.
summarise_over_seeds <- function(eval_tbl) {
  eval_tbl |>
    group_by(resolution, cv_structure, iteration, pa_ratio, dist_order) |>
    summarise(n_seeds = n(), n_pts = first(n),
              mean_pr_auc = mean(pr_auc, na.rm = TRUE), sd_pr_auc = sd(pr_auc, na.rm = TRUE),
              min_pr_auc = min(pr_auc, na.rm = TRUE), max_pr_auc = max(pr_auc, na.rm = TRUE),
              mean_orig_pr_auc = mean(orig_pr_auc, na.rm = TRUE), sd_orig_pr_auc = sd(orig_pr_auc, na.rm = TRUE),
              .groups = 'drop') |>
    arrange(desc(mean_pr_auc))
}

summary_all <- summarise_over_seeds(eval_all)
summary_270 <- summarise_over_seeds(eval_270)

dir.create(file.path(p2proj, 'results', 'tables'), showWarnings = FALSE, recursive = TRUE)
write.csv(eval_all, file.path(p2proj, 'results', 'tables', '2026-groundtruth-gridsearch-pointpredict-all.csv'), row.names = FALSE)
write.csv(eval_270, file.path(p2proj, 'results', 'tables', '2026-groundtruth-gridsearch-pointpredict-270m.csv'), row.names = FALSE)
write.csv(summary_all, file.path(p2proj, 'results', 'tables', '2026-groundtruth-gridsearch-summary-all.csv'), row.names = FALSE)
write.csv(summary_270, file.path(p2proj, 'results', 'tables', '2026-groundtruth-gridsearch-summary-270m.csv'), row.names = FALSE)

summary_all
