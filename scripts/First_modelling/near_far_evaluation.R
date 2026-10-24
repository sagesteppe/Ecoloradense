## Evaluate the near / far model pairing on the 2026 field ground truth - the first use of
## those data anywhere in the choice of models.
##   near = the iteration-1 production model of that CV design (PA 1:1, original pool)
##   far  = the high-PA pick for that design and resolution (select_far_ratio.R: one-SE rule
##          on holdout background at 1% prevalence; median replicate at the chosen ratio)
##   switch = near inside the near model's area of applicability (CAST::aoa), far outside
## Each is scored on all 2026 points, and near vs far are compared within bins of the near
## model's dissimilarity index and of distance to the nearest known occurrence - where does
## the far model start to win?
##
## Usage (working directory scripts/): Rscript First_modelling/near_far_evaluation.R
##
## OUTCOME (2026-10-09) - negative; the near/far line of work was dropped here.
##   - 99-100% of 2026 points fall inside the near model's AOA, so the AOA cannot separate
##     near from far on these data and 'switch' is effectively 'near'.
##   - The 'far' picks beat the production near models overall and in most bins, but only
##     because the near models (1:1, original pool) over-predict (Cox intercept -0.3 to -2.2).
##     That is the field-absence effect already shown by highPA_grid.R, not a near/far effect.
##   - > 2.5 km from known occurrences (n = 49) both do poorly and far is worse.
##   - Near and far differ in pool, pipeline and ratio at once, so the comparison is confounded.

suppressPackageStartupMessages({
  library(sf)
  library(tidyverse)
  library(terra)
  library(ranger)
  library(CAST)
})
source('functions.R')
p2proj <- PROJ_ROOT
out_dir <- file.path(PROJ_ROOT, 'results_highPA_fieldabs', 'tables')

## 1. 2026 ground-truth points (as in predict_gridsearch_to_points.R) ----------------------
gt <- read.csv(file.path(p2proj, 'data', 'GroundTruthing', '2026_groundtruth_points.csv')) |>
  mutate(Presence = as.integer(Plants > 0)) |>
  st_as_sf(coords = c('Long', 'Lat'), crs = 4326) |>
  st_transform(32613)
visited <- st_read(file.path(p2proj, 'data', 'GroundTruthPts', 'visited_subset.gpkg'),
                   layer = 'visited_points', quiet = TRUE) |>
  st_transform(st_crs(gt)) |>
  mutate(Presence = as.integer(if_any(c(Prsnc_M, Prsnc_J, Prsnc_S), ~ replace_na(. > 0, FALSE))))
gt <- bind_rows(select(gt, Presence), st_sf(Presence = visited$Presence, geometry = st_geometry(visited))) |>
  filter(!is.na(Presence))
known_occ <- st_read(file.path(p2proj, 'data', 'collections', 'occurrences_coloradense', 'occurrences.shp'),
                     quiet = TRUE) |> st_transform(32613)
gt$dist_known_occ <- as.numeric(st_distance(gt, known_occ[st_nearest_feature(gt, known_occ), ], by_element = TRUE))

## 2. The pairs ----------------------------------------------------------------------------
near <- tibble(
  design = rep(c('spatial', 'classic'), each = 3),
  resolution = rep(c('3arc', '1arc', '1-3arc'), 2),
  near_root = rep(c('results', 'results_classicsplit'), each = 3),
  near_seed = c(1017, 1059, 1021, 1013, 1020, 1061)
)
far <- read.csv(file.path(out_dir, 'far-ratio-median-models.csv')) |>
  filter(target_prevalence == 0.01) |>
  transmute(design, resolution, far_ratio = chosen_ratio, far_seed = median_replicate,
            far_root = ifelse(design == 'spatial', 'results_highPA_fieldabs', 'Results_random_iter1a_high_pa'))
pairs <- inner_join(near, far, by = c('design', 'resolution'))

score <- function(truth, p){
  if(length(unique(truth)) < 2) return(tibble(n = length(truth), brier_scaled = NA_real_, cox_intercept = NA_real_,
                                              cox_slope = NA_real_, roc_auc = NA_real_))
  bind_cols(tibble(n = length(truth)),
            presenceScores(truth, p)$wide |> select(brier_scaled, cox_intercept, cox_slope, roc_auc))
}

## 3. Predict, compute the near model's AOA at the points, score -------------------------
covs_by_res <- map(set_names(unique(pairs$resolution)), function(res)
  terra::extract(rastReader(paste0('dem_', res), file.path('..', 'data', 'spatial', 'processed')), vect(gt)))

results <- pmap(pairs, function(design, resolution, near_root, near_seed, far_ratio, far_seed, far_root){
  message(sprintf('%s %s: near seed %d, far 1:%s seed %d', design, resolution, near_seed, far_ratio, far_seed))
  near_f <- fname_build(resolution, 1, 1,         0, near_seed)
  far_f  <- fname_build(resolution, 1, far_ratio, 0, far_seed)
  covs <- covs_by_res[[resolution]]
  ok <- complete.cases(covs)
  d <- tibble(truth = gt$Presence[ok], dist = gt$dist_known_occ[ok])

  d$p_near <- predict(readRDS(file.path(PROJ_ROOT, near_root, 'models', paste0(near_f, '.rds'))), data = covs[ok, ])$predictions[, '1']
  d$p_far  <- predict(readRDS(file.path(PROJ_ROOT, far_root,  'models', paste0(far_f,  '.rds'))), data = covs[ok, ])$predictions[, '1']

  # the near model's trainDI: reuse the one aoaSurfaces() saved if the backfill got there
  tdi_path <- file.path(PROJ_ROOT, near_root, 'modelsTune', paste0(near_f, '-trainDI.rds'))
  near_tune <- readRDS(file.path(PROJ_ROOT, near_root, 'modelsTune', paste0(near_f, '.rds')))
  tdi <- if(file.exists(tdi_path)) readRDS(tdi_path) else CAST::trainDI(near_tune, verbose = FALSE)
  a <- CAST::aoa(covs[ok, ], near_tune, trainDI = tdi, verbose = FALSE)
  d$DI <- a$DI; d$AOA <- a$AOA
  d$p_switch <- ifelse(d$AOA == 1, d$p_near, d$p_far)

  overall <- map_dfr(c(near = 'p_near', far = 'p_far', switch = 'p_switch'), ~ score(d$truth, d[[.x]]), .id = 'model')

  # near vs far within bins - where does far start to win?
  d <- d |> mutate(
    aoa_bin  = ifelse(AOA == 1, 'inside near AOA', 'outside near AOA'),
    di_bin   = cut(DI / tdi$threshold, c(-Inf, 0.5, 1, 2, Inf), labels = c('DI < 0.5 thr', '0.5-1 thr', '1-2 thr', '> 2 thr')),
    dist_bin = cut(dist, c(-Inf, 270, 450, 2500, Inf), labels = c('< 270 m', '270-450 m', '450 m - 2.5 km', '> 2.5 km')))
  binned <- map_dfr(c('aoa_bin', 'di_bin', 'dist_bin'), function(b)
    d |> group_by(bin = .data[[b]]) |>
      group_modify(~ bind_rows(near = score(.x$truth, .x$p_near), far = score(.x$truth, .x$p_far), .id = 'model')) |>
      ungroup() |> mutate(binning = b, n_presence = NA))

  list(overall = mutate(overall, design = design, resolution = resolution, far_ratio = far_ratio,
                        pct_inside_near_aoa = round(100 * mean(d$AOA == 1), 1)),
       binned = mutate(binned, design = design, resolution = resolution),
       points = mutate(d, design = design, resolution = resolution))
})

overall <- map_dfr(results, 'overall')
binned  <- map_dfr(results, 'binned')
write.csv(overall, file.path(out_dir, 'near-far-2026-overall.csv'), row.names = FALSE)
write.csv(binned,  file.path(out_dir, 'near-far-2026-binned.csv'), row.names = FALSE)
write.csv(map_dfr(results, 'points'), file.path(out_dir, 'near-far-2026-points.csv'), row.names = FALSE)

options(width = 200)
message('\n=== 2026 field ground truth: near vs far vs switch ===')
print(as.data.frame(select(overall, design, resolution, far_ratio, pct_inside_near_aoa, model, n, brier_scaled, cox_intercept, cox_slope, roc_auc)),
      digits = 3, row.names = FALSE)
message('\n=== scaled Brier, near vs far, by bin ===')
binned |> select(design, resolution, binning, bin, model, n, brier_scaled) |>
  pivot_wider(names_from = model, values_from = brier_scaled) |>
  as.data.frame() |> print(digits = 3, row.names = FALSE)
