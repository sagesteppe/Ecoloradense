## Case-weight smoke test: refit the high-ratio far pick with a reweighting scheme on
## presences, Platt-recalibrate on an unweighted 15% carve-out of Train, and score the
## result on the 2026 field points. One pair only - spatial x 1-3arc - against the
## existing W0 fit the near-far pipeline already selected.
##
##   WKDE - inverse local env-density on PCA scores (envPCWeights)
##   WMAH - Mahalanobis distance from the training-presence centroid, raw covariates
##          (envMahalanobisWeights)
##
## Pair: design = spatial, resolution = 1-3arc, target_prevalence = 0.01
##   -> chosen_ratio 16, median_replicate seed 5001 (far-ratio-median-models.csv)
##
## Writes:
##   results_highPA_fieldabs/models/1-3arc-Iteration1-PA1:16DO:0-Seed5001-{SCHEME}.rds
##   results_highPA_fieldabs/modelsTune/.../-{SCHEME}.rds and -{SCHEME}-platt.rds
##   results_highPA_fieldabs/tables/far-{SCHEME}-smoke-2026.csv
##
## Usage (working directory scripts/):
##   Rscript First_modelling/wkde_smoke.R           # defaults to WKDE
##   Rscript First_modelling/wkde_smoke.R WMAH

suppressPackageStartupMessages({
  library(sf)
  library(tidyverse)
  library(terra)
  library(caret)
  library(ranger)
})
source('functions.R')
source('final_sdm/Partitioner.R')
source('final_sdm/partition_entrypoints.R')

p2proj <- PROJ_ROOT
p2proc <- file.path('..', 'data', 'spatial', 'processed')
p2dat  <- file.path('..', 'data', 'Data4modelling')
results_root <- file.path(PROJ_ROOT, 'results_highPA_fieldabs')
out_dir <- file.path(results_root, 'tables')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

args <- commandArgs(trailingOnly = TRUE)
scheme <- if(length(args) >= 1) args[1] else 'WKDE'
stopifnot('first argument must be one of WKDE, WMAH' = scheme %in% c('WKDE', 'WMAH'))

resolution <- 10         # 1-3arc
ratio      <- 16
seed       <- 5001
iteration  <- 1

## 1. Reconstruct the training data exactly as highPA_grid.R does ----------------------
candidates <- st_read(file.path(p2dat, 'iter1-pa-highPA-candidates.gpkg'), quiet = TRUE)
reviewed   <- st_read(file.path(p2dat, 'iter1-pa-highPA-reviewed.gpkg'), quiet = TRUE)
abs <- bind_rows(
  st_read(file.path(p2dat, 'iter1-pa.gpkg'), quiet = TRUE) |> select(Occurrence, ID),
  reviewed |> select(Occurrence, ID))

pres_file <- c('90' = '90m', '30' = '30m', '10' = '3m')[[as.character(resolution)]]
x <- st_read(file.path(p2dat, paste0(pres_file, '-presence-iter1.gpkg')), quiet = TRUE) |>
  rename(Occurrence = Presenc) |>
  st_as_sf() |>
  mutate(field_abs = Occurrence == 0)
x <- bind_rows(x, mutate(abs, field_abs = FALSE))
x <- filter(x, st_is(x, 'POINT')) |>
  select(Occurrence, field_abs) |>
  mutate(.id = seq_len(n()))
abs_idx <- which(x$Occurrence == 0)
x$abs_k <- NA_integer_
x$abs_k[abs_idx] <- seq_along(abs_idx)

## 2. Reuse the cached spatial split (split 1) ------------------------------------------
split_path <- file.path(results_root, 'tables', paste0(res_string(resolution), '-split.rds'))
stopifnot('cached split not found - run highPA_grid.R first' = file.exists(split_path))
split <- readRDS(split_path)

## 3. Fit the weighted model -----------------------------------------------------------
message(sprintf('Fitting %s: ratio 1:%d, seed %d, resolution %s',
                scheme, ratio, seed, res_string(resolution)))
t0 <- Sys.time()
fit <- distOrder_PAratio_simulator(
  x = x, distOrder = 0, PAratio = ratio, resolution = resolution, seed = seed,
  predict_surface = FALSE, se_prediction = FALSE,
  split = split, fixed_col = 'field_abs',
  results_root = results_root, iteration = iteration, train_split = 0.8, p2proc = p2proc,
  weight_scheme = scheme)
message(sprintf('Fit done in %.1f min.', as.numeric(difftime(Sys.time(), t0, units = 'mins'))))

## 4. Score on 2026 ground truth, W0 vs WKDE -------------------------------------------
gt <- read.csv(file.path(p2proj, 'data', 'GroundTruthing', '2026_groundtruth_points.csv')) |>
  mutate(Presence = as.integer(Plants > 0)) |>
  st_as_sf(coords = c('Long', 'Lat'), crs = 4326) |>
  st_transform(32613)
visited <- st_read(file.path(p2proj, 'data', 'GroundTruthPts', 'visited_subset.gpkg'),
                   layer = 'visited_points', quiet = TRUE) |>
  st_transform(st_crs(gt)) |>
  mutate(Presence = as.integer(if_any(c(Prsnc_M, Prsnc_J, Prsnc_S), ~ replace_na(. > 0, FALSE))))
gt <- bind_rows(select(gt, Presence),
                st_sf(Presence = visited$Presence, geometry = st_geometry(visited))) |>
  filter(!is.na(Presence))
known_occ <- st_read(file.path(p2proj, 'data', 'collections', 'occurrences_coloradense',
                               'occurrences.shp'), quiet = TRUE) |> st_transform(32613)
gt$dist_known_occ <- as.numeric(
  st_distance(gt, known_occ[st_nearest_feature(gt, known_occ), ], by_element = TRUE))

covs <- terra::extract(rastReader(paste0('dem_', res_string(resolution)), p2proc), vect(gt))
ok <- complete.cases(covs)
truth <- gt$Presence[ok]
dist  <- gt$dist_known_occ[ok]

load_pred <- function(tag){
  f <- sprintf('%s-Iteration%d-PA1:%dDO:0-Seed%d%s.rds',
               res_string(resolution), iteration, ratio, seed,
               if(tag == 'W0') '' else paste0('-', tag))
  mdl <- readRDS(file.path(results_root, 'models', f))
  p <- predict(mdl, data = covs[ok, ])$predictions[, '1']
  if(tag != 'W0'){
    platt <- readRDS(file.path(results_root, 'modelsTune',
                               sub('\\.rds$', '-platt.rds', f)))
    p <- platt$apply(p)
  }
  p
}

p_W0  <- load_pred('W0')
p_new <- load_pred(scheme)

score <- function(truth, p){
  if(length(unique(truth)) < 2) return(tibble(n = length(truth), brier_scaled = NA_real_,
                                              cox_intercept = NA_real_, cox_slope = NA_real_, roc_auc = NA_real_))
  bind_cols(tibble(n = length(truth)),
            presenceScores(truth, p)$wide |>
              select(brier_scaled, cox_intercept, cox_slope, roc_auc))
}

overall <- bind_rows(
  mutate(score(truth, p_W0),  model = 'W0',    bin = 'all'),
  mutate(score(truth, p_new), model = scheme,  bin = 'all'))

dist_bin <- cut(dist, c(-Inf, 270, 450, 2500, Inf),
                labels = c('< 270 m', '270-450 m', '450 m - 2.5 km', '> 2.5 km'))
binned <- map_dfr(levels(dist_bin), function(b){
  i <- which(dist_bin == b)
  bind_rows(
    mutate(score(truth[i], p_W0[i]),  model = 'W0',   bin = b),
    mutate(score(truth[i], p_new[i]), model = scheme, bin = b))
})

tab <- bind_rows(mutate(overall, binning = 'overall'),
                 mutate(binned,  binning = 'dist_known_occ')) |>
  mutate(design = 'spatial', resolution = res_string(resolution),
         ratio = ratio, seed = seed)
write.csv(tab, file.path(out_dir, sprintf('far-%s-smoke-2026.csv', scheme)), row.names = FALSE)

options(width = 180)
message(sprintf('\n=== W0 vs %s on 2026 field ground truth (far pick, spatial x 1-3arc) ===', scheme))
print(as.data.frame(select(tab, binning, bin, model, n, brier_scaled, cox_intercept, cox_slope, roc_auc)),
      digits = 3, row.names = FALSE)
