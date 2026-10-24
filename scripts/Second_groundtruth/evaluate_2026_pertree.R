## Per-tree 2026 ground-truth evaluation.
##
## For every ranger model under <tree>/models/*.rds in each Results_*_iter1a_high_pa tree,
## predict at the 2026 ground-truth points and write a per-model eval CSV into
## <tree>/tables/, mirroring modeller()'s fit-time eval schema (long form:
## metric/estimator/estimate/std_error, plus resolution/iteration columns).
##
## Writes four files per model, crossing a 270 m near-population distance filter with a
## presence-referenced interpolation filter (gt points inside the model's applicability
## domain under CAST::aoa() with a PRESENCE-only reference cloud - the extrapolation subset
## is the same model asked a question it has no training signal for, so scoring on the
## interp subset isolates "model knew enough to answer"):
##   <tree>/tables/<model>-2026gt.csv              - all gt points
##   <tree>/tables/<model>-2026gt-270m.csv         - gt points within 270 m of a known occ
##   <tree>/tables/<model>-2026gt-interp.csv       - gt points in-domain (Interpolation)
##   <tree>/tables/<model>-2026gt-270m-interp.csv  - both filters applied
##
## Skips any (tree, model) whose four output CSVs already exist, so an interrupted run
## resumes cleanly. The interp AOA call needs the caret `train` object in modelsTune/
## (not the ranger model used for probability predictions), via presenceTrainDI() in
## suitability_functions.R; a model missing from modelsTune/ has its two interp CSVs
## written with NA metrics and n = 0.
##
## Usage (working directory scripts/):
##   Rscript Second_groundtruth/evaluate_2026_pertree.R                   # all 4 trees
##   Rscript Second_groundtruth/evaluate_2026_pertree.R Results_LMPO_iter1a_high_pa
##   Rscript Second_groundtruth/evaluate_2026_pertree.R <tree1> <tree2> ...

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(tidyverse)
  library(ranger)
  library(CAST)
})

p2proj <- Find(dir.exists, c('/media/steppe/hdd/EriogonumColoradenseTaxonomy', '~/Documents/Ecoloradense'))
source(file.path(p2proj, 'scripts', 'functions.R'))  # rastReader(), presenceScores(), presenceTrainDI()
p2proc <- file.path(p2proj, 'data', 'spatial', 'processed')

default_trees <- c('Results_LMPO_iter1a_high_pa',
                    'Results_random_iter1a_high_pa',
                    'Results_spatialknn_iter1a_high_pa',
                    'Results_spatialmapping_iter1a_high_pa')
args  <- commandArgs(trailingOnly = TRUE)
trees <- if(length(args) >= 1) args else default_trees
missing_trees <- trees[!dir.exists(file.path(p2proj, trees))]
if(length(missing_trees))
  stop('tree(s) not found: ', paste(missing_trees, collapse = ', '))

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
  mutate(Presence = as.integer(if_any(c(Prsnc_M, Prsnc_J, Prsnc_S),
                                        ~ replace_na(. > 0, FALSE))))
visited_gt <- st_sf(Presence = visited_gt$Presence, geometry = st_geometry(visited_gt))

gt <- bind_rows(gt |> select(Presence), visited_gt)

known_occ <- st_read(file.path(p2proj, 'data', 'collections', 'occurrences_coloradense',
                                'occurrences.shp'), quiet = TRUE) |>
  st_transform(32613)
gt$dist_known_occ <- as.numeric(
  st_distance(gt, known_occ)[cbind(seq_len(nrow(gt)),
                                   st_nearest_feature(gt, known_occ))])

## 2. Per-resolution predictor extraction (cache; the predictor stack doesn't change with
##    split/PAratio/seed, so one terra::extract() per resolution covers every model).

parse_resolution <- function(model_id) sub('-Iteration.*$', '', model_id)
parse_iteration  <- function(model_id) as.integer(sub('.*-Iteration([0-9]+)-PA.*', '\\1', model_id))

point_covs <- new.env(parent = emptyenv())
get_covs <- function(res_tag){
  if(!exists(res_tag, envir = point_covs, inherits = FALSE)){
    rast_dat <- rastReader(paste0('dem_', res_tag), p2proc)
    assign(res_tag, terra::extract(rast_dat, vect(gt)), envir = point_covs)
  }
  get(res_tag, envir = point_covs, inherits = FALSE)
}

## 3. Score the gt subset selected by `keep`, in modeller()'s long-form eval schema ----

score_subset <- function(pred, keep, res_tag, iter){
  keep <- keep & !is.na(pred) & !is.na(gt$Presence)
  if(sum(keep) == 0 || length(unique(gt$Presence[keep])) < 2)
    return(tibble(metric = 'pr_auc', estimator = 'binary',
                  estimate = NA_real_, std_error = NA_real_,
                  resolution = res_tag, iteration = iter, n = sum(keep)))
  scores <- presenceScores(truth = gt$Presence[keep], prob = pred[keep])
  scores$metrics |>
    mutate(resolution = res_tag, iteration = iter, n = sum(keep))
}

## Classify every gt point as in-domain (Interpolation) under a presence-referenced AOA
## built from the caret `train` object in modelsTune/. Returns a logical vector of length
## nrow(gt); NA where the caret model is missing or AOA can't be computed, propagated as
## FALSE into the interp subset so those models score as n = 0 interp rather than erroring.
classify_interp <- function(tree, model_id, covs){
  caret_path <- file.path(p2proj, tree, 'modelsTune', paste0(model_id, '.rds'))
  if(!file.exists(caret_path)) return(rep(NA, nrow(gt)))
  m <- readRDS(caret_path)
  tdi <- tryCatch(presenceTrainDI(m, LPD = FALSE), error = function(e) NULL)
  if(is.null(tdi)) return(rep(NA, nrow(gt)))
  nd <- covs[, tdi$variables, drop = FALSE]
  a <- tryCatch(CAST::aoa(nd, trainDI = tdi, verbose = FALSE), error = function(e) NULL)
  if(is.null(a)) return(rep(NA, nrow(gt)))
  as.logical(a$AOA)
}

## 4. Walk each tree, write per-model CSVs ---------------------------------------------

for(tree in trees){
  models_dir <- file.path(p2proj, tree, 'models')
  tables_dir <- file.path(p2proj, tree, 'tables')
  dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)
  model_files <- list.files(models_dir, pattern = '\\.rds$', full.names = TRUE)
  message('\n=== ', tree, ' (', length(model_files), ' models) ===')
  if(!length(model_files)) next

  for(i in seq_along(model_files)){
    f <- model_files[i]
    model_id <- sub('\\.rds$', '', basename(f))
    out_all        <- file.path(tables_dir, paste0(model_id, '-2026gt.csv'))
    out_270        <- file.path(tables_dir, paste0(model_id, '-2026gt-270m.csv'))
    out_interp     <- file.path(tables_dir, paste0(model_id, '-2026gt-interp.csv'))
    out_270_interp <- file.path(tables_dir, paste0(model_id, '-2026gt-270m-interp.csv'))
    if(all(file.exists(c(out_all, out_270, out_interp, out_270_interp)))) next

    res_tag <- parse_resolution(model_id)
    iter    <- parse_iteration(model_id)
    covs    <- get_covs(res_tag)
    rf      <- readRDS(f)
    pred    <- predict(rf, data = covs, type = 'response')$predictions[, '1']

    if(!file.exists(out_all))
      write.csv(score_subset(pred, rep(TRUE, nrow(gt)), res_tag, iter),
                out_all, row.names = FALSE)
    if(!file.exists(out_270))
      write.csv(score_subset(pred, gt$dist_known_occ <= 270, res_tag, iter),
                out_270, row.names = FALSE)

    if(!file.exists(out_interp) || !file.exists(out_270_interp)){
      interp <- classify_interp(tree, model_id, covs)
      interp_keep <- !is.na(interp) & interp
      if(!file.exists(out_interp))
        write.csv(score_subset(pred, interp_keep, res_tag, iter),
                  out_interp, row.names = FALSE)
      if(!file.exists(out_270_interp))
        write.csv(score_subset(pred, interp_keep & gt$dist_known_occ <= 270, res_tag, iter),
                  out_270_interp, row.names = FALSE)
    }

    if(i %% 25 == 0 || i == length(model_files))
      message('  ', i, ' / ', length(model_files))
  }
}
