## Per-tree 2026 ground-truth evaluation.
##
## For every ranger model under <tree>/models/*.rds in each Results_*_iter1a_high_pa tree,
## predict at the 2026 ground-truth points and write a per-model eval CSV into
## <tree>/tables/, mirroring modeller()'s fit-time eval schema (long form:
## metric/estimator/estimate/std_error, plus resolution/iteration columns).
##
## Writes two files per model, so a 270 m near-population subset is scored alongside the
## full point set the same way digitize_evaluate_2026_groundtruth.R does:
##   <tree>/tables/<model>-2026gt.csv       - all 2026 gt points
##   <tree>/tables/<model>-2026gt-270m.csv  - points within 270 m of a known occurrence
##
## Skips any (tree, model) whose both output CSVs already exist, so an interrupted run
## resumes cleanly.
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
})

p2proj <- Find(dir.exists, c('/media/steppe/hdd/EriogonumColoradenseTaxonomy', '~/Documents/Ecoloradense'))
source(file.path(p2proj, 'scripts', 'functions.R'))  # rastReader(), presenceScores()
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

## 3. Score one model, within an optional distance-to-known-occurrence filter ----------

score_model <- function(rf_model, covs, keep_dist, res_tag, iter){
  pred <- predict(rf_model, data = covs, type = 'response')$predictions[, '1']
  keep <- keep_dist & !is.na(pred) & !is.na(gt$Presence)
  if(sum(keep) == 0 || length(unique(gt$Presence[keep])) < 2)
    return(tibble(metric = 'pr_auc', estimator = 'binary',
                  estimate = NA_real_, std_error = NA_real_,
                  resolution = res_tag, iteration = iter, n = sum(keep)))
  scores <- presenceScores(truth = gt$Presence[keep], prob = pred[keep])
  scores$metrics |>
    mutate(resolution = res_tag, iteration = iter, n = sum(keep))
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
    out_all  <- file.path(tables_dir, paste0(model_id, '-2026gt.csv'))
    out_270  <- file.path(tables_dir, paste0(model_id, '-2026gt-270m.csv'))
    if(file.exists(out_all) && file.exists(out_270)) next

    res_tag <- parse_resolution(model_id)
    iter    <- parse_iteration(model_id)
    covs    <- get_covs(res_tag)
    rf      <- readRDS(f)

    if(!file.exists(out_all))
      write.csv(score_model(rf, covs, rep(TRUE, nrow(gt)), res_tag, iter),
                out_all, row.names = FALSE)
    if(!file.exists(out_270))
      write.csv(score_model(rf, covs, gt$dist_known_occ <= 270, res_tag, iter),
                out_270, row.names = FALSE)

    if(i %% 25 == 0 || i == length(model_files))
      message('  ', i, ' / ', length(model_files))
  }
}
