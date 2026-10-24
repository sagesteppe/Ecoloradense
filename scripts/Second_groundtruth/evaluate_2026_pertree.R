## Per-tree 2026 ground-truth evaluation, stratified by applicability-domain x distance.
##
## For every ranger model under <tree>/models/*.rds in each Results_*_iter1a_high_pa tree,
## predict at the 2026 ground-truth points and write ONE long-form CSV per model:
##   <tree>/tables/<model>-2026gt.csv
##
## Each file has one row per (metric, SurfaceArea, Distance) cell. The Distance axis is
## split only within the Interpolation subset; a handful of Extrapolation points fall
## within 270 m too but that cell is always too thin to score, so Extrapolation is pooled:
##   SurfaceArea = Interpolation, Distance = <=270m
##   SurfaceArea = Interpolation, Distance = >270m
##   SurfaceArea = Extrapolation, Distance = all
## Three mutually-exclusive subsets of the gt points, each scored with presenceScores()
## (pr_auc, roc_auc, brier, Cox calibration, DCA scalars).
##
## The AOA step needs the caret `train` object in modelsTune/ (not the ranger model used
## for probability predictions), via presenceTrainDI() in suitability_functions.R. If AOA
## can't be computed for a model (missing caret file / trainDI error / aoa() error), the
## SurfaceArea column is left as NA and all four cells report n = 0 - the first failure of
## each kind is messaged so a broken AOA path is visible.
##
## Skips any (tree, model) whose CSV already exists, so an interrupted run resumes cleanly.
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

## 3. Score one (SurfaceArea, Distance) subset in modeller()'s long-form eval schema ----

score_subset <- function(pred, keep, res_tag, iter, surface_area, distance){
  keep <- keep & !is.na(pred) & !is.na(gt$Presence)
  base_cols <- list(resolution = res_tag, iteration = iter, n = sum(keep),
                    SurfaceArea = surface_area, Distance = distance)
  if(sum(keep) == 0 || length(unique(gt$Presence[keep])) < 2)
    return(tibble(metric = 'pr_auc', estimator = 'binary',
                  estimate = NA_real_, std_error = NA_real_, !!!base_cols))
  scores <- presenceScores(truth = gt$Presence[keep], prob = pred[keep])
  bind_cols(scores$metrics, as_tibble(base_cols))
}

## Classify every gt point as in-domain (Interpolation) under a presence-referenced AOA
## built from the caret `train` object in modelsTune/. Returns a logical vector of length
## nrow(gt); NA where the caret model is missing or AOA can't be computed. The first
## failure of each kind per run is messaged so a broken AOA path is visible rather than
## silently filling every row with n = 0.
interp_warned <- new.env(parent = emptyenv())
warn_once <- function(kind, model_id, msg){
  if(exists(kind, envir = interp_warned, inherits = FALSE)) return(invisible())
  assign(kind, TRUE, envir = interp_warned)
  message('  [interp ', kind, '] ', model_id, ': ', msg)
}
classify_interp <- function(tree, model_id, covs){
  caret_path <- file.path(p2proj, tree, 'modelsTune', paste0(model_id, '.rds'))
  if(!file.exists(caret_path)){
    warn_once('missing_caret', model_id, caret_path)
    return(rep(NA, nrow(gt)))
  }
  m <- readRDS(caret_path)
  tdi <- tryCatch(presenceTrainDI(m, LPD = FALSE), error = function(e) e)
  if(inherits(tdi, 'error')){
    warn_once('presenceTrainDI_error', model_id, conditionMessage(tdi))
    return(rep(NA, nrow(gt)))
  }
  missing_vars <- setdiff(tdi$variables, names(covs))
  if(length(missing_vars)){
    warn_once('missing_vars', model_id,
              paste(length(missing_vars), 'tdi$variables missing from covs:',
                    paste(head(missing_vars, 5), collapse = ', ')))
    return(rep(NA, nrow(gt)))
  }
  nd <- covs[, tdi$variables, drop = FALSE]
  a <- tryCatch(CAST::aoa(nd, trainDI = tdi, verbose = FALSE), error = function(e) e)
  if(inherits(a, 'error')){
    warn_once('aoa_error', model_id, conditionMessage(a))
    return(rep(NA, nrow(gt)))
  }
  as.logical(a$AOA)
}

## 4. Walk each tree, one long CSV per model -------------------------------------------

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
    out <- file.path(tables_dir, paste0(model_id, '-2026gt.csv'))
    if(file.exists(out)) next

    res_tag <- parse_resolution(model_id)
    iter    <- parse_iteration(model_id)
    covs    <- get_covs(res_tag)
    rf      <- readRDS(f)
    pred    <- predict(rf, data = covs, type = 'response')$predictions[, '1']

    interp <- classify_interp(tree, model_id, covs)   # TRUE / FALSE / NA
    is_near <- gt$dist_known_occ <= 270
    cells <- list(
      list(mask = !is.na(interp) &  interp & is_near,  sa = 'Interpolation', d = '<=270m'),
      list(mask = !is.na(interp) &  interp & !is_near, sa = 'Interpolation', d = '>270m'),
      list(mask = !is.na(interp) & !interp,            sa = 'Extrapolation', d = 'all')
    )
    rows <- map_dfr(cells, function(c)
      score_subset(pred, c$mask, res_tag, iter, c$sa, c$d))
    write.csv(rows, out, row.names = FALSE)

    if(i %% 25 == 0 || i == length(model_files))
      message('  ', i, ' / ', length(model_files))
  }
}
