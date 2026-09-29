## Fixed presence:absence ratio grid, 1:1 - 1:50, on the extended absence pool.
##
## The adaptive search (adaptive_PAratio_search(), Modelling.Rmd) only reached ~1:4, the
## cap of the original absence pool, and reweighted holdout scoring
## (summarize_iter_calibration_reweighted.R) found the best-calibrated ratio at the top
## of that range at realistic landscape prevalences. This grid pushes past it:
##   1:1 - 1:4 (anchors, refit on the new pool so they're comparable),
##   1:5 - 1:10 by 1, 1:12 - 1:20 by 2, 1:25 - 1:50 by 5.
##
## Three split designs, mirroring the original search's three results trees:
##   spatial - n_splits independent class-stratified knndm splits per resolution
##             (spatial_class_split()), each cached (results_highPA/tables/
##             {res}-split{k}.rds; split 1 is {res}-split.rds) and each run through the
##             whole ratio grid -> results_highPA/
##   classic - modeller()'s random caret::createDataPartition() split, redrawn for every
##             fit from its seed, so the holdout spans all populations
##             -> results_highPA_classicsplit/
##   lopo    - leave-one-population-out (population_lopo_search()) at only a few ratios
##             (lopo_ratios below): every population with >= 3 presences takes a turn
##             as the whole test set -> results_highPA_populationsplit/
## spatial and classic use the same seeds at every ratio (paired design). modeller()
## skips any fit whose model is already on disk, so an interrupted run resumes.
##
## Same modeller() pipeline as every other production fit (Boruta, knndm CV, caret
## tuning, class weights), predict_surface = FALSE - holdout evaluation only. Each tree is
## laid out like results/, so summarize_iter_calibration_reweighted.R <tree> and
## predict_gridsearch_to_points.R can read it.
##
## Needs data/Data4modelling/iter1-pa-highPA-reviewed.gpkg - the output of
## generate_highPA_absences.R after manual review in QGIS.
##
## Usage (working directory scripts/): resolution, then optional mode (default spatial):
##   Rscript First_modelling/highPA_grid.R 90             # 3 arc-second, spatial split
##   Rscript First_modelling/highPA_grid.R 90 classic
##   Rscript First_modelling/highPA_grid.R 90 lopo
##   (30 = 1 arc-second, 10 = 1/3 arc-second - the slowest)
## Optional 3rd argument: number of spatial splits (default 5; ignored for classic/lopo).
## Optional 4th argument: seeds per ratio per split (default 3). For classic every seed is
##   its own random split, so e.g. `90 classic 1 15` gives 15 shuffles per ratio.

suppressPackageStartupMessages({
  library(sf)
  library(tidyverse)
  library(terra)
  library(caret)
  library(ranger)
})
source('functions.R')

args <- commandArgs(trailingOnly = TRUE)
resolution <- as.numeric(args[1])
stopifnot('first argument must be a resolution: 90, 30 or 10' = resolution %in% c(90, 30, 10))
mode <- if(length(args) >= 2) args[2] else 'spatial'
stopifnot('second argument must be spatial, classic or lopo' = mode %in% c('spatial', 'classic', 'lopo'))
# default 1: repeated knndm draws came out near-identical (split 2's test set shared 85%
# of split 1's at 3 arc-second), so extra spatial splits add cost, not independent draws
n_splits <- if(length(args) >= 3) as.integer(args[3]) else 1
n_seeds <- if(length(args) >= 4) as.integer(args[4]) else 3
if(mode != 'spatial') n_splits <- 1 # classic redraws its split every fit; lopo has its own folds

ratios <- c(1:4, 5:10, seq(12, 20, by = 2), seq(25, 50, by = 5))
lopo_ratios <- c(4, 20) # current production ratio vs. well past it
# Split k's knndm draw is seeded 5000 + resolution + 1000 * (k - 1), and its model seeds are
# 5000 + 100 * (k - 1) + 1..n_seeds: modeller()'s file stem carries the seed but not the
# split, so seeds must be unique across splits or split k's fits would reload split 1's.
# Split 1 (split seed 5000 + resolution, model seeds 5001-5003) is the original single-
# split run, so its fits are reused as-is.
split_seed_of <- function(k) 5000 + resolution + 1000 * (k - 1)
seeds_of <- function(k) 5000 + 100 * (k - 1) + seq_len(n_seeds)
# Keep every manually placed field absence (the absences stored in each resolution's
# presence file, median ~160-355 m from a presence) in every fit, and subsample only the
# background absences (median ~30 km away) to reach each ratio. Without this a
# landscape-wide subsample dilutes them away - at 1:4 on the extended pool only ~12 of the
# 247 3 arc-second field absences survive vs ~1/4 of all absences on the original pool - so
# raising the ratio would change both the implied prevalence and how much near-population
# absence the model sees. FALSE reproduces the original (diluted) runs, kept in the
# un-suffixed results_highPA* trees for comparison.
keep_field_absences <- TRUE
root_sfx <- if(keep_field_absences) '_fieldabs' else ''
results_root <- file.path(PROJ_ROOT, paste0(c(spatial = 'results_highPA',
                                              classic = 'results_highPA_classicsplit',
                                              lopo    = 'results_highPA_populationsplit')[[mode]], root_sfx))
p2proc <- file.path('..', 'data', 'spatial', 'processed')
p2dat  <- file.path('..', 'data', 'Data4modelling')
iteration <- 1

for(d in c('models', 'modelsTune', 'test_data', 'tables', 'evaluations', 'suitability_maps'))
  dir.create(file.path(results_root, d), showWarnings = FALSE, recursive = TRUE)

## 1. Absence pool: original + manually reviewed extension --------------------------------
candidates <- st_read(file.path(p2dat, 'iter1-pa-highPA-candidates.gpkg'), quiet = TRUE)
reviewed   <- st_read(file.path(p2dat, 'iter1-pa-highPA-reviewed.gpkg'), quiet = TRUE)
stopifnot('reviewed file contains IDs not in the candidates - wrong file?' =
            all(reviewed$ID %in% candidates$ID))

# record what the review removed - the provenance the manuscript's methods need
removed <- candidates |>
  filter(!ID %in% reviewed$ID) |>
  st_drop_geometry() |>
  select(ID, dist_presence_m)
write.csv(removed, file.path(p2dat, 'iter1-pa-highPA-removed.csv'), row.names = FALSE)
message(sprintf('Manual review removed %d of %d candidates (%.1f%%).',
                nrow(removed), nrow(candidates), 100 * nrow(removed) / nrow(candidates)))

abs <- bind_rows(
  st_read(file.path(p2dat, 'iter1-pa.gpkg'), quiet = TRUE) |> select(Occurrence, ID),
  reviewed |> select(Occurrence, ID)
)

## 2. Modelling data, as in Modelling.Rmd ---------------------------------------------------
pres_file <- c('90' = '90m', '30' = '30m', '10' = '3m')[[as.character(resolution)]]
x <- st_read(file.path(p2dat, paste0(pres_file, '-presence-iter1.gpkg')), quiet = TRUE) |>
  rename(Occurrence = Presenc) |>
  st_as_sf() |>
  mutate(field_abs = Occurrence == 0) # the presence file's own absences are the field ones
x <- bind_rows(x, mutate(abs, field_abs = FALSE))
x <- filter(x, st_is(x, 'POINT')) |>
  select(Occurrence, field_abs) |>
  mutate(.id = seq_len(n()))
fixed_col <- if(keep_field_absences) 'field_abs' else NULL
message(sprintf('%d field absences%s.', sum(x$field_abs),
                if(keep_field_absences) ' kept in every fit' else ' subsampled with the background'))

max_ratio <- sum(x$Occurrence == 0) / sum(x$Occurrence == 1)
message(sprintf('%s: %d presences, %d absences -> max achievable ratio 1:%.1f.',
                res_string(resolution), sum(x$Occurrence == 1), sum(x$Occurrence == 0), max_ratio))
if(max_ratio < max(ratios))
  stop('The reviewed pool supports only 1:', round(max_ratio, 1), ' at this resolution; ',
       max(ratios), ' was requested. Generate more candidates, or lower the grid.')

areas <- st_read(file.path('..', 'data', 'GroundTruthPts', 'ground_truth_median-VISITED.gpkg'),
                 layer = 'areas', quiet = TRUE) |>
  st_transform(st_crs(x))

## 3. Spatial splits: n_splits independent knndm draws, each cached -------------------------
# The train/test draw has been the single largest source of holdout variance for this
# rare, spatially clustered presence class - more than PA ratio or seed - so the grid is
# repeated over several draws rather than trusting one.
get_spatial_split <- function(k){
  sfx <- if(k == 1) '' else k # split 1 keeps the original single-split file name
  split_path <- file.path(results_root, 'tables', paste0(res_string(resolution), '-split', sfx, '.rds'))
  # reuse the un-suffixed tree's split if it exists, so the field-absence runs differ from
  # the diluted runs only in absence sampling (same points, same .id order)
  base_split <- file.path(PROJ_ROOT, 'results_highPA', 'tables', basename(split_path))
  if(!file.exists(split_path) && root_sfx != '' && file.exists(base_split))
    file.copy(base_split, split_path)
  if(file.exists(split_path)){
    split <- readRDS(split_path)
    stopifnot('cached split was built on different data - delete it to rebuild' =
                identical(sort(c(split$train_id, split$test_id)), x$.id))
  } else {
    rast_dat <- rastReader(paste0('dem_', res_string(resolution)), p2proc)
    set.seed(split_seed_of(k)) # spatial_class_split() doesn't seed itself
    split <- spatial_class_split(x, rast_dat, train_split = 0.8, test_tolerance = 0.1)
    saveRDS(split, split_path)
  }

  # Which populations' presences landed in test. A single knndm draw can hold out an
  # entire population (the iteration-1 1/3 arc-second split put all of Cochetopa in
  # test and none in train), which turns the holdout into a leave-that-population-out
  # test - worth knowing before reading anything into these holdout scores.
  prs <- x[x$Occurrence == 1, ]
  prs$Site <- areas$Site[st_nearest_feature(prs, areas)]
  prs$set <- factor(ifelse(prs$.id %in% split$test_id, 'test', 'train'), levels = c('test', 'train'))
  split_by_site <- as.data.frame.matrix(table(prs$Site, prs$set))
  split_by_site$pct_test <- round(100 * split_by_site$test / rowSums(split_by_site[, c('test', 'train')]))
  write.csv(split_by_site, file.path(results_root, 'tables',
                                     paste0(res_string(resolution), '-split', sfx, '-by-site.csv')))
  message('Split ', k, ' - presences in train/test by population:')
  print(split_by_site)
  all_test <- rownames(split_by_site)[split_by_site$train == 0 & split_by_site$test >= 5]
  if(length(all_test) > 0)
    warning('Split ', k, ': population(s) held out entirely in test: ', paste(all_test, collapse = ', '),
            ' - the holdout is effectively leave-population-out for them.', immediate. = TRUE)
  split
}

## 4. Fit ---------------------------------------------------------------------------------
log_path <- file.path(results_root, 'tables', paste0(res_string(resolution), '-highPA-', mode, '-log.csv'))
log_row <- function(row){
  write.table(row, log_path, sep = ',', row.names = FALSE,
              col.names = !file.exists(log_path), append = file.exists(log_path))
}

if(mode %in% c('spatial', 'classic')){
  # split-major order: each split's whole ratio grid finishes before the next split
  # starts, so partial results are always complete grids on some number of splits
  for(k in seq_len(n_splits)){
    # classic: NULL -> modeller() draws its own random split per fit
    split <- if(mode == 'spatial') get_spatial_split(k) else NULL
    for(r in ratios){
      for(s in seeds_of(k)){
        t0 <- Sys.time()
        fit <- tryCatch(
          distOrder_PAratio_simulator(
            x = x, distOrder = 0, PAratio = r, resolution = resolution, seed = s,
            predict_surface = FALSE, se_prediction = FALSE, split = split, fixed_col = fixed_col,
            results_root = results_root, iteration = iteration, train_split = 0.8, p2proc = p2proc),
          error = function(e){ message('FAILED 1:', r, ' seed ', s, ': ', conditionMessage(e)); NULL })
        mins <- as.numeric(difftime(Sys.time(), t0, units = 'mins'))
        row <- data.frame(
          resolution = res_string(resolution), mode = mode,
          split = if(mode == 'spatial') k else NA, requested_ratio = r, seed = s,
          fitted_ratio = if(is.null(fit)) NA else fit$PAratio,
          pr_auc = if(is.null(fit)) NA else fit$metrics$estimate[fit$metrics$metric == 'pr_auc'],
          minutes = round(mins, 1), finished = format(Sys.time()))
        log_row(row)
        message(sprintf('split %s, 1:%s seed %d done in %.1f min (PR-AUC %.3f)', row$split, r, s, mins, row$pr_auc))
      }
    }
  }
} else {
  # One fold per population at each lopo ratio. population_lopo_search() seeds fold i as
  # base_seed + i; the ratio is in modeller()'s file stem, so folds at different ratios
  # don't collide. Its per-population table (with each fold's seed recoverable as
  # 6000 + row) is written to results_highPA_populationsplit/tables/.
  for(r in lopo_ratios){
    t0 <- Sys.time()
    lopo <- population_lopo_search(
      x = select(x, -.id), resolution = resolution, distOrder = 0, iteration = iteration,
      p2proc = p2proc, PAratio = r, areas = areas, base_seed = 6000, results_root = results_root,
      # with 18k absences, a landscape-wide subsample at 1:4 leaves small zones with no
      # test absences; subsample training absences only, test on each zone's full set
      ratio_train_only = TRUE, fixed_col = fixed_col)
    mins <- as.numeric(difftime(Sys.time(), t0, units = 'mins'))
    log_row(cbind(resolution = res_string(resolution), mode = mode, requested_ratio = r,
                  seed = 6000 + seq_len(nrow(lopo$per_population)), lopo$per_population,
                  minutes_all_folds = round(mins, 1), finished = format(Sys.time())))
    message(sprintf('LOPO 1:%s done: %d folds in %.1f min.', r, nrow(lopo$per_population), mins))
  }
}
