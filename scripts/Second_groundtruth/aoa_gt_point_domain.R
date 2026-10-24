## Interpolation vs extrapolation of every grid-search model AT the 2026 ground-truth points.
##
## Instead of building a full presence-referenced AOA raster per model (slow, and there are
## 500+ fitted models), predict the AOA only where it is needed - at the ~283 ground-truth
## points - and return the binary in/out verdict per point. Mirrors
## predict_gridsearch_to_points.R (same GT points, same per-resolution predictor extraction,
## same model walk), but:
##   - loads the caret `train` object (modelsTune/, needed for the DI reference cloud and
##     variable weights), not the ranger model used for the probability predictions;
##   - builds a PRESENCE-only trainDI (presenceTrainDI(); pseudo-absences are background, not
##     occurrence data, so the applicability domain is defined by the known niche);
##   - runs CAST::aoa() on the point rows and reports AOA (1 = inside) -> 'Interpolation',
##     0 = outside -> 'Extrapolation', with the raw DI alongside.
## The 270 m near-population subset is just a filter on dist_known_occ downstream - not a
## separate run.
##
## Pairs with the 2026 field metrics (aggregate_2026gt_metrics.R / predict_gridsearch_to_
## points.R): join on (tree, model, point) to split those scores by interpolation vs
## extrapolation. The full presence-referenced SURFACES, when wanted for maps, come from
## aoaSurfaces(reference = 'presence') - the slow path this script avoids.
##
## Usage (working directory scripts/): Rscript Second_groundtruth/aoa_gt_point_domain.R

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(tidyverse)
  library(CAST)
})

p2proj <- '/media/steppe/hdd/EriogonumColoradenseTaxonomy'
p2proc <- file.path(p2proj, 'data', 'spatial', 'processed')
source(file.path(p2proj, 'scripts', 'functions.R'))   # rastReader(), presenceTrainDI()

## 1. Ground-truth points (identical construction to predict_gridsearch_to_points.R) -------

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

## 2. Locate every fitted grid-search model (caret objects in modelsTune/) -----------------

cv_roots <- c(
  'results'                               = 'spatial knndm',
  'results_classicsplit'                  = 'classic',
  'results_populationsplit'               = 'LOPO',
  'results_populationgroupsplit'          = 'leave-pop-group-out',
  'Results_spatialknn_iter1a_high_pa'     = 'spatial knndm, high-PA pool',
  'results_highPA_fieldabs'               = 'spatial knndm, high-PA pool + field abs',
  'results_highPA_populationsplit_fieldabs' = 'LOPO, high-PA pool + field abs',
  'Results_random_iter1a_high_pa'         = 'random, high-PA pool + field abs',
  'Results_spatialmapping_iter1a_high_pa' = 'spatial_mapping (Partitioner), high-PA pool + field abs',
  'Results_LMPO_iter1a_high_pa'           = 'LMPO (Partitioner), high-PA pool + field abs'
)
# caret models live in modelsTune/; exclude the -trainDI / -platt side files
caret_models <- function(root)
  grep('-(trainDI|platt)\\.rds$', invert = TRUE, value = TRUE,
       list.files(file.path(p2proj, root, 'modelsTune'), pattern = '\\.rds$', full.names = TRUE))
cv_roots <- cv_roots[vapply(names(cv_roots), function(r) length(caret_models(r)) > 0, logical(1))]

parse_model_id <- function(f) {
  id <- sub('\\.rds$', '', basename(f))
  tibble(
    path = f, model = id,
    resolution = sub('-Iteration.*$', '', id),
    iteration  = as.integer(sub('.*-Iteration([0-9]+)-PA.*', '\\1', id)),
    pa_ratio   = sub('.*-PA(.+)DO:.*', '\\1', id),
    dist_order = sub('.*DO:([^-]+)-Seed.*', '\\1', id),
    seed       = as.integer(sub('.*-Seed([0-9]+)$', '\\1', id))
  )
}

model_files <- map_dfr(names(cv_roots), function(root)
  bind_cols(map_dfr(caret_models(root), parse_model_id),
            tree = root, cv_structure = unname(cv_roots[root])))

## 3. Point-level predictor values once per resolution (reused by every model) -------------

res_labels <- c('3arc' = '3 arc-second', '1arc' = '1 arc-second', '1-3arc' = '1/3 arc-second')
point_covs_by_res <- map(unique(model_files$resolution), function(res)
  terra::extract(rastReader(paste0('dem_', res), p2proc), vect(gt))) |>
  set_names(unique(model_files$resolution))

## 4. Presence-referenced AOA verdict at each GT point -------------------------------------

classify_model <- function(path, tree, model, resolution, iteration, pa_ratio, dist_order,
                           seed, cv_structure, ...) {
  covs <- point_covs_by_res[[resolution]]
  m    <- readRDS(path)                        # caret train object
  tdi  <- presenceTrainDI(m, LPD = FALSE)       # reference cloud = presences only
  nd   <- covs[, tdi$variables, drop = FALSE]   # tdi$variables handles any Boruta pruning
  a    <- CAST::aoa(nd, trainDI = tdi, verbose = FALSE)
  aoa  <- as.integer(a$AOA)
  tibble(
    tree = tree, model = model, resolution = unname(res_labels[resolution]),
    iteration = iteration, pa_ratio = pa_ratio, dist_order = dist_order, seed = seed,
    cv_structure = cv_structure,
    point_id = seq_len(nrow(gt)), Presence = gt$Presence, dist_known_occ = gt$dist_known_occ,
    DI = as.numeric(a$DI), aoa = aoa,
    domain = ifelse(is.na(aoa), NA_character_,
                    ifelse(aoa == 1L, 'Interpolation', 'Extrapolation')))
}

points_all <- pmap_dfr(model_files, classify_model)

## 5. Per-model summary (how much of the field set each model interpolates vs extrapolates) -

summary_by_model <- points_all |>
  group_by(tree, cv_structure, model, resolution, iteration, pa_ratio, dist_order, seed) |>
  summarise(n_pts = sum(!is.na(aoa)),
            n_interp = sum(aoa == 1L, na.rm = TRUE),
            n_extrap = sum(aoa == 0L, na.rm = TRUE),
            frac_interp = n_interp / n_pts,
            n_interp_270 = sum(aoa == 1L & dist_known_occ <= 270, na.rm = TRUE),
            n_extrap_270 = sum(aoa == 0L & dist_known_occ <= 270, na.rm = TRUE),
            .groups = 'drop') |>
  arrange(desc(frac_interp))

out <- file.path(p2proj, 'results', 'tables')
dir.create(out, showWarnings = FALSE, recursive = TRUE)
write.csv(points_all,       file.path(out, '2026-groundtruth-aoa-domain-points.csv'),  row.names = FALSE)
write.csv(summary_by_model, file.path(out, '2026-groundtruth-aoa-domain-summary.csv'), row.names = FALSE)

message(sprintf('wrote %d point-verdicts across %d models to %s',
                nrow(points_all), nrow(model_files), out))
print(head(summary_by_model, 10) |> as.data.frame())
