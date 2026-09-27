## Digitize the raw 2026 opportunistic ground-truthing points and evaluate
## the fitted suitability surfaces (3, 1, and 1/3 arc-second) against them.
## Any resolution missing its -Pr.tif under results/suitability_maps/ is
## skipped rather than erroring, since not all resolutions are fit yet.

library(sf)
library(terra)
library(tidyverse)
library(yardstick)

p2proj <- Find(dir.exists, c('/media/steppe/hdd/EriogonumColoradenseTaxonomy', '~/Documents/Ecoloradense'))
source(file.path(p2proj, 'scripts', 'functions.R')) # presenceScores()

## 1. Digitize -----------------------------------------------------------

gt_raw <- read.csv(file.path(p2proj, 'data', 'GroundTruthing', '2026_groundtruth_points.csv')) |>
  rename(Plot = `Plot..if.tissue.`) |>
  mutate(Date = as.Date(Date, format = '%m/%d/%y'),
         Presence = as.integer(Plants > 0))

gt <- gt_raw |>
  st_as_sf(coords = c('Long', 'Lat'), crs = 4326) |>
  st_transform(32613)

st_write(gt, file.path(p2proj, 'data', 'GroundTruthing', '2026_groundtruth_points.gpkg'),
          append = FALSE, quiet = TRUE)

## 1b. Combine with visited random ground-truth points --------------------
## `visited_subset.gpkg` (built in Visited_Random_Points.qmd) holds the
## random base_pts/os_pts that field crews have since visited. A point is
## Present if any of its three visit columns (Prsnc_M/J/S) recorded a
## nonzero count.

visited_gt <- st_read(file.path(p2proj, 'data', 'GroundTruthPts', 'visited_subset.gpkg'),
                       layer = 'visited_points', quiet = TRUE) |>
  st_transform(st_crs(gt)) |>
  mutate(Presence = as.integer(if_any(c(Prsnc_M, Prsnc_J, Prsnc_S), ~ replace_na(. > 0, FALSE))))

## visited_gt's geometry column is named `geom` (from the gpkg) while gt's is
## `geometry` (from st_as_sf) - bind_rows() does not reconcile differing sfc
## column names, so rebuild visited_gt on a `geometry` column before binding.
visited_gt <- st_sf(Presence = visited_gt$Presence, geometry = st_geometry(visited_gt))

gt <- bind_rows(
  gt |> select(Presence),
  visited_gt
)

## 2. Evaluate against available suitability surfaces ---------------------

known_occ <- st_read(file.path(p2proj, 'data', 'collections', 'occurrences_coloradense', 'occurrences.shp'),
                      quiet = TRUE) |>
  st_transform(32613)

gt <- gt |>
  mutate(dist_known_occ = as.numeric(st_distance(gt, known_occ)[cbind(seq_len(nrow(gt)), st_nearest_feature(gt, known_occ))]))

suit_rasters <- list.files(file.path(p2proj, 'results', 'suitability_maps'),
                            pattern = '-Pr\\.tif$', full.names = TRUE)

res_labels <- c('3arc' = '3 arc-second', '1arc' = '1 arc-second', '1-3arc' = '1/3 arc-second')

## PR-AUC plus Brier / Cox calibration / decision-curve scores (presenceScores()); the
## full per-model decision curves go to their own tables.
evaluate_surface <- function(f, within_m = NULL) {
  model <- sub('-Pr\\.tif$', '', basename(f))
  res_tag <- sub('-Iteration.*$', '', model)
  iteration <- as.integer(sub('.*-Iteration([0-9]+)-.*', '\\1', model))
  base <- tibble(model = model, resolution = unname(res_labels[res_tag]), iteration = iteration)
  r <- rast(f)
  pts <- if (is.null(within_m)) gt else filter(gt, dist_known_occ <= within_m)
  pred <- terra::extract(r, vect(pts))[, 2]
  keep <- !is.na(pred) & !is.na(pts$Presence)
  if (sum(keep) == 0 || length(unique(pts$Presence[keep])) < 2) {
    return(list(metrics = bind_cols(base, n = sum(keep), pr_auc = NA_real_), dca = NULL))
  }
  scores <- presenceScores(truth = pts$Presence[keep], prob = pred[keep])
  list(metrics = bind_cols(base, n = sum(keep), scores$wide),
       dca = bind_cols(base, scores$dca))
}

if (length(suit_rasters) == 0) {
  message('No suitability rasters found under results/suitability_maps/ - nothing to evaluate yet.')
  res_all <- res_270 <- list(list(metrics = tibble(model = character(), resolution = character(),
                                                    iteration = integer(), n = integer(), pr_auc = double()),
                                  dca = NULL))
} else {
  res_all <- map(suit_rasters, evaluate_surface)
  res_270 <- map(suit_rasters, evaluate_surface, within_m = 270)
}
eval_all <- map_dfr(res_all, 'metrics')
eval_270 <- map_dfr(res_270, 'metrics')

dir.create(file.path(p2proj, 'results', 'tables'), showWarnings = FALSE, recursive = TRUE)
write.csv(eval_all, file.path(p2proj, 'results', 'tables', '2026-groundtruth-evaluation-all.csv'), row.names = FALSE)
write.csv(eval_270, file.path(p2proj, 'results', 'tables', '2026-groundtruth-evaluation-270m.csv'), row.names = FALSE)
write.csv(map_dfr(res_all, 'dca'), file.path(p2proj, 'results', 'tables', '2026-groundtruth-dca-curve-all.csv'), row.names = FALSE)
write.csv(map_dfr(res_270, 'dca'), file.path(p2proj, 'results', 'tables', '2026-groundtruth-dca-curve-270m.csv'), row.names = FALSE)

eval_all
eval_270
