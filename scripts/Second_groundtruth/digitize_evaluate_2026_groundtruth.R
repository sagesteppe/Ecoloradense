## Digitize the raw 2026 opportunistic ground-truthing points and evaluate
## the fitted suitability surfaces (3, 1, and 1/3 arc-second) against them.
## Any resolution missing its -Pr.tif under results/suitability_maps/ is
## skipped rather than erroring, since not all resolutions are fit yet.

library(sf)
library(terra)
library(tidyverse)
library(yardstick)

p2proj <- '/home/sagesteppe/Documents/Ecoloradense'

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

evaluate_surface <- function(f, within_m = NULL) {
  res_tag <- sub('-Iteration.*$', '', basename(f))
  r <- rast(f)
  pts <- if (is.null(within_m)) gt else filter(gt, dist_known_occ <= within_m)
  pred <- terra::extract(r, vect(pts))[, 2]
  keep <- !is.na(pred)
  if (sum(keep) == 0 || length(unique(pts$Presence[keep])) < 2) {
    return(tibble(resolution = unname(res_labels[res_tag]), n = sum(keep), pr_auc = NA_real_))
  }
  df <- tibble(truth = factor(pts$Presence[keep], levels = c(0, 1)), pred = pred[keep])
  tibble(resolution = unname(res_labels[res_tag]), n = nrow(df),
         pr_auc = yardstick::pr_auc(df, truth, pred, event_level = 'second')$.estimate)
}

if (length(suit_rasters) == 0) {
  message('No suitability rasters found under results/suitability_maps/ - nothing to evaluate yet.')
  eval_all <- eval_270 <- tibble(resolution = character(), n = integer(), pr_auc = double())
} else {
  eval_all <- map_dfr(suit_rasters, evaluate_surface)
  eval_270 <- map_dfr(suit_rasters, evaluate_surface, within_m = 270)
}

dir.create(file.path(p2proj, 'results', 'tables'), showWarnings = FALSE, recursive = TRUE)
write.csv(eval_all, file.path(p2proj, 'results', 'tables', '2026-groundtruth-evaluation-all.csv'), row.names = FALSE)
write.csv(eval_270, file.path(p2proj, 'results', 'tables', '2026-groundtruth-evaluation-270m.csv'), row.names = FALSE)

eval_270
