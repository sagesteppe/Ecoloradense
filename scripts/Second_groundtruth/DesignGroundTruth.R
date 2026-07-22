library(terra)
library(sf)
library(tidyverse)
library(clhs)
source('functions.R')

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')

# Generates ~n_total ground-truth sampling points that traverse three (or four)
# correlated axes of the occupancy logistic regression:
#   1. Habitat suitability  (above sensitivity threshold -> 0.99)
#   2. Euclidean distance to nearest known population (capped at euclid_cap_m)
#   3. Least-cost distance to nearest known population (anisotropic preferred)
#   4. Model uncertainty — SE from infinitesimal jackknife (optional; include
#      when se_prediction = TRUE has been run for this version)
#
# Points are selected via conditioned Latin Hypercube Sampling (cLHS), run
# separately on each of two geographic clusters to guarantee both range areas
# (CB area and Cocheotopa Dome) are represented.
#
# Run after CostDistances.R has written outputs for the target version.

# =========================================================================
# Parameters — edit here
# =========================================================================

ver          <- '1-3arc-Iteration1-PA1:2.7DO:0'   # model version to use
cost_thresh  <- 'spec_sens'   # threshold type for cost distance surfaces
n_total      <- 100           # total ground-truth points
euclid_cap_m <- 2400          # ~1.5 miles; hard cap on euclidean distance
n_candidates <- 5000          # candidate pool drawn from suitable cells
clhs_iter    <- 10000         # simulated annealing iterations for cLHS

# Trail feasibility
# trail_buffer_m: hard buffer — candidates outside this distance from any
#   trail are dropped before cLHS. Applied only to the non-Cocheotopa cluster;
#   Cocheotopa Dome is open terrain and does not require trail access.
# cocheotopa_bbox: xmin, xmax, ymin, ymax in the same CRS as the suitability
#   raster (UTM 32613). The k-means cluster whose centroid falls inside this
#   box is treated as Cocheotopa and receives no trail constraint.
trail_path       <- '../data/hikingTrails/trails.shp'
trail_buffer_m   <- 150
cocheotopa_bbox  <- c(xmin = , xmax = , ymin = , ymax = )  # fill in UTM 32613 coords

# =========================================================================
# File paths
# =========================================================================

p2res  <- '../results'
p2data <- '../data/Data4modelling'

suit_path  <- file.path(p2res, 'suitability_maps', paste0(ver, '-Pr.tif'))
aoa_path   <- file.path(p2res, 'suitability_maps', paste0(ver, '-AOA.tif'))
thresh_path <- file.path(p2res, 'threshold_masks',  paste0(ver, '-thresholds.tif'))
se_path    <- file.path(p2res, 'suitability_maps', paste0(ver, '-SE.tif'))

aniso_path <- file.path(p2res, 'cost_distances',
                        paste0(ver, '-', cost_thresh, '-costDistAniso.tif'))
iso_path   <- file.path(p2res, 'cost_distances',
                        paste0(ver, '-', cost_thresh, '-costDist.tif'))

stopifnot(
  'Suitability raster not found'  = file.exists(suit_path),
  'AOA raster not found'          = file.exists(aoa_path),
  'Threshold raster not found'    = file.exists(thresh_path)
)

if (!file.exists(aniso_path) && !file.exists(iso_path)) {
  stop(
    'No cost distance file found for version: ', ver, ' / threshold: ', cost_thresh,
    '\nRun CostDistances.R first.'
  )
}

# =========================================================================
# Load rasters
# =========================================================================

suit_r   <- terra::rast(suit_path)
aoa_r    <- terra::rast(aoa_path)
thresh_r <- terra::rast(thresh_path)

# Use sensitivity layer as the candidate pool boundary (most inclusive threshold)
sens_r <- thresh_r[['sensitivity']]

cost_r <- if (file.exists(aniso_path)) {
  message('Using anisotropic cost distance.')
  terra::rast(aniso_path)
} else {
  message('Anisotropic cost distance not yet available; using isotropic.')
  terra::rast(iso_path)
}

use_se <- file.exists(se_path)
if (use_se) {
  se_r <- terra::rast(se_path)
  message('SE layer found — model uncertainty included as 4th cLHS axis.')
} else {
  message('No SE layer found; proceeding on 3 axes. Re-run with se_prediction = TRUE to add uncertainty axis.')
}

# =========================================================================
# Known presences — source of euclidean distances
# 3m-presence-iter1.gpkg used throughout as the canonical presence record
# (same source as patchAttributes uses for occupied patch identification)
# =========================================================================

known_pops <- sf::st_read(
  file.path(p2data, '3m-presence-iter1.gpkg'), quiet = TRUE
) |>
  dplyr::filter(Presenc == 1) |>
  sf::st_transform(terra::crs(suit_r))

# =========================================================================
# Candidate pool — random sample from AOA-masked suitable cells
# =========================================================================

# Mask to AOA then to sensitivity threshold so all candidates are both within
# the model's inference area and predicted suitable at the most inclusive level
domain <- terra::mask(suit_r, aoa_r,   maskvalues = 0) |>
          terra::mask(sens_r,           maskvalues = NA)

candidates <- terra::spatSample(
  domain, size = n_candidates,
  method = 'random', na.rm = TRUE, xy = TRUE, as.df = TRUE
) |>
  setNames(c('x', 'y', 'suit')) |>
  sf::st_as_sf(coords = c('x', 'y'), crs = terra::crs(suit_r), remove = FALSE)

message(nrow(candidates), ' raw candidates drawn.')

# =========================================================================
# Extract predictor values at candidates
# =========================================================================

candidates$cost_dist <- terra::extract(cost_r, candidates)[[2]]

near_idx <- sf::st_nearest_feature(candidates, known_pops)
candidates$euclid_dist <- as.numeric(
  sf::st_distance(candidates, known_pops[near_idx, ], by_element = TRUE)
)

if (use_se) {
  candidates$se <- terra::extract(se_r, candidates)[[2]]
}

# =========================================================================
# Trail feasibility — compute distance to nearest trail for each candidate
# Applied as a hard per-cluster filter inside the cLHS loop, not here.
# =========================================================================

use_trails <- !is.null(trail_path) && file.exists(trail_path)

if (use_trails) {
  trails <- sf::st_read(trail_path, quiet = TRUE) |>
    sf::st_transform(sf::st_crs(candidates))
  near_trail_idx  <- sf::st_nearest_feature(candidates, trails)
  candidates$trail_dist <- as.numeric(
    sf::st_distance(candidates, trails[near_trail_idx, ], by_element = TRUE)
  )
  message('Trail distances computed (buffer: ', trail_buffer_m, ' m).')
} else {
  message('No trail layer supplied — feasibility filter disabled.')
}

# =========================================================================
# Filter: euclidean cap + drop NA extractions
# (trail buffer applied per-cluster below; Cocheotopa is exempt)
# =========================================================================

keep_cols <- c('suit', 'euclid_dist', 'cost_dist', if (use_se) 'se')

candidates_filt <- candidates |>
  sf::st_drop_geometry() |>
  dplyr::filter(
    euclid_dist <= euclid_cap_m,
    dplyr::if_all(dplyr::all_of(keep_cols), \(v) !is.na(v))
  )

message(nrow(candidates_filt), ' candidates after filtering to euclid_cap and removing NAs.')

if (nrow(candidates_filt) < n_total * 3) {
  warning(
    'Candidate pool is small relative to n_total. Consider increasing ',
    'n_candidates or relaxing euclid_cap_m.'
  )
}

# =========================================================================
# Split into two geographic clusters — ensures both range areas are sampled
# (CB area and Cocheotopa Dome are ~50 km apart; a single-domain cLHS
# risks concentrating in the larger/denser cluster)
# =========================================================================

km <- kmeans(candidates_filt[, c('x', 'y')], centers = 2, nstart = 25)
candidates_filt$cluster <- km$cluster

# Identify Cocheotopa cluster — whichever centroid falls inside cocheotopa_bbox
cocheto_cluster <- which(
  km$centers[, 'x'] >= cocheotopa_bbox[['xmin']] &
  km$centers[, 'x'] <= cocheotopa_bbox[['xmax']] &
  km$centers[, 'y'] >= cocheotopa_bbox[['ymin']] &
  km$centers[, 'y'] <= cocheotopa_bbox[['ymax']]
)
if (length(cocheto_cluster) != 1L) {
  stop(
    'cocheotopa_bbox matched ', length(cocheto_cluster),
    ' cluster centroid(s) — adjust the bounding box so exactly one centroid falls inside it.'
  )
}
message('Cocheotopa cluster: ', cocheto_cluster, ' (no trail constraint).')

cluster_sizes <- table(candidates_filt$cluster)
# Allocate proportionally to cluster size, minimum 20 points per cluster
alloc <- round(n_total * prop.table(cluster_sizes))
alloc <- pmax(alloc, 20L)
alloc <- round(alloc * (n_total / sum(alloc)))  # rescale to n_total

message('Cluster allocation: ', paste(names(alloc), alloc, sep = '=', collapse = ', '))

# =========================================================================
# cLHS — run independently per cluster
# x and y are included so geographic spread is baked into the sample;
# they are treated as additional axes to balance, not just coordinates
# =========================================================================

clhs_axes <- c(keep_cols, 'x', 'y')

selected <- vector('list', 2)
for (k in seq_along(alloc)) {

  sub <- candidates_filt[candidates_filt$cluster == k, ]

  if (use_trails && k != cocheto_cluster) {
    sub <- dplyr::filter(sub, trail_dist <= trail_buffer_m)
    message(
      'Cluster ', k, ': trail buffer (', trail_buffer_m, ' m) applied — ',
      nrow(sub), ' candidates remain.'
    )
  } else {
    message('Cluster ', k, ' (Cocheotopa): no trail constraint — ', nrow(sub), ' candidates.')
  }

  if (nrow(sub) < alloc[k]) {
    stop(
      'Cluster ', k, ' has only ', nrow(sub), ' candidates after trail filter ',
      'but needs ', alloc[k], ' points. Increase n_candidates or trail_buffer_m.'
    )
  }

  message('Cluster ', k, ': running cLHS (n=', alloc[k], ', pool=', nrow(sub), ')...')

  idx <- clhs::clhs(
    x        = sub[, clhs_axes],
    size     = alloc[k],
    iter     = clhs_iter,
    progress = FALSE,
    simple   = TRUE
  )
  selected[[k]] <- sub[idx, ]
}

# =========================================================================
# Assemble output
# =========================================================================

ground_truth <- dplyr::bind_rows(selected) |>
  dplyr::mutate(
    pt_id     = paste0('GT-', sprintf('%03d', dplyr::row_number())),
    cost_type = ifelse(file.exists(aniso_path), 'anisotropic', 'isotropic'),
    version   = ver,
    .before   = 1
  ) |>
  sf::st_as_sf(coords = c('x', 'y'), crs = terra::crs(suit_r))

dir.create(file.path(p2res, 'GroundTruthSampling'), showWarnings = FALSE)

out_path <- file.path(
  p2res, 'GroundTruthSampling',
  paste0(ver, '-', cost_thresh, '-groundTruth.gpkg')
)
sf::st_write(ground_truth, out_path, delete_dsn = TRUE, quiet = TRUE)
message(nrow(ground_truth), ' ground-truth points written to:\n  ', out_path)

# =========================================================================
# Quick coverage check — verify axes are well-spread across selected points
# =========================================================================

summary_cols <- c(keep_cols, if (use_trails) 'trail_dist')

ground_truth |>
  sf::st_drop_geometry() |>
  dplyr::summarise(
    dplyr::across(
      dplyr::all_of(summary_cols),
      list(min = min, med = median, max = max),
      .names = '{.col}_{.fn}'
    )
  ) |>
  tidyr::pivot_longer(everything()) |>
  print(n = Inf)
