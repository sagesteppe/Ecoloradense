## Evaluate the two threshold rules (sensitivity vs. spec_sens) against the
## field-walked transect data.
##
## Core idea: each transect is a straight, distance-ordered walk (`dist_m`
## stations) that starts on confirmed core ground and crosses BOTH threshold
## boundaries on the same line. So instead of scoring each station as a
## right/wrong point-in-polygon hit, find where the plants actually stop
## along the line and compare that distance to where each rule's polygon
## boundary crosses the same line - a signed, metres-scale offset per rule,
## not just an accuracy number.
##
## Data: both the transect design and the threshold surface are the
## `1-3arc-Iteration1-PA1:1DO:0-Seed1061` model run - `transect_placement_
## Seed1061.R` cuts its transects directly from that model's own threshold
## raster, and the field-walked stations in visited_subset.gpkg's
## `visited_transects` layer coincide with that design's station points at
## 0m (exact match, checked across all 338 visited stations).
##
## Result (confirmed against field records, 2026-09-20): both threshold
## rules overshoot the true population edge - `sensitivity` markedly less so
## (bias -2.1m, MAE 9.1m, n=6 usable transects) than `spec_sens` (bias
## -12.7m, MAE 13.5m). The near-total absence outside those 6 transects
## (0/24 core-fallback stations, 0/119 exploratory-stage stations, 26/33
## core-disagreement transects left-censored even at their nominally
## confirmed-core stations) is genuine field absence, not incomplete
## searching - both threshold cuts substantially overpredict actual
## population geography.

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(terra)
})

p2proj <- '/media/steppe/hdd/EriogonumColoradenseTaxonomy'

VER <- '1-3arc-Iteration1-PA1:1DO:0-Seed1061'
TRANSECTS_PATH    <- file.path(p2proj, 'results', 'GroundTruthSampling', paste0(VER, '-transects.gpkg'))
THRESHOLD_TIF      <- file.path(p2proj, 'results', 'threshold_masks',      paste0(VER, '-thresholds.tif'))
SAMPLE_AREAS_PATH  <- file.path(p2proj, 'data', 'hikingTrails', 'GroundTruch_Areas.gpkg')
VISITED_PATH       <- file.path(p2proj, 'data', 'GroundTruthPts', 'visited_subset.gpkg')

## ---- 1. load polygons + station points -------------------------------------
## Polygons are built directly from THRESHOLD_TIF (same polygonize -> dedupe ->
## min-area-filter -> union recipe as transect_placement_Seed1061.R), not read
## from a pre-saved threshold_layers.gpkg, so swapping VER is the only change
## needed to re-point this at a different model run.

build_threshold_polygons <- function(tif_path, sample_areas_path) {
  test_rast <- terra::rast(tif_path)[[c('spec_sens', 'sensitivity')]]
  sample_areas <- st_read(sample_areas_path, quiet = TRUE) |>
    st_transform(terra::crs(test_rast)) |>
    vect()

  patches_by_layer <- lapply(names(test_rast), function(nm) {
    mask(test_rast[[nm]], sample_areas) |>
      as.polygons() |>
      st_as_sf() |>
      st_cast('POLYGON')
  })
  test_v <- bind_rows(setNames(patches_by_layer, names(test_rast)), .id = 'layer')
  test_v <- test_v[!duplicated(st_equals(test_v)), ]

  min_area <- 2 * prod(terra::res(test_rast))
  test_v <- test_v |>
    mutate(area = as.numeric(st_area(geometry))) |>
    filter(area >= min_area)

  list(
    sensitivity = test_v |> filter(layer == 'sensitivity') |> st_union() |> st_as_sf(),
    spec_sens   = test_v |> filter(layer == 'spec_sens')   |> st_union() |> st_as_sf()
  )
}

thresholds  <- build_threshold_polygons(THRESHOLD_TIF, SAMPLE_AREAS_PATH)
sensitivity <- thresholds$sensitivity
spec_sens   <- thresholds$spec_sens

stations <- bind_rows(
  st_read(TRANSECTS_PATH, layer = 'core',        quiet = TRUE),
  st_read(TRANSECTS_PATH, layer = 'exploratory', quiet = TRUE)
)

## ---- 2. reconstruct per-transect geometry from its stations ----------------
## transect_placement_Seed1061.R writes stations but not the line itself -
## direction and length are recoverable from the two extreme stations
## (min/max dist_m) of each transect, so there's no need to re-derive dx/dy
## from scratch.

transect_lines <- stations |>
  st_drop_geometry() |>
  bind_cols(st_coordinates(stations)) |>
  group_by(transect, patch_ID, stage, type) |>
  summarise(
    x0 = X[which.min(dist_m)], y0 = Y[which.min(dist_m)],
    x1 = X[which.max(dist_m)], y1 = Y[which.max(dist_m)],
    length_m = max(dist_m),
    .groups = 'drop'
  ) |>
  mutate(
    dx = (x1 - x0) / length_m,
    dy = (y1 - y0) / length_m
  )

## ---- 3. crossing distance: where a rule's boundary crosses this transect ---
## Mirrors exit_distance() in transect_placement_Seed1061.R. A transect can
## clip a polygon boundary more than once near a corner - take the outermost
## crossing (max distance from the start), consistent with "the edge is the
## outermost point this rule calls suitable along this line".
crossing_distance <- function(x0, y0, dx, dy, length_m, poly, crs) {
  ray <- st_sfc(st_linestring(rbind(c(x0, y0), c(x0 + length_m * dx, y0 + length_m * dy))), crs = crs)
  hit <- suppressWarnings(st_intersection(ray, st_boundary(st_geometry(poly))))
  if (length(hit) == 0 || all(st_is_empty(hit))) return(NA_real_)
  pts <- st_coordinates(hit)
  max(sqrt((pts[, 'X'] - x0)^2 + (pts[, 'Y'] - y0)^2))
}

transect_lines <- transect_lines |>
  rowwise() |>
  mutate(
    cross_sensitivity = crossing_distance(x0, y0, dx, dy, length_m, sensitivity, st_crs(stations)),
    cross_spec_sens   = crossing_distance(x0, y0, dx, dy, length_m, spec_sens,   st_crs(stations))
  ) |>
  ungroup()

## ---- 4. observed edge distance from field presence/absence -----------------
## Walks stations outward (increasing dist_m) and returns the midpoint of the
## LAST present -> absent transition - "last" because a transect can be
## patchy near the true edge (present/absent/present), and the design intent
## is the population's outermost extent, not the first gap encountered.
## Two failure modes get flagged rather than silently dropped:
##   - "right_censored": every station present - the true edge is beyond
##     this transect's 150m cap, offset is a lower bound only.
##   - "left_censored": every station absent, including the nominally
##     confirmed-core start - no transition to anchor an edge estimate on,
##     so it can't contribute to the offset average (see the module-level
##     Result note: this was confirmed the majority case in practice).
observed_edge_distance <- function(dist_m, presence) {
  ord <- order(dist_m)
  d <- dist_m[ord]; p <- presence[ord]
  if (all(p == 1)) return(list(edge = NA_real_, flag = 'right_censored'))
  if (all(p == 0)) return(list(edge = NA_real_, flag = 'left_censored'))
  transitions <- which(p[-length(p)] == 1 & p[-1] == 0)
  if (length(transitions) == 0) return(list(edge = NA_real_, flag = 'no_present_to_absent_step'))
  k <- max(transitions)
  flag <- if (length(transitions) > 1) 'non_monotonic' else 'ok'
  list(edge = mean(c(d[k], d[k + 1])), flag = flag)
}

## ---- 5. field data: visited stations, joined onto the design's own IDs -----
## visited_subset.gpkg's `visited_transects` layer carries patch_ID + dist_m
## but not `transect` (a patch can hold several transects, so patch_ID alone
## doesn't identify one) - recovered here via nearest-feature match against
## `stations`, which is exact (0m, and dist_m agrees for all 338 visited
## stations) since both trace back to the same design.
load_field_data <- function(visited_path, stations) {
  visited <- st_read(visited_path, layer = 'visited_transects', quiet = TRUE) |>
    filter(visited) |>
    st_transform(st_crs(stations))

  nn <- st_nearest_feature(visited, stations)
  stopifnot('Visited station did not match a design station at 0m' =
              all(as.numeric(st_distance(visited, stations[nn, ], by_element = TRUE)) == 0))

  visited$transect <- stations$transect[nn]

  ## drop visited's own patch_ID/source/visited - transect_lines' patch_ID
  ## (joined in downstream via `transect`) is the one to trust, and keeping
  ## both around just collides on the join.
  visited |>
    st_drop_geometry() |>
    mutate(presence = as.integer(if_any(c(Prsnc_M, Prsnc_J, Prsnc_S), ~ replace_na(. > 0, FALSE)))) |>
    select(transect, dist_m, presence)
}

## ---- 6. compute offsets -----------------------------------------------------
evaluate_transect_boundary <- function(visited_path = VISITED_PATH) {
  field <- load_field_data(visited_path, stations)

  by_transect <- field |>
    group_by(transect) |>
    summarise(edge_info = list(observed_edge_distance(dist_m, presence)), .groups = 'drop') |>
    mutate(
      observed_edge = map_dbl(edge_info, 'edge'),
      flag          = map_chr(edge_info, 'flag')
    ) |>
    select(-edge_info) |>
    left_join(transect_lines, by = 'transect')

  ## type == 'disagreement' only: fallback-type transects are placed on
  ## agreed-suitable ground by construction and never cross either boundary,
  ## so cross_sensitivity/cross_spec_sens are structurally NA for them - they
  ## answer "is this agreed-suitable ground actually occupied", not "which
  ## rule's edge is closer", and get scored separately below instead of
  ## diluting n with guaranteed NAs.
  offsets <- by_transect |>
    filter(flag == 'ok', type == 'disagreement') |>
    mutate(
      offset_sensitivity = observed_edge - cross_sensitivity,
      offset_spec_sens   = observed_edge - cross_spec_sens
    )

  ## Sign convention: positive offset = the rule's boundary sits INSIDE the
  ## true edge (the population extends farther than that rule says); negative
  ## = the rule overshoots past where plants actually stop.
  summary_by_rule <- offsets |>
    filter(stage == 'core') |>
    pivot_longer(c(offset_sensitivity, offset_spec_sens), names_to = 'rule', values_to = 'offset_m') |>
    mutate(rule = sub('^offset_', '', rule)) |>
    group_by(rule) |>
    summarise(
      n = sum(!is.na(offset_m)),
      bias_m = mean(offset_m, na.rm = TRUE),
      mae_m  = mean(abs(offset_m), na.rm = TRUE),
      .groups = 'drop'
    )

  ## Fallback-type transects: no boundary comparison possible, but still a
  ## useful check - do plants actually occur on ground both rules already
  ## agree is suitable? Confirmed 0% here (see the module-level Result note)
  ## - undercuts trusting either rule's *interior*, independent of where its
  ## edge falls.
  fallback_occupancy <- by_transect |>
    filter(type == 'fallback') |>
    left_join(field |> group_by(transect) |> summarise(any_present = any(presence == 1), .groups = 'drop'),
              by = 'transect') |>
    group_by(stage) |>
    summarise(n = n(), occupied_rate = mean(any_present), .groups = 'drop')

  ## Secondary sanity check: cheap per-station point-in-polygon accuracy,
  ## reusing the same joined data - not the primary metric (throws away the
  ## distance-ordering the transect design paid for) but a quick cross-check
  ## that the offset numbers above aren't an artefact of the crossing-distance
  ## logic.
  station_hits <- stations |>
    inner_join(field, by = c('transect', 'dist_m'))
  station_hits$in_sensitivity <- lengths(st_intersects(station_hits, sensitivity)) > 0
  station_hits$in_spec_sens   <- lengths(st_intersects(station_hits, spec_sens))   > 0
  station_hits <- station_hits |>
    st_drop_geometry() |>
    pivot_longer(c(in_sensitivity, in_spec_sens), names_to = 'rule', values_to = 'predicted') |>
    mutate(rule = sub('^in_', '', rule), predicted = as.integer(predicted))

  accuracy_by_rule <- station_hits |>
    group_by(rule) |>
    summarise(
      n = n(),
      accuracy = mean(predicted == presence),
      sensitivity = sum(predicted == 1 & presence == 1) / sum(presence == 1),
      specificity = sum(predicted == 0 & presence == 0) / sum(presence == 0),
      .groups = 'drop'
    )

  ## Exploratory-stage transects have no disagreement band to calibrate
  ## against (their patches failed the `populated` filter in
  ## transect_placement_Seed1061.R) - report detection rate only, kept
  ## separate from the core-stage boundary comparison above rather than
  ## folded into it.
  exploratory_detection <- field |>
    left_join(select(transect_lines, transect, stage, patch_ID), by = 'transect') |>
    filter(stage == 'exploratory') |>
    group_by(patch_ID) |>
    summarise(any_present = any(presence == 1), .groups = 'drop')

  list(
    by_transect = by_transect,
    offsets = offsets,
    summary_by_rule = summary_by_rule,
    fallback_occupancy = fallback_occupancy,
    accuracy_by_rule = accuracy_by_rule,
    exploratory_detection = exploratory_detection,
    censoring = count(by_transect, flag)
  )
}

result <- evaluate_transect_boundary()
print(result$summary_by_rule)
print(result$fallback_occupancy)
print(result$accuracy_by_rule)
print(result$censoring)
