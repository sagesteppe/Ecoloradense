## =============================================================================
## Transect placement for population-boundary verification
## Steps 4-6: dominant edge (OBB) -> tilted rays across the flanks ->
##            rank by threshold disagreement -> transects as start + bearing + seq
##
## All geometry is PLANAR (UTM): locally exact at 100 m and lets st_intersection
## run on the hulls without spherical-geometry distortion. geosphere is used only
## at the end to translate the finished planar direction into a compass HEADING a
## person reads and dials - rhumb (constant-heading) bearing, then declination.
##
## Inputs you bring from steps 1-3:
##   hulls   : list of 3 concave-hull polygons (your threshold rules), any order
##   origins : sf POINTS to cast from (pixel centres inside the shape, or a subsample)
##   occ     : (optional) sf POINTS of occurrences / high-count cells, for the
##             inside anchor in step 6
## Requires sf >= 1.0-0 (GEOS >= 3.9) for st_minimum_rotated_rectangle().
## =============================================================================
suppressPackageStartupMessages({ 
  library(sf); 
  library(dplyr); 
  library(geosphere); 
  library(terra); 
  library(tidyverse);
  library(tmap)}
)

setwd('~/Documents/Ecoloradense/scripts')
tr_path<- file.path('..', 'results', 'threshold_masks', '3arc-Iteration1-PA1_3.6DO_0-thresholds.tif')
test_rast <- terra::rast(tr_path)

## from thresholds 
## use only: 'sensitivity' and 'spec_sens' NOT 'equal_spec_sens'
test_rast <- test_rast[[c('spec_sens', 'sensitivity')]]

## and restrict to ground truth areas. 
sample_areas = st_read(file.path('..', 'data', 'hikingTrails', 'GroundTruch_Areas.gpkg')) |>
  st_transform(terra::crs(test_rast)) |>
  vect()

## add known presences to be able to tag polygons
pres <- st_read(
  file.path('..', 'data', 'Data4modelling', '3m-presence-iter1.gpkg'),
  quiet = T
) |>
  filter(Presenc == 1) |>
  select(geom)


## polygonize each threshold layer separately (masked to sampleable ground), tagged by
## source layer, then stack. as.polygons() on both layers at once would dissolve by the
## unique COMBINATION of the two layers' values - we need each layer's patches kept
## distinct so identical patches across layers can be deduped below.
patches_by_layer <- lapply(names(test_rast), function(nm) {
  sf_p <- mask(test_rast[[nm]], sample_areas) |>
    as.polygons() |>
    sf::st_as_sf() |>
    st_cast('POLYGON')
  st_sf(layer = nm, geometry = st_geometry(sf_p))
})

patches_by_layer <- setNames(patches_by_layer, names(test_rast))  # lookup used below

test_v <- do.call(rbind, patches_by_layer) |>
  mutate(ID = row_number())

## remove patches that are identical across the two layers: if spec_sens and sensitivity
## agree exactly on a patch's footprint, it's the same ground - test it once.
nrow(test_v)
test_v <- test_v[!duplicated(st_equals(test_v)), ]
nrow(test_v)

## subset to polygons >= 2 pixels of THIS raster. A flat m^2 number is meaningless once
## sites get pooled across products at different resolutions (1arc vs 3arc) - deriving
## it from res() keeps "2 pixels" the same concept regardless of source grid.
min_area <- 2 * prod(terra::res(test_rast))
###  13339.58 = 2 pixels of 3arc DEM resolution. 
test_v <- test_v |>
  mutate(area = as.numeric(st_area(geometry))) |>
  filter(area >= min_area)
nrow(test_v) 


## identify the difference between patches... 
test_v |>
  group_by(layer) |>
  summarize(area = sum(area))


## areas that overlap
sensitivity = test_v |>
  filter(layer == 'sensitivity') |>
  summarize(geometry = st_union(geometry))

spec_sens = test_v |>
  filter(layer == 'spec_sens') |>
  summarize(geometry = st_union(geometry))

ggplot() + 
  geom_sf(data = sensitivity) + 
  geom_sf(data = spec_sens, fill = 'red')

## obtain only areas where the two methods diverge in 
# classifying habitat
if(st_area(spec_sens) > st_area(sensitivity)){
  diff = st_difference(spec_sens, sensitivity)
} else{
  diff = st_difference(sensitivity, spec_sens)
}
## split back into pieces
diff = st_cast(diff, 'POLYGON')
diff = mutate(diff, diff_ID = row_number())

## tag sites with known occurrences within 100m. 
test_v$populated = lengths(st_intersects(test_v, st_buffer(pres, 100)))>0

## align the diffs back to their source... this captures diffs from the same patch 
diff = bind_cols(
  diff, 
  test_v[
  unlist(st_covered_by(diff, test_v)),c('ID', 'populated')
] |>
  rename(src_ID = ID) |>
  st_drop_geometry()

)


## identify the first

DECLINATION <- 8.4   # degrees EAST, 2026 - magnetic N is 8.4* E of true N.

## ---- start point: step inside the crossing band by inner_margin (confirmed ground)
transect_start <- function(x, y, dx, dy, mid, inner_margin = 20) {
  s0 <- max(mid - inner_margin, 0)
  c(x + s0 * dx, y + s0 * dy)                     # planar start (UTM)
}

## ---- stations: planar, just a seq along the unit direction from the start -----
stations_along <- function(start, dir, length_m = 100, step_m = 10) {
  s <- seq(0, length_m, by = step_m)
  data.frame(dist_m = s, x = start[1] + s * dir[1], y = start[2] + s * dir[2])
}

## ---- compass heading: planar dir -> true rhumb bearing -> magnetic + note ------
## geosphere only here. bearingRhumb = constant compass heading (not great-circle).
## East declination: magnetic = true − decl  ("east is least"). Stores BOTH; the
## true bearing is the source of record, magnetic is derived and dated.
transect_headings <- function(tr, utm_crs, decl = DECLINATION, length_m = 100) {
  true_b <- mapply(function(x, y, dx, dy) {
    ends <- rbind(c(x, y), c(x + length_m * dx, y + length_m * dy))
    ll   <- st_coordinates(st_transform(
              st_sfc(st_point(ends[1, ]), st_point(ends[2, ]), crs = utm_crs), 4326))
    bearingRhumb(ll[1, ], ll[2, ])
  }, tr$x, tr$y, tr$dx, tr$dy)

  tr$bearing_true <- round(true_b %% 360, 1)                 # source of record
  tr$declination  <- decl
  tr$decl_year    <- format(Sys.Date(), "%Y")
  tr$bearing_mag  <- round((true_b - decl) %% 360, 1)        # dial THIS on the compass
  tr$heading_note <- sprintf(
    "walk %05.1f\u00b0 MAG  (%05.1f\u00b0 true, decl +%.1f\u00b0E %s)",
    tr$bearing_mag, tr$bearing_true, decl, tr$decl_year)
  tr
}

## ---- final approach: raster-native edges, no curve-fitting --------------------
## The loess/vector-boundary work above assumed the two threshold layers CROSS;
## they don't (see the diff/ scratch above - `sensitivity` nests entirely inside
## `spec_sens`). So "the edge" is just: cells of the smaller, higher-confidence
## layer (agreed by BOTH rules) that touch a cell outside it. A raster grid gives
## the outward direction for free from its own 8 neighbours - no curve fitting.

core_layer <- if (st_area(spec_sens) > st_area(sensitivity)) "sensitivity" else "spec_sens"
core_rast  <- test_rast[[core_layer]]
core_bool  <- !is.na(as.vector(values(core_rast)))

## edge cells: inside the core layer, with >=1 of the 8 neighbours outside it
cand_cells <- which(core_bool)
adj        <- terra::adjacent(core_rast, cells = cand_cells, directions = 8)
is_edge    <- apply(adj, 1, function(nb) any(!core_bool[nb], na.rm = TRUE))
edge_cells <- cand_cells[is_edge]
length(edge_cells)

## direction per edge cell: unit vector from the core toward its non-core
## neighbours, straight from the grid offsets ("easy directionality" - no smoothing)
offsets <- expand.grid(dr = -1:1, dc = -1:1)
offsets <- offsets[!(offsets$dr == 0 & offsets$dc == 0), ]

edge_dir <- t(vapply(edge_cells, function(cel) {
  rc      <- terra::rowColFromCell(core_rast, cel)
  nb_cell <- terra::cellFromRowCol(core_rast, rc[1] + offsets$dr, rc[2] + offsets$dc)
  outside <- !core_bool[nb_cell]
  outside[is.na(outside)] <- FALSE
  if (!any(outside)) return(c(NA_real_, NA_real_))       # symmetric neighbourhood, no clear way out
  v <- colSums(cbind(offsets$dc, -offsets$dr)[outside, , drop = FALSE])  # +col = east (+x), +row = south (-y)
  v / sqrt(sum(v^2))
}, numeric(2)))

edge_pts <- st_as_sf(
  data.frame(terra::xyFromCell(core_rast, edge_cells),
             dx = edge_dir[, 1], dy = edge_dir[, 2], cell = edge_cells),
  coords = c("x", "y"), crs = st_crs(test_v), remove = FALSE
) |>
  filter(!is.na(dx))   # drop the rare degenerate cell with no clear outward side

## tag each edge point with its source patch, restricted to populated patches
core_patches <- filter(test_v, layer == core_layer, populated)
edge_pts <- st_join(edge_pts, core_patches[, "ID"], join = st_intersects) |>
  filter(!is.na(ID))

## a candidate transect is only usable if it stays on sampleable ground for its
## whole 100 m - a good number of them were straying outside sample_areas
sample_areas_u <- st_union(st_as_sf(sample_areas))

valid_transect <- function(x, y, dx, dy, length_m = 100, step_m = 10) {
  stn <- stations_along(c(x, y), c(dx, dy), length_m, step_m)
  pts <- st_as_sf(stn, coords = c("x", "y"), crs = st_crs(test_v))
  all(lengths(st_intersects(pts, sample_areas_u)) > 0)
}

## per patch: draw untried edges (the edge is the 50 m mark - start is stepped
## 50 m back into the core) until 5 pass the containment check above, capped at
## max_iter draws total so a patch with mostly out-of-bounds edges gives up
## rather than looping forever
max_iter <- 10

tr_list <- lapply(unique(edge_pts$ID), function(pid) {
  pool <- filter(edge_pts, ID == pid) |>
    st_drop_geometry() |>
    transmute(cell, dx, dy, x = x - 50 * dx, y = y - 50 * dy)
  ord  <- sample(nrow(pool))   # draw order, without replacement
  keep <- integer(0)
  i <- 1
  while (length(keep) < 5 && i <= min(max_iter, length(ord))) {
    cand <- ord[i]
    if (valid_transect(pool$x[cand], pool$y[cand], pool$dx[cand], pool$dy[cand]))
      keep <- c(keep, cand)
    i <- i + 1
  }
  if (length(keep) == 0) return(NULL)
  mutate(pool[keep, ], patch_ID = pid)
})

tr <- bind_rows(tr_list)
tr <- transect_headings(tr, utm_crs = st_crs(test_v))

stns <- do.call(rbind, lapply(seq_len(nrow(tr)), function(k)
  cbind(patch_ID = tr$patch_ID[k], transect = k,
        stations_along(c(tr$x[k], tr$y[k]), c(tr$dx[k], tr$dy[k])))))

# # deliverables: a start-waypoint + bearing table is the field sheet; lines optional
#write.csv(
#  tr[, c("transect","x","y","bearing_true","bearing_mag","heading_note")], 
#  "field_sheet.csv"
#)

st_write(
  st_as_sf(stns, coords = c("x","y"), crs = 32613), 
  "transects.gpkg",
  append = F
)