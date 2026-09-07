## =============================================================================
## Transect placement for population-boundary verification
## Steps 4-6: dominant edge (OBB) -> cardinal (N/S/E/W) rays across the flanks ->
##            rank by threshold disagreement -> transects as start + bearing + seq
##
## Inputs you bring from steps 1-3:
##   hulls   : list of 3 concave-hull polygons (your threshold rules), any order
##   origins : sf POINTS to cast from (pixel centres inside the shape, or a subsample)
##   occ     : (optional) sf POINTS of occurrences / high-count cells, for the
##             inside anchor in step 6
## =============================================================================
suppressPackageStartupMessages({ 
  library(sf); 
  library(dplyr); 
  library(geosphere); 
  library(terra); 
  library(tidyverse);
  library(leaflet)
  }
)
set.seed(27)

ver <- '1-3arc-Iteration1-PA1:1DO:0-Seed1061'
SITE_CROP_BUFFER_M <- 300   # shared: polygonize step below and the boundaries() step further down

tr_path<- file.path('..','..',  'results', 'threshold_masks', paste0(ver, '-thresholds.tif'))
test_rast <- terra::rast(tr_path)

## from thresholds 
## use only: 'sensitivity' and 'spec_sens' NOT 'equal_spec_sens'
test_rast <- test_rast[[c('spec_sens', 'sensitivity')]]

## and restrict to ground truth areas.
## Cochetopa run: known-good area (2024 "Cochetopa Dome" field data already
## exists there) - scoped to just this Site for this pass.
sample_areas = st_read(file.path('..', '..', 'data', 'hikingTrails', 'GroundTruch_Areas.gpkg')) |>
  #filter(Site == 'Cochetopa') |>
  st_transform(terra::crs(test_rast)) |>
  vect()

## add known presences to be able to tag polygons
pres <- st_read(
  file.path('..', '..', 'data', 'Data4modelling', '3m-presence-iter1.gpkg'),
  quiet = T
) |>
  filter(Presenc == 1) |>
  select(geom)


## polygonize each threshold layer separately (masked to sampleable ground), tagged by
## source layer, then stack. as.polygons() on both layers at once would dissolve by the
## unique COMBINATION of the two layers' values - we need each layer's patches kept
## distinct so identical patches across layers can be deduped below.
##
## Cached to disk: mask()+as.polygons() on the full 1-3arc domain (~330M cells,
## the union of just the 27 Sites' extents is still ~177M) is the actual slow
## part of this script (tens of minutes) - the boundaries()-crash fix further
## down was a different bottleneck. Cropped per-Site here, same reasoning as
## the boundaries() loop, and the result is cached so re-running this script
## while iterating on the transect-placement rules below doesn't repeat it.
patches_cache_path <- file.path('..', '..', 'results', 'GroundTruthSampling',
                                 paste0(ver, '-patches_by_layer_cache.rds'))

if (file.exists(patches_cache_path)) {

  message('Loading cached polygonized patches: ', patches_cache_path)
  test_v <- readRDS(patches_cache_path)

} else {

  patches_by_layer <- list()
  for (s in seq_len(nrow(sample_areas))) {
    site_ext <- terra::ext(terra::buffer(sample_areas[s, ], SITE_CROP_BUFFER_M))
    for (nm in names(test_rast)) {
      lyr_site <- try(terra::crop(test_rast[[nm]], site_ext), silent = TRUE)
      if (inherits(lyr_site, "try-error") || all(is.na(terra::values(lyr_site)))) next
      ## mask against THIS Site only (sample_areas[s,]), not all 27 - nearby
      ## Sites' buffered crop windows can overlap, and masking against the
      ## full sample_areas would let a patch straddling two such windows get
      ## polygonized twice (clipped differently each time by each Site's own
      ## crop extent), breaking the 1:1 patch <-> diff join further down.
      sf_p <- mask(lyr_site, sample_areas[s, ]) |>
        as.polygons() |>
        sf::st_as_sf() |>
        st_cast('POLYGON')
      if (nrow(sf_p) == 0) next
      patches_by_layer[[length(patches_by_layer) + 1]] <- st_sf(layer = nm, geometry = st_geometry(sf_p))
    }
  }

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

  saveRDS(test_v, patches_cache_path)
  message('Cached polygonized patches to: ', patches_cache_path)
}
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


## write both threshold layers out for GUI GIS inspection, one gpkg two layers
threshold_layers_path <- file.path('..', '..', 'results', 'GroundTruthSampling', paste0(ver, '-threshold_layers.gpkg'))
### The ':' characters in ver (e.g. "PA1:1DO:0") make GDAL's driver-guessing
### misread the path as a "driver:connection" string - name the driver explicitly.
st_write(sensitivity,
  threshold_layers_path,
  layer = 'sensitivity',
  delete_dsn = TRUE,
  quiet = TRUE,
  driver = 'GPKG'
)

st_write(spec_sens,
  threshold_layers_path,
  layer = 'spec_sens',
  append = TRUE,
  quiet = TRUE,
  driver = 'GPKG'
)

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

## align the diffs back to their source... this captures diffs from the same patch.
## st_join(), not bind_cols(unlist(st_covered_by(...))): the latter assumes every
## diff polygon covers exactly one test_v patch, which per-Site cropping can
## occasionally violate at a crop seam (a diff sliver touching two adjacent
## patches, or briefly touching none) - st_join() handles 0-or-many matches
## without erroring; ties are broken by keeping the first match. (diagnostic
## only - `diff` isn't consumed by anything past this point.)
diff = sf::st_join(diff, test_v[, c('ID', 'populated')], join = st_covered_by, left = TRUE) |>
  dplyr::rename(src_ID = ID)
n_before <- nrow(diff)
diff = diff |> dplyr::group_by(diff_ID) |> dplyr::slice(1) |> dplyr::ungroup()
if (nrow(diff) < n_before) {
  message(sprintf('%d diff polygon(s) covered by >1 patch (crop-seam artifact); kept first match.',
                   n_before - nrow(diff)))
}
if (any(is.na(diff$src_ID))) {
  message(sprintf('%d diff polygon(s) covered by no patch; src_ID left NA.', sum(is.na(diff$src_ID))))
}


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
transect_headings <- function(tr, utm_crs, decl = DECLINATION) {
  true_b <- mapply(function(x, y, dx, dy, length_m) {
    ends <- rbind(c(x, y), c(x + length_m * dx, y + length_m * dy))
    ll   <- st_coordinates(st_transform(
              st_sfc(st_point(ends[1, ]), st_point(ends[2, ]), crs = utm_crs), 4326))
    bearingRhumb(ll[1, ], ll[2, ])
  }, tr$x, tr$y, tr$dx, tr$dy, tr$length_m)

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
## the outward direction for free from its own 4 cardinal neighbours (N/S/E/W) -
## no curve fitting, and every transect walks a pure compass direction.

core_layer <- if (st_area(spec_sens) > st_area(sensitivity)) "sensitivity" else "spec_sens"
core_full  <- test_rast[[core_layer]]

#core_full <- terra::fillHoles(core_full, nearest = FALSE)

## raster-native from here down, but per-Site: terra::boundaries() on the full
## 1-3arc domain (~330M cells, almost all NA outside the ground-truth areas)
## segfaults - see git history (b0a39fc/840c7cd), where this ran scoped to a
## single Site ('Cochetopa') and, even unscoped, only ever ran at 3arc
## resolution (~90m, ~81x fewer cells than 1-3arc), never the full 1-3arc
## domain. Instead of dropping the other 26 Sites, crop to each Site in turn
## (small +300m buffer so a patch near a Site's edge doesn't get a false
## boundary cell from the crop window itself) and combine the edge points.
## core_bool/core_edge stay small per-Site SpatRasters; nothing here ever
## materializes the full domain.

## direction per edge cell: 4-connected (N/S/E/W) neighbours only, so every
## candidate transect direction is a pure cardinal - no diagonal/tilted rays.
## A cell exposed on two cardinal sides (a corner) yields one candidate per side.
cardinal <- data.frame(
  dr = c(-1,  1,  0, 0),   # N, S, W, E neighbour offsets (row, col)
  dc = c( 0,  0, -1, 1),
  dx = c( 0,  0, -1, 1),   # +col = east (+x)
  dy = c( 1, -1,  0, 0)    # -row = north (+y)
)

edge_rows <- list()
cell_offset <- 0   # keeps `cell` distinguishable across Sites - informational only

for (s in seq_len(nrow(sample_areas))) {

  site_ext  <- terra::ext(terra::buffer(sample_areas[s, ], SITE_CROP_BUFFER_M))
  core_rast <- try(terra::crop(core_full, site_ext), silent = TRUE)
  if (inherits(core_rast, "try-error") || all(is.na(terra::values(core_rast)))) next

  core_bool <- !is.na(core_rast)
  core_edge <- terra::boundaries(core_rast, inner = TRUE, directions = 4)
  edge_pts_v <- terra::as.points(terra::ifel(core_edge == 1, core_edge, NA))
  if (length(edge_pts_v) == 0) next
  edge_cells <- terra::cellFromXY(core_edge, terra::crds(edge_pts_v))

  edge_rows_site <- lapply(edge_cells, function(cel) {
    rc      <- terra::rowColFromCell(core_rast, cel)
    nb_cell <- terra::cellFromRowCol(core_rast, rc[1] + cardinal$dr, rc[2] + cardinal$dc)
    outside <- !core_bool[nb_cell][[1]]     # core_bool[cells] -> 1-col data.frame; unwrap to vector
    outside[is.na(outside)] <- FALSE        # off-grid neighbour: not an outward edge
    if (!any(outside)) return(NULL)         # interior cell, or exposed only diagonally
    xy <- terra::xyFromCell(core_rast, cel)
    data.frame(x = xy[1], y = xy[2], dx = cardinal$dx[outside], dy = cardinal$dy[outside],
               cell = cel + cell_offset)
  })

  edge_rows <- c(edge_rows, edge_rows_site)
  cell_offset <- cell_offset + terra::ncell(core_rast)
}
length(Filter(Negate(is.null), edge_rows))

edge_pts_raw <- st_as_sf(
  do.call(rbind, edge_rows),
  coords = c("x", "y"), crs = st_crs(test_v), remove = FALSE
)

## tag each edge point with its source patch, restricted to populated patches -
## core verification budget only goes to patches with a known occurrence
## within 100m (see `populated` above). Patches without one are handled as a
## separate, flagged "exploratory" stage further down, not silently dropped.
core_patches <- filter(test_v, layer == core_layer, populated)
edge_pts <- st_join(edge_pts_raw, core_patches[, "ID"], join = st_intersects) |>
  filter(!is.na(ID))

## a candidate transect is only usable if it stays on sampleable ground for its
## whole length - a good number of them were straying outside sample_areas
sample_areas_u <- st_union(st_as_sf(sample_areas))

valid_transect <- function(x, y, dx, dy, length_m, step_m = 10) {
  stn <- stations_along(c(x, y), c(dx, dy), length_m, step_m)
  pts <- st_as_sf(stn, coords = c("x", "y"), crs = st_crs(test_v))
  all(lengths(st_intersects(pts, sample_areas_u)) > 0)
}

## the OTHER threshold layer - the larger, looser one that fully contains
## core_layer. Stopping the outward leg at the core edge only clears the
## strict rule; ground just past it can still read "suitable" under the loose
## rule, which is why transects were landing entirely inside threshold.
outer_poly <- if (core_layer == "sensitivity") spec_sens else sensitivity

INWARD_MARGIN    <- 10    # start this far inside the core edge, on confirmed ground -
## dropped from 50m: median spec_sens patch width at 1-3arc is ~6m (99.9% of
## patches under 50m wide), so a 50m inward step almost always overshot clean
## through the patch. 10m still exceeds the median but is far more of these
## patches' actual footprint than 50m ever was.
CLEAR_MARGIN     <- 15    # walk this far past the outer threshold once it's cleared
MAX_TRANSECT_LEN <- 150   # hard cap on total transect length (field logistics)
MAX_OUTWARD <- MAX_TRANSECT_LEN - INWARD_MARGIN - CLEAR_MARGIN  # cap on the disagreement-
## band search so INWARD_MARGIN + out_dist + CLEAR_MARGIN can never exceed MAX_TRANSECT_LEN

## distance along (dx,dy) from (x,y) at which the ray exits `poly` - i.e. how
## far you have to walk from the strict edge to clear the loose threshold too.
## NA if it doesn't clear within max_reach (disagreement band too wide here).
exit_distance <- function(x, y, dx, dy, poly, max_reach = MAX_OUTWARD) {
  ray    <- st_sfc(st_linestring(rbind(c(x, y), c(x + max_reach * dx, y + max_reach * dy))),
                    crs = st_crs(test_v))
  inside <- st_intersection(ray, st_geometry(poly))
  if (length(inside) == 0 || all(st_is_empty(inside))) return(0)
  pts <- st_coordinates(inside)
  out <- max(sqrt((pts[, "X"] - x)^2 + (pts[, "Y"] - y)^2))
  if (out >= max_reach) return(NA_real_)   # never actually exited within the cap
  out
}

## per-patch transect quota: proportional-to-size with a floor, same logic as
## the region allocation in sdm_groundtruth_design_3arc.R (weight on
## area^size_exponent so a few huge patches don't swallow the whole budget,
## floor guarantees every patch still gets sampled, ceiling caps the biggest).
MIN_PER_PATCH  <- 3     # floor: every patch gets at least this many transects
MAX_PER_PATCH  <- 12    # ceiling: no single patch gets more than this many
SIZE_EXPONENT  <- 0.5   # allocate on area^exponent, not raw area

patch_ids  <- unique(edge_pts$ID)
patch_area <- setNames(core_patches$area[match(patch_ids, core_patches$ID)], patch_ids)

N_TRANSECTS_BUDGET <- 5 * length(patch_ids)  # total budget across all patches
## (keeps the old flat-5-per-patch total as the overall budget, just reallocated by size)

budget_for_prop <- max(N_TRANSECTS_BUDGET - MIN_PER_PATCH * length(patch_ids), 0)
weight <- patch_area^SIZE_EXPONENT
prop   <- weight / sum(weight)
n_patch <- MIN_PER_PATCH + round(prop * budget_for_prop)
n_patch <- pmin(n_patch, MAX_PER_PATCH)
n_patch <- setNames(n_patch, patch_ids)

## fallback for patches with no reachable disagreement band: a patch that both
## rules agree is suitable well past MAX_OUTWARD in every direction never
## produces a candidate with a finite exit_distance(), so the main loop below
## keeps 0 transects for it - MIN_PER_PATCH's floor was never actually
## guaranteed. Fill the floor with fixed-length transects confirming
## agreed-suitable ground instead of a disagreement band.
FALLBACK_LEN <- INWARD_MARGIN + CLEAR_MARGIN

## per patch: evaluate every edge (the edge is the core-strict boundary; the
## outward leg is extended per-candidate until it clears the loose threshold
## too, so every kept transect runs confirmed-core -> disagreement -> truly
## outside), keep the ones that pass both the exit and containment checks,
## then take `target` of them, the ones most spread out across this patch -
## not just the first ones drawn - so coverage isn't clumped in one corner.
## Shared across the core and exploratory stages below - only the patch pool
## and its target count differ between them.
place_patch_transects <- function(pid, target, edge_pts_pool) {
  pool <- filter(edge_pts_pool, ID == pid) |>
    st_drop_geometry() |>
    transmute(cell, dx, dy, edge_x = x, edge_y = y)

  pool$out_dist <- mapply(exit_distance, pool$edge_x, pool$edge_y, pool$dx, pool$dy,
                           MoreArgs = list(poly = outer_poly))
  pool$len <- INWARD_MARGIN + pool$out_dist + CLEAR_MARGIN
  pool$sx  <- pool$edge_x - INWARD_MARGIN * pool$dx
  pool$sy  <- pool$edge_y - INWARD_MARGIN * pool$dy

  eligible <- which(!is.na(pool$out_dist))
  eligible <- Filter(function(k)
    valid_transect(pool$sx[k], pool$sy[k], pool$dx[k], pool$dy[k], pool$len[k]),
    eligible)

  ## spread coverage across the patch via greedy sequential selection, not a
  ## static rank: a candidate's raw average distance to EVERY other eligible
  ## point favours whatever sits at the geometric extreme of the patch (the
  ## habitat edge/tail), so the top-k by that score cluster right back
  ## together out there. Instead: seed with the most isolated candidate,
  ## then repeatedly add whichever remaining candidate is, on average,
  ## farthest from the ones already picked - each new pick is chosen relative
  ## to what's already accepted, not to the whole pool.
  if (length(eligible) <= 1 || target <= 1) {
    keep <- head(eligible, target)
  } else {
    coords  <- as.matrix(pool[eligible, c("edge_x", "edge_y")])
    full_d  <- as.matrix(dist(coords))              # indices are into `eligible`
    n_elig  <- length(eligible)

    avg_dist_all <- rowSums(full_d) / (n_elig - 1)
    picked    <- which.max(avg_dist_all)            # 1) seed: most isolated overall
    remaining <- setdiff(seq_len(n_elig), picked)

    while (length(picked) < min(target, n_elig) && length(remaining) > 0) {
      avg_to_picked <- rowMeans(full_d[remaining, picked, drop = FALSE])
      nxt <- remaining[which.max(avg_to_picked)]     # 2/3) farthest from accepted so far
      picked    <- c(picked, nxt)
      remaining <- setdiff(remaining, nxt)
    }
    keep <- eligible[picked]
  }
  lens <- pool$len[keep]
  type <- rep("disagreement", length(keep))

  ## floor guarantee: this patch's disagreement band never resolved (or wasn't
  ## wide enough to hit MIN_PER_PATCH) - fall back to fixed-length transects on
  ## untried edges, checked only for staying on sampleable ground. Floor is
  ## capped at `target` too - a patch trimmed below MIN_PER_PATCH by the
  ## per-Site cap must not get re-inflated back up here.
  floor_n <- min(MIN_PER_PATCH, target)
  if (length(keep) < floor_n) {
    fb_ord <- sample(setdiff(seq_len(nrow(pool)), keep))
    j <- 1
    while (length(keep) < floor_n && j <= length(fb_ord)) {
      cand <- fb_ord[j]
      sx <- pool$edge_x[cand] - INWARD_MARGIN * pool$dx[cand]
      sy <- pool$edge_y[cand] - INWARD_MARGIN * pool$dy[cand]
      if (valid_transect(sx, sy, pool$dx[cand], pool$dy[cand], FALLBACK_LEN)) {
        keep <- c(keep, cand)
        lens <- c(lens, FALLBACK_LEN)
        type <- c(type, "fallback")
      }
      j <- j + 1
    }
  }
  if (length(keep) == 0) return(NULL)
  mutate(pool[keep, ],
         patch_ID = pid, length_m = lens, type = type,
         x = edge_x - INWARD_MARGIN * dx, y = edge_y - INWARD_MARGIN * dy) |>
    select(-edge_x, -edge_y, -out_dist, -len, -sx, -sy)
}

## ---- STAGE 2 setup: flagged (unpopulated) patches --------------------------
## Patches with no known occurrence within 100m get no proportional-to-size
## share of the core budget above - sending crews to purely model-predicted
## ground with no supporting record risks a fishing expedition. Instead: flag
## them, and place just the guaranteed floor on each, kept separate (own
## `stage` tag, own output layer) so it's a conscious add-on, never silently
## blended into the core field sheet.
EXPLORATORY_PER_PATCH <- MIN_PER_PATCH  # flat floor only, no size-proportional share

core_patches_flagged <- filter(test_v, layer == core_layer, !populated)
edge_pts_flagged <- st_join(edge_pts_raw, core_patches_flagged[, "ID"], join = st_intersects) |>
  filter(!is.na(ID))
patch_ids_flagged <- unique(edge_pts_flagged$ID)
n_patch_flagged <- setNames(rep(EXPLORATORY_PER_PATCH, length(patch_ids_flagged)), patch_ids_flagged)

## ---- per-Site cap: no Site gets more than MAX_PER_SITE transects total ----
## (core + exploratory combined), mirroring max_per_region in
## sdm_groundtruth_design_3arc.R. Patch quotas are allocated per-patch above,
## but a field crew visits per-Site, so the ceiling has to apply at that
## level too - tag each patch with the Site polygon it falls in, then trim.
MAX_PER_SITE <- 15

site_v <- st_read(file.path('..', '..', 'data', 'hikingTrails', 'GroundTruch_Areas.gpkg'), quiet = TRUE) |>
  st_transform(st_crs(test_v)) |>
  select(Site)
patches_joined <- st_join(bind_rows(core_patches, core_patches_flagged), site_v,
                           join = st_intersects, largest = TRUE)
patch_site <- setNames(patches_joined$Site, patches_joined$ID)

## trim one unit at a time from whichever patch currently holds the largest
## quota in the over-budget Site - exploratory quotas go first (lowest
## priority), then core - until that Site's total is at/under the cap.
## Excess is just dropped, not redistributed elsewhere (same policy as the
## region cap in sdm_groundtruth_design_3arc.R).
for (s in unique(patch_site)) {
  core_ids <- names(n_patch)[patch_site[names(n_patch)] == s]
  flag_ids <- names(n_patch_flagged)[patch_site[names(n_patch_flagged)] == s]
  over <- sum(n_patch[core_ids], n_patch_flagged[flag_ids]) - MAX_PER_SITE
  while (over > 0 && any(n_patch_flagged[flag_ids] > 0)) {
    k <- flag_ids[which.max(n_patch_flagged[flag_ids])]
    n_patch_flagged[k] <- n_patch_flagged[k] - 1
    over <- over - 1
  }
  while (over > 0 && any(n_patch[core_ids] > 0)) {
    k <- core_ids[which.max(n_patch[core_ids])]
    n_patch[k] <- n_patch[k] - 1
    over <- over - 1
  }
}

tr_core <- bind_rows(lapply(patch_ids, function(pid)
  place_patch_transects(pid, n_patch[as.character(pid)], edge_pts))) |>
  mutate(stage = "core")

tr_flagged <- bind_rows(lapply(patch_ids_flagged, function(pid)
  place_patch_transects(pid, n_patch_flagged[as.character(pid)], edge_pts_flagged))) |>
  mutate(stage = "exploratory")

tr <- bind_rows(tr_core, tr_flagged)   # core rows first: crossing-check below favours core on ties
tr <- transect_headings(tr, utm_crs = st_crs(test_v))

## no two transects may cross - among any crossing pair, drop the later one
tr_lines <- st_sfc(
  lapply(seq_len(nrow(tr)), function(k)
    st_linestring(rbind(c(tr$x[k], tr$y[k]),
                         c(tr$x[k] + tr$length_m[k] * tr$dx[k],
                           tr$y[k] + tr$length_m[k] * tr$dy[k])))),
  crs = st_crs(test_v)
)

crosses <- st_crosses(tr_lines, tr_lines, sparse = FALSE)
pairs   <- which(upper.tri(crosses) & crosses, arr.ind = TRUE)

drop <- logical(nrow(tr))
for (r in seq_len(nrow(pairs))) {
  i <- pairs[r, 1]; j <- pairs[r, 2]
  if (!drop[i]) drop[j] <- TRUE
}
if (any(drop)) message(sum(drop), " transect(s) dropped for crossing another transect")
tr <- tr[!drop, ]

stns <- do.call(rbind, lapply(seq_len(nrow(tr)), function(k)
  cbind(patch_ID = tr$patch_ID[k], transect = k, type = tr$type[k], stage = tr$stage[k],
        stations_along(c(tr$x[k], tr$y[k]), c(tr$dx[k], tr$dy[k]), length_m = tr$length_m[k]))))

# # deliverables: a start-waypoint + bearing table is the field sheet; lines optional
#write.csv(
#  tr[, c("transect","x","y","bearing_true","bearing_mag","heading_note")],
#  "field_sheet.csv"
#)

transects_path <- file.path('..', '..', 'results', 'GroundTruthSampling', paste0(ver, "-transects.gpkg"))
### driver named explicitly - see note above on ver's ':' characters.
st_write(
  st_as_sf(filter(stns, stage == "core"), coords = c("x","y"), crs = 32613),
  transects_path, layer = "core", delete_dsn = TRUE, quiet = TRUE, driver = 'GPKG'
)
if (any(stns$stage == "exploratory")) {
  st_write(
    st_as_sf(filter(stns, stage == "exploratory"), coords = c("x","y"), crs = 32613),
    transects_path, layer = "exploratory", append = TRUE, quiet = TRUE, driver = 'GPKG'
  )
} else {
  message("No flagged (unpopulated) patches produced exploratory transects.")
}

