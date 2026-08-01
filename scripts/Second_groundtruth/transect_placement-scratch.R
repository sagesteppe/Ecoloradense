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
  library(tidyverse)
})

#setwd('./scripts/Second_groundtruth')
tr_path<- file.path('..', '..', 'results', 'threshold_masks', '1-3arc-Iteration1-PA1:2.7DO:0-thresholds.tif')
test_rast <- terra::rast(tr_path)

## from thresholds 
## use only: 'sensitivity' and 'spec_sens' NOT 'equal_spec_sens'
test_rast <- test_rast[[c('spec_sens', 'sensitivity')]]

## and restrict to ground truth areas. 
sample_areas = st_read(file.path('..', '..', 'data', 'hikingTrails', 'GroundTruch_Areas.gpkg')) |>
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

diff_pop = filter(diff, populated == 'TRUE')
table(diff_pop$src_ID)

## view the diffs on some major sites. 
bounds = st_bbox(filter(test_v, ID == 164))
ggplot() + 
  geom_sf(data = filter(test_v, ID == 164)) + 
  geom_sf(
    data = filter(diff, src_ID == 164), 
    fill = 'blue',
    alpha = 0.5) +
 # geom_sf(data = priority_targets, fill = 'pink') + 
  geom_sf(data = pres) +
  coord_sf(
    xlim = c(bounds['xmin'], bounds['xmax']), 
    ylim = c(bounds['ymin'], bounds['ymax'])
  ) + 
  theme_minimal() 

bounds = st_bbox(filter(test_v, ID == 17))
ggplot() + 
  geom_sf(data = filter(test_v, ID == 17)) + 
  geom_sf(
    data = filter(diff, src_ID == 17), 
    fill = 'blue',
    alpha = 0.5)+ 
  geom_sf(data = pres) +
  geom_sf(data = priority_targets, fill = 'pink') + 
  coord_sf(
    xlim = c(bounds['xmin'], bounds['xmax']), 
    ylim = c(bounds['ymin'], bounds['ymax'])
  ) +
  theme_minimal()

bounds = st_bbox(filter(test_v, ID == 49))
ggplot() + 
  geom_sf(data = filter(test_v, ID == 49)) + 
  geom_sf(
    data = filter(diff, src_ID == 49), 
    fill = 'blue',
    alpha = 0.5 )+ 
  geom_sf(data = priority_targets, fill = 'pink') + 
  geom_sf(data = pres) +
  coord_sf(
    xlim = c(bounds['xmin'], bounds['xmax']), 
    ylim = c(bounds['ymin'], bounds['ymax'])
  ) +
  theme_minimal()

## identify the edges most adjacent to the known occurrences 
## ((use this approach to accomodate occurrences outside of the thresholds))


## identify points near the diffs -- highest priority sampling...
## but edge effects may vary.
priority_targets = diff_pop[ lengths(st_is_within_distance(diff_pop, pres, 100))>0,]

ggplot() + 
  geom_sf(data = priority_targets)
###########################################################################