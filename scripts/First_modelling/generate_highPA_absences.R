## Generate candidate background absences to extend the iteration-1 absence pool far
## enough to fit presence:absence ratios up to 1:50 (see highPA_grid.R).
##
## The existing pool (iter1-pa.gpkg, 1071 gDist background points, plus each
## resolution's own field absences) caps the achievable ratio at ~1:4 - 1:6 depending on
## resolution; distOrder_PAratio_simulator() silently caps any higher request. Reweighted
## holdout scoring (summarize_iter_calibration_reweighted.R) put the best-calibrated
## ratio at the top of the tested range at realistic (1-3%) landscape prevalences, so the
## grid needs to extend well past it.
##
## Same method as the original pool (GenerateAbsences.Rmd, "Iteration 1 absences"):
## sdm::background(method = 'gDist') on the 1 arc-second stack around the 30m presences.
## Sized so the pool supports 1:52 at the resolution with the most presences (1/3
## arc-second), leaving headroom for the manual review to remove points in plausibly
## suitable habitat and still reach 1:50.
##
## Workflow:
##   1. Run this script -> data/Data4modelling/iter1-pa-highPA-candidates.gpkg
##   2. In QGIS, delete candidates that plausibly fall in suitable habitat (as was done
##      for the original pool) and save the result as
##      data/Data4modelling/iter1-pa-highPA-reviewed.gpkg (same layer/columns).
##   3. highPA_grid.R reads the reviewed file, records which candidate IDs were removed
##      (data/Data4modelling/iter1-pa-highPA-removed.csv, for the manuscript), and
##      refuses to run if too few remain for 1:50.
##
## Run with the working directory set to scripts/ - e.g.
## `Rscript First_modelling/generate_highPA_absences.R`.

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
})
source('functions.R')
set.seed(52)

target_ratio <- 52
p2proc <- file.path('..', 'data', 'spatial', 'processed')
p2dat  <- file.path('..', 'data', 'Data4modelling')

## 1. How many new points are needed ----------------------------------------------------
existing_pool <- st_read(file.path(p2dat, 'iter1-pa.gpkg'), quiet = TRUE)

pool_need <- lapply(c('90m', '30m', '3m'), function(r){
  p <- st_read(file.path(p2dat, paste0(r, '-presence-iter1.gpkg')), quiet = TRUE)
  p <- p[st_is(p, 'POINT'), ]
  data.frame(resolution = r, n_presence = sum(p$Presenc == 1),
             n_field_abs = sum(p$Presenc == 0),
             needed = target_ratio * sum(p$Presenc == 1) - sum(p$Presenc == 0) - nrow(existing_pool))
}) |> bind_rows()
print(pool_need)
n_new <- max(pool_need$needed)
message(sprintf('Generating %d new candidate absences (1:%d at the finest resolution).', n_new, target_ratio))

## 2. Sample, same method as the original pool ------------------------------------------
arc1 <- rastReader('dem_1arc', p2proc)
m30 <- st_read(file.path(p2dat, '30m-presence-iter1.gpkg'), quiet = TRUE) |>
  filter(Presenc == 1, st_is(geom, 'POINT'))
all_known <- bind_rows(
  # drop the empty GEOMETRYCOLLECTION rows - GEOS treats distance to an empty geometry
  # as 0, so st_is_within_distance() would flag every candidate
  st_read(file.path(p2dat, '3m-presence-iter1.gpkg'), quiet = TRUE) |>
    filter(st_is(geom, 'POINT')) |> select(),
  existing_pool |> select()
)

# oversample: some draws are dropped below for landing on/near existing points
cand <- sdm::background(arc1, n = ceiling(n_new * 1.1), method = 'gDist', sp = vect(m30)) |>
  select(x, y) |>
  st_as_sf(coords = c('x', 'y'), crs = 32613)

# drop candidates within one 1 arc-second cell of an existing presence/absence, then
# thin back to the number needed
near_known <- lengths(st_is_within_distance(cand, all_known, dist = 30)) > 0

#ggplot() + 
#  geom_sf(data = cand) +
#  geom_sf(data = all_known, col = 'red')

cand <- cand[!near_known, ]
if(nrow(cand) < n_new) warning('Only ', nrow(cand), ' candidates after dropping near-duplicates; ',
                               n_new, ' were wanted - rerun with a larger oversample.')
cand <- cand[sample(seq_len(nrow(cand)), min(n_new, nrow(cand))), ]

## 3. Annotate for review and write --------------------------------------------------------
# IDs continue after the original pool's so the two can be combined without collisions;
# dist_presence_m is there to style/filter by in QGIS (points nearest known presences
# are the ones most worth a look).
presences <- st_read(file.path(p2dat, '3m-presence-iter1.gpkg'), quiet = TRUE) |>
  filter(Presenc == 1, st_is(geom, 'POINT'))
cand <- cand |>
  mutate(
    Occurrence = 0,
    ID = max(existing_pool$ID) + seq_len(n()),
    dist_presence_m = round(as.numeric(st_distance(cand, presences[st_nearest_feature(cand, presences), ],
                                                   by_element = TRUE)))
  )

out <- file.path(p2dat, 'iter1-pa-highPA-candidates.gpkg')
st_write(cand, out, append = FALSE, quiet = TRUE)
message(sprintf('Wrote %d candidates to %s (IDs %.0f-%.0f). Review in QGIS and save as iter1-pa-highPA-reviewed.gpkg.',
                nrow(cand), out, min(cand$ID), max(cand$ID)))
print(summary(cand$dist_presence_m))
