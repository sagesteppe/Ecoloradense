## Build 3m-presence-iter2.gpkg from the collated presences file.
##
## The 2026 iteration ("iter2") of presences is defined as:
##   iter1 field-verified presences (346)
## + 2026 opportunistic ground-truth points (96)
## + 2026 transect points (13)
## - the 7 flagged bad geocodes identified in flag_bad_geocodes_2026.R
##
## data/Data4modelling/all_presences_collated.gpkg already satisfies this: the
## collate step dropped the flagged rows (nearest flagged point is >180 m from
## any collated point). This script verifies that invariant and writes the
## canonical iter2 presences file.
##
## Produces: data/Data4modelling/3m-presence-iter2.gpkg

suppressMessages({
  library(sf)
  library(dplyr)
})

p2proj <- rprojroot::find_root(rprojroot::is_git_root)

coll_path <- file.path(p2proj, 'data', 'Data4modelling', 'all_presences_collated.gpkg')
flg_path  <- file.path(p2proj, 'data', 'GroundTruthing', '2026_flagged_bad_geocodes.gpkg')
out_path  <- file.path(p2proj, 'data', 'Data4modelling', '3m-presence-iter2.gpkg')

coll <- st_read(coll_path, quiet = TRUE)
flg  <- st_read(flg_path,  quiet = TRUE)

min_dist_m <- st_distance(st_transform(flg, 3857),
                          st_transform(coll, 3857)) |>
  as.numeric() |>
  matrix(nrow = nrow(flg)) |>
  apply(1, min)

if (any(min_dist_m < 1)) {
  stop('Flagged bad geocodes still present in collated file; aborting.')
}

iter2 <- coll |>
  mutate(point_id = paste0('iter2_', row_number())) |>
  select(point_id, source, date, site, geom)

st_write(iter2, out_path, append = FALSE, quiet = TRUE)

message(sprintf('wrote %d iter2 presences to %s', nrow(iter2), out_path))
print(table(iter2$source))
