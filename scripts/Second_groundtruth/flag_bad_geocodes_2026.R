## Flag historic presence records whose geocodes the 2026 opportunistic
## ground-truthing revealed to be wrong (e.g. plants not found nearby,
## coordinates falling in clearly unsuitable terrain on revisit).
##
## Source of truth for *which* records are bad:
##   data/GroundTruthing/PublicPresenceRecordsToRemoveBasedOn2026GroundTruth.csv
## Those rows carry a `fid` that references the GPKG feature id in the
## 'known_pres' layer of ~/Documents/Ecolo_groundtruth/ecolo_gt2/ground_truth_median.gpkg
## (the working QGIS project used to review the 2026 field data), NOT any ID
## column inside this repo's modelling data. This script resolves that fid to
## a geometry, then spatially matches it back onto 3m-presence-iter1.gpkg -
## the historic-record rows there carry no OBJECT/ID of their own (both NA),
## so geometry is the only reliable join key.
##
## This does NOT edit 3m-presence-iter1.gpkg or build iter2 - iter2 is built
## once the rest of the 2026 field data is uploaded. This script's only job is
## to produce a durable, spatial flag file that the iter2-assembly step can
## anti-join against, and to make noise: these are records whose accuracy is
## reported in the manuscript (SDManuscript.Rmd - "16 herbarium records were
## removed due to low geolocation quality", Data Acquisition) - this is a
## SEPARATE, additional batch that slipped through that earlier manual
## review, caught only now via field revisits. That count/sentence likely
## needs a follow-up edit once this batch is confirmed; deliberately not
## touched here since it's a scientific-content call, not a data-wiring one.

library(sf)
library(dplyr)

p2proj <- '/home/sagesteppe/Documents/Ecoloradense'
gt_project <- path.expand('~/Documents/Ecolo_groundtruth/ecolo_gt2/ground_truth_median.gpkg')

to_remove <- read.csv(file.path(p2proj, 'data', 'GroundTruthing',
                                 'PublicPresenceRecordsToRemoveBasedOn2026GroundTruth.csv')) |>
  mutate(fid = as.integer(fid))

known_pres <- st_read(gt_project, layer = 'known_pres', fid_column_name = 'fid', quiet = TRUE) |>
  mutate(fid = as.integer(fid))

missing_fids <- setdiff(to_remove$fid, known_pres$fid)
if (length(missing_fids) > 0) {
  stop('fid(s) in PublicPresenceRecordsToRemoveBasedOn2026GroundTruth.csv not found in known_pres: ',
       paste(missing_fids, collapse = ', '))
}

bad_pts <- known_pres[known_pres$fid %in% to_remove$fid, ] |>
  left_join(to_remove, by = 'fid')

iter1 <- st_read(file.path(p2proj, 'data', 'Data4modelling', '3m-presence-iter1.gpkg'), quiet = TRUE)
stopifnot(st_crs(bad_pts) == st_crs(iter1))

nn <- st_nearest_feature(bad_pts, iter1)
dist_m <- as.numeric(st_distance(bad_pts, iter1[nn, ], by_element = TRUE))

## an exact-geometry match is expected (these points were copied, not
## re-digitized) - anything else means iter1 was rebuilt/coordinates moved
## since this CSV was written, and needs a human look rather than a silent
## fuzzy match.
tol_m <- 0.01
if (any(dist_m > tol_m)) {
  bad <- which(dist_m > tol_m)
  stop('No exact-geometry match in 3m-presence-iter1.gpkg for fid(s): ',
       paste(bad_pts$fid[bad], collapse = ', '),
       ' (nearest match ', round(dist_m[bad], 2), 'm away) - investigate before flagging.')
}

flagged <- bad_pts |>
  mutate(
    iter1_row = nn,
    iter1_Type = iter1$Type[nn],
    dist_to_iter1_m = dist_m
  ) |>
  select(fid, date, site, notes, Presenc, iter1_row, iter1_Type, dist_to_iter1_m)

out_path <- file.path(p2proj, 'data', 'GroundTruthing', '2026_flagged_bad_geocodes.gpkg')
st_write(flagged, out_path, append = FALSE, quiet = TRUE)

message(strrep('!', 72))
message('BAD GEOCODES FLAGGED FOR REMOVAL - ', nrow(flagged), ' historic presence record(s)')
message('Confirmed via 2026 field revisits; matched at 0m to rows in 3m-presence-iter1.gpkg.')
message('Written to: ', out_path)
message('These are NOT yet removed from 3m-presence-iter1.gpkg or any iter2 dataset -')
message('exclude them (anti-join on geometry) when 3m-presence-iter2.gpkg is assembled.')
message('')
message('MANUSCRIPT IMPACT: SDManuscript.Rmd Data Acquisition (~line 149) reports')
message('"16 herbarium records were removed due to low geolocation quality" from manual')
message('review - this batch slipped through that review and was only caught by field')
message('revisits. That sentence/count likely needs a follow-up once this is finalized.')
message(strrep('!', 72))
print(st_drop_geometry(flagged))
