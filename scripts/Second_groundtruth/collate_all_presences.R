## Collates every known presence record - the original herbarium/iNat occurrence
## compilation, the iteration-1 field-verified presences, and both 2026 ground-truthing
## presence sources - into a single gpkg. Meant as a reference layer for hand-drawing a
## new population/"areas" polygon dataset: the grouping that will define spatial-block CV
## folds for the final models (permuting on POPULATION, not individual site, since sites
## cluster tightly within a population - see the cv_structure note in
## Second_groundtruth/predict_gridsearch_to_points.R).

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(purrr)
})

p2proj <- '/media/steppe/hdd/EriogonumColoradenseTaxonomy'

## rep()s source/date/site to nrow(geom) rather than relying on st_sf() recycling, which
## errors instead of recycling a length-1 attribute against a 0-row geometry set (the
## visited-points source below turns out to have zero presences).
make_pts <- function(source, date, site, geom) {
  n <- length(geom)
  st_sf(source = rep(source, length.out = n), date = rep(date, length.out = n),
        site = rep(site, length.out = n), geometry = geom)
}

## 1. Iteration-1 field-verified presences (through 2024), minus records flagged bad by
## 2026 ground-truthing. PublicPresenceRecordsToRemoveBasedOn2026GroundTruth.csv's `fid`
## references a feature id in an *external* QGIS project's gpkg, NOT iter1's own OBJECT
## column (those numbers coincidentally overlap but point at different rows entirely -
## see Second_groundtruth/flag_bad_geocodes_2026.R's header) - its already-resolved
## output (2026_flagged_bad_geocodes.gpkg, matched to iter1 by geometry there) is the
## reliable way to identify these rows here.
iter1_all <- st_read(file.path(p2proj, 'data', 'Data4modelling', '3m-presence-iter1.gpkg'), quiet = TRUE)
bad_geocodes <- st_read(file.path(p2proj, 'data', 'GroundTruthing', '2026_flagged_bad_geocodes.gpkg'), quiet = TRUE)

nn <- st_nearest_feature(bad_geocodes, iter1_all)
dist_m <- as.numeric(st_distance(bad_geocodes, iter1_all[nn, ], by_element = TRUE))
stopifnot('2026_flagged_bad_geocodes.gpkg no longer matches 3m-presence-iter1.gpkg by exact geometry' =
            all(dist_m < 0.01))

iter1 <- iter1_all[-nn, ] |> dplyr::filter(Presenc == 1)
iter1_pts <- make_pts('3m-presence-iter1', as.character(iter1$Date), as.character(iter1$Site), st_geometry(iter1))

## 2. 2026 opportunistic ground-truthing presences
gt_raw <- read.csv(file.path(p2proj, 'data', 'GroundTruthing', '2026_groundtruth_points.csv')) |>
  dplyr::mutate(Presence = as.integer(Plants > 0))
opp_2026 <- gt_raw |>
  dplyr::filter(Presence == 1) |>
  st_as_sf(coords = c('Long', 'Lat'), crs = 4326) |>
  st_transform(32613)
opp_2026_pts <- make_pts('2026-opportunistic', as.character(opp_2026$Date), NA_character_, st_geometry(opp_2026))

## 3. 2026 visited random points (GRTS design draws that turned out occupied)
visited <- st_read(file.path(p2proj, 'data', 'GroundTruthPts', 'visited_subset.gpkg'),
                    layer = 'visited_points', quiet = TRUE) |>
  dplyr::mutate(Presence = as.integer(dplyr::if_any(c(Prsnc_M, Prsnc_J, Prsnc_S), ~ replace_na(. > 0, FALSE)))) |>
  dplyr::filter(Presence == 1) |>
  st_transform(32613)
visited_pts <- make_pts('2026-visited', NA_character_, as.character(visited$Site), st_geometry(visited))

## 4. 2026 transect stations (core "known" + exploratory "other" transects) that turned
## out occupied - the field-filled counterpart of results/GroundTruthSampling/*-transects.gpkg's
## planned stations, entered directly into ground_truth_median-VISITED.gpkg.
transect_pres <- purrr::map_dfr(c('transects_known', 'transects_other'), function(lyr) {
  st_read(file.path(p2proj, 'data', 'GroundTruthPts', 'ground_truth_median-VISITED.gpkg'),
          layer = lyr, quiet = TRUE) |>
    dplyr::mutate(Presence = as.integer(dplyr::if_any(c(Prsnc_M, Prsnc_J, Prsnc_S), ~ replace_na(. > 0, FALSE))),
                  layer = lyr) |>
    dplyr::filter(Presence == 1)
}) |>
  st_transform(32613)
transect_pts <- make_pts('2026-transect', NA_character_, paste0(transect_pres$layer, '-patch', transect_pres$patch_ID),
                          st_geometry(transect_pres))

## occurrences_coloradense/occurrences.shp (the raw herbarium/iNat compilation, used
## elsewhere as `known_occ`) is deliberately NOT included here: 100 of its 125 points are
## exact (0m) duplicates of an iter1 presence, and the other 21 are old (1951-2012),
## imprecisely-geolocated herbarium records - the kind manually screened out of the
## modelling presence set for low geolocation quality (see flag_bad_geocodes_2026.R's
## header re: "16 herbarium records removed due to low geolocation quality"). Neither
## adds anything trustworthy for drawing population/spatial-block boundaries.

all_presences <- dplyr::bind_rows(iter1_pts, opp_2026_pts, visited_pts, transect_pts)

message(sprintf('%d total presence points collated:', nrow(all_presences)))
print(table(all_presences$source))

out_path <- file.path(p2proj, 'data', 'Data4modelling', 'all_presences_collated.gpkg')
st_write(all_presences, out_path, append = FALSE, quiet = TRUE)
message('Wrote ', out_path)
