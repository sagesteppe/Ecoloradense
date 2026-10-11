## Generates a template for the user to VERIFY the mapping between training
## population codes (`Lctn_bb`: XTAL, BV, CD, HB, ...) and 2026 field site names
## (Crystal, Bellview, Coche, ...).
##
## Why this is needed: the two naming systems don't join. Only 'Frigid' matches
## by spelling, and even that is 851 m from its nearest training population
## (HB), so the name match is itself suspect. Everything downstream that needs
## population identity currently falls back to a PROXIMITY guess:
##   - validate_against_2026_groundtruth.R's KNOWN/NEW regime split
##   - densityModeller()'s per-population arithmetic mean baseline
##   - knn_density_null()'s within-population variant
## all of which collapse to global/across-population behaviour without a real
## crosswalk. Proximity is a reasonable stand-in but it cannot distinguish "the
## same population resurveyed" from "a different population 300 m away", and
## only the user knows which.
##
## Writes TWO templates so the mapping can be checked from both directions - a
## one-way list hides both unmatched training populations and 2026 sites that
## are genuinely new:
##   1. ..._2026_to_training_TEMPLATE.csv - one row per 2026 site, with its 3
##      nearest training populations as candidates.
##   2. ..._training_to_2026_TEMPLATE.csv - one row per training population, so
##      populations with no 2026 visit are visible rather than absent.
##
## Fill in the columns prefixed VERIFY_ and leave the rest as context.
## Coordinates are WGS84 decimal degrees for easy pasting into a map.
##
## Run: Rscript Census_uncertainty/make_population_crosswalk_template.R

suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(terra)
})

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')
source('functions.R')

bn <- '1-3arc-Iteration1-PA1:1'
p2proc <- file.path(PROJ_ROOT, 'data', 'spatial', 'processed')
out_dir <- file.path(PROJ_ROOT, 'data', 'GroundTruthing')
n_candidates <- 3

## 1. The two population sets ----------------------------------------------

train_sf <- wrapper(paste0(bn, 'DO:0-Seed1021-Pr.tif'), return_early = TRUE)

# Earthslide is a separate taxon and Upper Crossing a bad geocode, but BOTH are
# kept here: the point of this template is to verify identity, and the user may
# want to record what those sites actually correspond to rather than have them
# silently missing.
gt <- sf::st_read(file.path(PROJ_ROOT, 'data', 'GroundTruthing', '2026_groundtruth_points.gpkg'),
                  quiet = TRUE) |>
  dplyr::filter(!is.na(Plants))

dem <- terra::rast(file.path(p2proc, 'dem_1-3arc', 'dem.tif'))

#' Collapse a plot-level sf to one row per population, with summary context.
#' @param x sf of plots; @param id_col grouping column; @param count_col counts.
summarize_pops <- function(x, id_col, count_col){
  x$.id <- x[[id_col]]
  x$.ct <- x[[count_col]]
  x$.elev <- terra::extract(dem, terra::vect(x))[, 2]

  centroids <- x |>
    dplyr::group_by(.id) |>
    dplyr::summarize(.groups = 'drop') |>
    sf::st_centroid()

  stats <- x |>
    sf::st_drop_geometry() |>
    dplyr::group_by(.id) |>
    dplyr::summarize(
      n_plots     = dplyr::n(),
      n_occupied  = sum(.ct > 0),
      mean_count  = round(mean(.ct), 1),
      median_count = stats::median(.ct),
      max_count   = max(.ct),
      elev_m      = round(mean(.elev, na.rm = TRUE)),
      .groups = 'drop'
    )

  ll <- centroids |> sf::st_transform(4326) |> sf::st_coordinates()
  dplyr::left_join(centroids, stats, by = '.id') |>
    dplyr::mutate(lon = round(ll[, 1], 5), lat = round(ll[, 2], 5))
}

train_pops <- summarize_pops(train_sf, 'Lctn_bb', 'Prsnc_All')
gt_sites   <- summarize_pops(gt, 'Site', 'Plants')

message(sprintf('%d training populations, %d 2026 field sites.',
                nrow(train_pops), nrow(gt_sites)))

## 2. Nearest-neighbour candidates ------------------------------------------

#' For each row of `from`, the `k` nearest rows of `to`, flattened into
#' candidate_1..k / dist_1..k_m columns.
nearest_candidates <- function(from, to, k, prefix){
  d <- sf::st_distance(from, to)
  units(d) <- NULL
  out <- lapply(seq_len(nrow(from)), function(i){
    ord <- order(d[i, ])[seq_len(min(k, ncol(d)))]
    stats::setNames(
      as.list(c(rbind(to$.id[ord], round(d[i, ord])))),
      as.vector(rbind(paste0(prefix, '_', seq_along(ord)),
                      paste0('dist_', seq_along(ord), '_m')))
    )
  })
  dplyr::bind_rows(out)
}

## 3. Template 1: 2026 site -> training population --------------------------

t1 <- dplyr::bind_cols(
  sf::st_drop_geometry(gt_sites) |>
    dplyr::rename(Site_2026 = .id, n_plots_2026 = n_plots, n_occupied_2026 = n_occupied,
                  mean_count_2026 = mean_count, median_count_2026 = median_count,
                  max_count_2026 = max_count, elev_m_2026 = elev_m,
                  lon_2026 = lon, lat_2026 = lat),
  nearest_candidates(gt_sites, train_pops, n_candidates, 'cand_train_pop')
) |>
  dplyr::mutate(
    VERIFY_train_pop       = '',   # the Lctn_bb code this site IS, or NEW_POPULATION, or UNSURE
    VERIFY_same_population = '',   # yes / no / unsure
    VERIFY_notes           = ''
  ) |>
  dplyr::arrange(dist_1_m)

## 4. Template 2: training population -> 2026 site --------------------------

t2 <- dplyr::bind_cols(
  sf::st_drop_geometry(train_pops) |>
    dplyr::rename(Lctn_bb = .id, n_plots_train = n_plots, n_occupied_train = n_occupied,
                  mean_count_train = mean_count, median_count_train = median_count,
                  max_count_train = max_count, elev_m_train = elev_m,
                  lon_train = lon, lat_train = lat),
  nearest_candidates(train_pops, gt_sites, n_candidates, 'cand_2026_site')
) |>
  dplyr::mutate(
    VERIFY_2026_site  = '',   # the 2026 Site name, or NOT_VISITED_2026, or UNSURE
    VERIFY_notes      = ''
  ) |>
  dplyr::arrange(dist_1_m)

## 5. Write ----------------------------------------------------------------

f1 <- file.path(out_dir, 'population_crosswalk_2026_to_training_TEMPLATE.csv')
f2 <- file.path(out_dir, 'population_crosswalk_training_to_2026_TEMPLATE.csv')
write.csv(t1, f1, row.names = FALSE, na = '')
write.csv(t2, f2, row.names = FALSE, na = '')

cat('\n=== Template 1: 2026 site -> training population ===\n')
print(as.data.frame(t1[, c('Site_2026', 'n_plots_2026', 'mean_count_2026', 'elev_m_2026',
                           'cand_train_pop_1', 'dist_1_m', 'cand_train_pop_2', 'dist_2_m')]),
      row.names = FALSE)
cat('\n=== Template 2: training population -> 2026 site ===\n')
print(as.data.frame(t2[, c('Lctn_bb', 'n_plots_train', 'mean_count_train', 'elev_m_train',
                           'cand_2026_site_1', 'dist_1_m', 'cand_2026_site_2', 'dist_2_m')]),
      row.names = FALSE)

message('\nWritten: ', f1)
message('Written: ', f2)
message('Fill the VERIFY_ columns and hand back; nothing else needs editing.')
