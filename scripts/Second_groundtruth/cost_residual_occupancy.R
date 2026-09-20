## H3 (Occupancy Based Field Sampling) - tests whether Euclidean distance (E) and the
## region-orthogonalized cost-distance anomaly (R) from the ground-truth design scripts
## (sdm_groundtruth_design_1-3arc_Seed*.R) are actually useful predictors of observed
## occupancy at the 2026 ground-truth points. This is the analysis SDManuscript.Rmd's
## section 3.3 placeholder describes:
##   "A [...] logistic model identifies that Euclidean distance has X contributions,
##    while the cost-distance anomaly surface has X contribution."
## and the sub-population/new-population detection table just above it, using the same
## 270m / 450m distance bands as the Introduction's "Here we" list (extension < 270m,
## sub-population 270-450m, new population >= 450m).
##
## R is recomputed here (not read from the design script) since the GRTS output only
## keeps siteID/coords/presence, not the E/Cost/R covariates used to select it. Uses the
## Seed1021 realization of 1-3arc-Iteration1-PA1:1DO:0, the model version the design
## scripts built the sampling frame from; Seed1061 is the other available realization if
## a robustness check against the alternate model fit is wanted later.

suppressPackageStartupMessages({
  library(sf)
  library(terra)
  library(dplyr)
  library(mgcv)
  library(tidyr)
  library(readr)
})

p2proj <- '/media/steppe/hdd/EriogonumColoradenseTaxonomy'

cfg <- list(
  ver         = '1-3arc-Iteration1-PA1:1DO:0-Seed1021',
  cost_thresh = 'spec_sens',
  band_near   = 270,   # < this: extension of a known population
  band_far    = 450    # >= this: new population; in between: sub-population
)

## 1. Ground-truth points (same construction as digitize_evaluate_2026_groundtruth.R) ----

gt_raw <- read.csv(file.path(p2proj, 'data', 'GroundTruthing', '2026_groundtruth_points.csv')) |>
  rename(Plot = `Plot..if.tissue.`) |>
  mutate(Date = as.Date(Date, format = '%m/%d/%y'),
         Presence = as.integer(Plants > 0))

gt <- gt_raw |>
  st_as_sf(coords = c('Long', 'Lat'), crs = 4326) |>
  st_transform(32613)

visited_gt <- st_read(file.path(p2proj, 'data', 'GroundTruthPts', 'visited_subset.gpkg'),
                       layer = 'visited_points', quiet = TRUE) |>
  st_transform(st_crs(gt)) |>
  mutate(Presence = as.integer(if_any(c(Prsnc_M, Prsnc_J, Prsnc_S), ~ replace_na(. > 0, FALSE))))
visited_gt <- st_sf(Presence = visited_gt$Presence, geometry = st_geometry(visited_gt))

gt <- bind_rows(gt |> select(Presence), visited_gt)

## 2. Rebuild the design scripts' S / E / Cost / region stack for this model version -----

suit_path  <- file.path(p2proj, 'results', 'suitability_maps', paste0(cfg$ver, '-Pr.tif'))
aoa_path   <- file.path(p2proj, 'results', 'suitability_maps', paste0(cfg$ver, '-AOA.tif'))
aniso_path <- file.path(p2proj, 'results', 'cost_distances',
                        paste0(cfg$ver, '-', cfg$cost_thresh, '-costDistAniso.tif'))
iso_path   <- file.path(p2proj, 'results', 'cost_distances',
                        paste0(cfg$ver, '-', cfg$cost_thresh, '-costDist.tif'))

S_mean <- terra::rast(suit_path)
aoa_r  <- terra::rast(aoa_path)
Cost <- if (file.exists(aniso_path)) terra::rast(aniso_path) else terra::rast(iso_path)
terra::crs(Cost) <- terra::crs(S_mean)

region_v <- sf::st_read(file.path(p2proj, 'data', 'hikingTrails', 'GroundTruch_Areas.gpkg'), quiet = TRUE)
region_v$region_id <- as.integer(factor(region_v$Site))
region <- terra::rasterize(terra::vect(region_v), S_mean, field = 'region_id')

known_occ <- sf::st_read(file.path(p2proj, 'data', 'Data4modelling', '3m-presence-iter1.gpkg'), quiet = TRUE) |>
  dplyr::filter(Presenc == 1) |>
  sf::st_transform(terra::crs(S_mean))
E <- terra::distance(S_mean, terra::vect(known_occ))
names(E) <- 'E'

S_mean <- terra::mask(S_mean, aoa_r, maskvalues = 0) |> terra::mask(region)
E      <- terra::mask(E,      aoa_r, maskvalues = 0) |> terra::mask(region)
Cost   <- terra::mask(Cost,   aoa_r, maskvalues = 0) |> terra::mask(region)

cost_max      <- terra::global(Cost, 'max', na.rm = TRUE)[1, 1]
cost_gap_mask <- is.na(Cost) & !is.na(region)
if (terra::global(cost_gap_mask, 'sum', na.rm = TRUE)[1, 1] > 0) {
  Cost <- terra::ifel(cost_gap_mask, cost_max, Cost)
}

stk <- c(S_mean, E, Cost, region)
names(stk) <- c('S', 'E', 'Cost', 'region')

## Region-aware cost orthogonalization: R = within-region residual of Cost on E - same
## logic/guardrails as sdm_groundtruth_design_1-3arc_Seed1021.R.
df <- as.data.frame(stk, cells = TRUE, na.rm = TRUE) |>
  dplyr::mutate(region = factor(region, labels = levels(factor(region_v$region_id)))) |>
  dplyr::filter(!is.na(region)) |>
  dplyr::group_by(region) |>
  dplyr::group_modify(~{
    if (nrow(.x) >= 30 && length(unique(.x$E)) > 10) {
      mm <- tryCatch(gam(Cost ~ s(E, k = 6), data = .x, method = 'REML'),
                     error = function(e) lm(Cost ~ poly(E, 2), data = .x))
      .x$R <- resid(mm)
    } else .x$R <- as.numeric(scale(.x$Cost - .x$E))
    .x
  }) |>
  dplyr::ungroup()

## Rasterize R back onto the full grid so it (and E) can be sampled at the ground-truth
## point locations exactly like any other covariate.
R_rast <- terra::rast(S_mean)
R_rast[df$cell] <- df$R
names(R_rast) <- 'R'

## 3. Extract E/R at the ground-truth points, band-classify presences -------------------

pt_covs <- terra::extract(c(E, R_rast), terra::vect(gt))
pts <- gt |>
  dplyr::mutate(E = pt_covs$E, R = pt_covs$R) |>
  sf::st_drop_geometry()

## Points outside every named GroundTruth Site region have no region-relative R (masked
## to NA above, by design - see sdm_groundtruth_design_1-3arc_Seed1021.R's Details) and
## are dropped from the glm() below along with any other NA covariate.
message(sprintf('%d of %d ground-truth points fall inside a named region and have E/R.',
                sum(!is.na(pts$R)), nrow(pts)))

band_table <- pts |>
  dplyr::filter(Presence == 1, !is.na(E)) |>
  dplyr::mutate(band = dplyr::case_when(
    E < cfg$band_near ~ 'extension (<270m)',
    E < cfg$band_far  ~ 'sub-population (270-450m)',
    TRUE              ~ 'new population (>=450m)'
  )) |>
  dplyr::count(band, name = 'n_presences')

## 4. Logistic model of occupancy on E and R, plus each term's likelihood-ratio contribution

model_data <- pts |> tidyr::drop_na(Presence, E, R)

fit_null <- glm(Presence ~ 1,     family = binomial, data = model_data)
fit_E    <- glm(Presence ~ E,     family = binomial, data = model_data)
fit_R    <- glm(Presence ~ R,     family = binomial, data = model_data)
fit_full <- glm(Presence ~ E + R, family = binomial, data = model_data)

lrt_E_alone    <- anova(fit_null, fit_E,    test = 'Chisq')
lrt_R_alone    <- anova(fit_null, fit_R,    test = 'Chisq')
lrt_E_given_R  <- anova(fit_R,    fit_full, test = 'Chisq')
lrt_R_given_E  <- anova(fit_E,    fit_full, test = 'Chisq')

glm_summary <- as.data.frame(summary(fit_full)$coefficients) |>
  tibble::rownames_to_column('term') |>
  rename(estimate = Estimate, std_error = `Std. Error`, statistic = `z value`, p_value = `Pr(>|z|)`) |>
  mutate(model = 'full: Presence ~ E + R', .before = 1)

lrt_summary <- tibble(
  term = c('E (alone)', 'R (alone)', 'E (given R)', 'R (given E)'),
  deviance = c(lrt_E_alone$Deviance[2], lrt_R_alone$Deviance[2],
               lrt_E_given_R$Deviance[2], lrt_R_given_E$Deviance[2]),
  df = c(lrt_E_alone$Df[2], lrt_R_alone$Df[2], lrt_E_given_R$Df[2], lrt_R_given_E$Df[2]),
  p_value = c(lrt_E_alone$`Pr(>Chi)`[2], lrt_R_alone$`Pr(>Chi)`[2],
              lrt_E_given_R$`Pr(>Chi)`[2], lrt_R_given_E$`Pr(>Chi)`[2])
)

print(band_table)
print(glm_summary)
print(lrt_summary)

## 5. Write outputs -----------------------------------------------------------------------

dir.create(file.path(p2proj, 'results', 'tables'), showWarnings = FALSE, recursive = TRUE)
write.csv(band_table, file.path(p2proj, 'results', 'tables', '2026-groundtruth-occupancy-bands.csv'), row.names = FALSE)
write.csv(glm_summary, file.path(p2proj, 'results', 'tables', '2026-groundtruth-occupancy-glm-coefficients.csv'), row.names = FALSE)
write.csv(lrt_summary, file.path(p2proj, 'results', 'tables', '2026-groundtruth-occupancy-glm-lrt.csv'), row.names = FALSE)
