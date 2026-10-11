## Scores every Iteration1 count/density candidate fitted on the full training
## set (bn = 1-3arc-Iteration1-PA1:1), plus a full baseline ladder, against the
## 2026 field census (data/GroundTruthing/2026_groundtruth_points.gpkg).
##
## This is the single scoring path for the 2026 field data. It supersedes the
## separate summarize_2026_groundtruth_by_phase.R (merged in 2026-10-10); that
## script duplicated ~20% of this one verbatim, so the two were collapsed rather
## than kept in sync. Note the baseline ladder in section 3 is still a second
## implementation of the nulls inside densityModeller() (functions.R ~2234-2330)
## - deliberately left duplicated for now, since factoring it out would mean
## touching densityModeller() and regenerating its cached split-comparison
## outputs. Worth revisiting once the manuscript numbers settle.
##
## EXCLUSIONS (both are permanent - do not score these plots):
##   - Earthslide: reassigned to a separate taxon (E. lonchophyllum var.
##     tristichum); its data lives in 2026_Earthslide_is_tristichum.csv, not the
##     gpkg. Guarded below in case it is ever re-merged.
##   - Upper Crossing: the site was targeted off a bad iNaturalist geocode, so
##     the plant was never recorded there. All 5 plots are true zeros at a place
##     that was never a population, which hands free accuracy to any model that
##     predicts low - it measures the geocode error, not the model.
##
## TWO CUTS, and the second is the one that matters:
##   - All plots (177): reported, but NOT the headline. ~46% of plots are true
##     zeros, so for a rare species the MAE-optimal constant predictor is 0, and
##     "predict zero everywhere" beats most of the ladder by construction. A low
##     all-plot MAE mostly demonstrates that the species is rare.
##   - Occupied plots only (96; the `_occ` columns): given the plant is actually
##     present, how many are there? This is the question the count models exist
##     to answer, and the only cut on which they can be fairly judged.
##
## Elevational phase ('High'/'Low') is the natural break in site elevations: the
## ten scored sites cluster at 2751-2838 m and 3311-3700 m with a ~470 m gap.
## Phase matters because the two are on completely different count scales
## (occupied median 6 plants high vs 60.5 low), so their absolute MAEs are not
## comparable - hence the scale-free MdAPE columns.
##
## Extraction pattern (raster -> points, CRS) follows
## Second_groundtruth/digitize_evaluate_2026_groundtruth.R.
##
## Run: Rscript Census_uncertainty/validate_against_2026_groundtruth.R

suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(terra); library(brms); library(Metrics)
  library(parsnip); library(workflows); library(xgboost); library(bonsai); library(lightgbm)
})

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')
source('functions.R')
source('Census_uncertainty/density_bayes.R')

bn <- '1-3arc-Iteration1-PA1:1'
fp <- file.path('..', 'results', 'count_models')
p2proc <- file.path(PROJ_ROOT, 'data', 'spatial', 'processed')

## Phase B (per-elevational-phase) models produced by
## First_modelling/abundance_phaseB_run.R. If these files do not yet exist
## the Phase B rows are simply absent from the output; Phase A scoring is
## unaffected. See the driver script for the LOW_ELEV_POPS justification.
LOW_ELEV_POPS <- c('CD')
bn_phaseB <- c(High = '1-3arc-abundB-highelev', Low = '1-3arc-abundB-lowelev')

## Sites dropped before scoring, with the reason (see header).
excluded_sites <- c(
  'Earthslide'     = 'separate taxon (E. lonchophyllum var. tristichum)',
  'Upper Crossing' = 'site targeted off a bad iNaturalist geocode; never a population'
)

## Elevational phase membership (see header note on the ~470 m break).
phase_lookup <- c(
  'Coche Highway'  = 'Low',  'Antero'         = 'Low',  'Coche'    = 'Low',
  'Bellview 2'     = 'High', 'Bellview'       = 'High', 'Treasury' = 'High',
  'White Rock'     = 'High', 'Crystal'        = 'High', 'Frigid'   = 'High',
  'Virginia Basin' = 'High'
)

conformal_level <- 0.9

## VERIFIED population crosswalk (user, 2026-10-10): which 2026 field sites are
## resurveys of a training population, by identity rather than proximity. Any
## site absent from this list was a fresh survey of new ground.
##
## This replaces an earlier distance-based guess, which was WRONG on the single
## largest site: Frigid sits 851 m from training population HB and 1018 m from
## WMB, so proximity classed it as a resurvey, but it is in fact fresh ground -
## as is Treasury (493 m from YP). Note also that the training population
## literally named 'Frigid' is 2213 m from the 2026 Frigid site and is a
## different population, so the one name-identical pair in the data never
## matched either. Consequence worth carrying into the SAC discussion:
## distinct populations of this species occur <1 km apart, so geographic
## proximity does not imply population identity - which is exactly the
## assumption `knn_density_null()`'s within-population variant and distance-
## based spatial CV blocking rely on.
population_crosswalk <- c(
  'Crystal'    = 'XTAL',
  'Coche'      = 'CD',
  'Bellview'   = 'BV',
  'Bellview 2' = 'BV',
  # An extension of the existing Frigid population rather than new ground, per
  # the user - despite the 2213 m separation between the two centroids, which is
  # why proximity to the nearest population (HB, 851 m) was the wrong signal.
  'Frigid'     = 'Frigid'
)

## 1. Assemble the held-out covariate table -------------------------------

gt <- sf::st_read(file.path(PROJ_ROOT, 'data', 'GroundTruthing', '2026_groundtruth_points.gpkg'),
                  quiet = TRUE) |>
  dplyr::filter(!is.na(Plants))

for(s in names(excluded_sites)){
  hit <- which(tolower(gt$Site) == tolower(s))
  if(length(hit) > 0){
    message(sprintf('Dropping %d %s plot(s): %s.', length(hit), s, excluded_sites[[s]]))
    gt <- gt[-hit, ]
  }
}

unassigned <- setdiff(unique(gt$Site), names(phase_lookup))
if(length(unassigned) > 0){
  stop('Site(s) with no elevational phase assignment: ', paste(unassigned, collapse = ', '),
       ' - add them to phase_lookup.')
}

rast_dat <- rastReader('dem_1-3arc', p2proc)
suit_r <- terra::rast(file.path(PROJ_ROOT, 'results', 'suitability_maps',
                                paste0(bn, 'DO:0-Seed1021-Pr.tif')))

coords <- sf::st_coordinates(gt)
df <- dplyr::bind_cols(
  Observed = gt$Plants,
  Site = gt$Site,
  dplyr::select(terra::extract(rast_dat, terra::vect(gt)), -ID),
  Pr.SuitHab = terra::extract(suit_r, terra::vect(gt))[, 2]
)
# rast_dat carries per-cell Longitude/Latitude surfaces; overwrite with the
# point's exact coordinates, matching wrapper()'s convention (functions.R)
# rather than the raster's cell-binned approximation.
df$Longitude <- coords[, 1]
df$Latitude  <- coords[, 2]

n_before <- nrow(df)
df <- tidyr::drop_na(df)
df$point_id <- seq_len(nrow(df))
message(sprintf('%d of %d field plots retained after covariate extraction (%d dropped - outside raster extent or missing covariate).',
                nrow(df), n_before, n_before - nrow(df)))

df$Phase <- unname(phase_lookup[df$Site])
message('Plots by phase: ', paste(names(table(df$Phase)), table(df$Phase), sep = ' = ', collapse = ', '))
message(sprintf('Occupied (Observed > 0): %d of %d plots (%d high-elevation, %d low-elevation).',
                sum(df$Observed > 0), nrow(df),
                sum(df$Observed > 0 & df$Phase == 'High'),
                sum(df$Observed > 0 & df$Phase == 'Low')))

## Every scored candidate returns this shape. ScoredBy marks the training
## stack (pooled Phase A vs per-phase Phase B); downstream summaries group on
## it so the two stacks' rows never collapse together.
pred_row <- function(label, preds, scored_by = 'Phase A (pooled)'){
  data.frame(ScoredBy = scored_by, Model = label,
             point_id = df$point_id, Site = df$Site, Phase = df$Phase,
             Observed = df$Observed, Predicted = as.numeric(preds))
}

## 2. Score the ML candidates ----------------------------------------------

# Site/Phase/point_id/Observed ride along in `new_data`; parsnip ignores
# columns the recipe doesn't ask for.
score_ml <- function(model_suffix, label){
  f <- file.path(fp, 'models', paste0(bn, '-', model_suffix, '.rds'))
  if(!file.exists(f)){
    message('Missing model, skipped: ', f)
    return(NULL)
  }
  fit <- readRDS(f)$Model
  pred_row(label, stats::predict(fit, new_data = df)[[1]])
}

## Phase B: route each plot to the model trained on its own elevational
## phase; a missing phase model leaves that phase's plots NA (and warns once).
## The two phases together cover every plot - a plot with no phase assignment
## would have been caught upstream.
score_ml_phaseB <- function(model_suffix, label){
  out <- rep(NA_real_, nrow(df))
  any_hit <- FALSE
  for(ph in names(bn_phaseB)){
    f <- file.path(fp, 'models', paste0(bn_phaseB[[ph]], '-', model_suffix, '.rds'))
    if(!file.exists(f)){
      message(sprintf('Phase B %s-elev %s missing, those plots left NA.', ph, label))
      next
    }
    any_hit <- TRUE
    mask <- df$Phase == ph
    out[mask] <- as.numeric(stats::predict(readRDS(f)$Model, new_data = df[mask, ])[[1]])
  }
  if(!any_hit) return(NULL)
  pred_row(label, out, scored_by = 'Phase B (per-phase)')
}

ml_candidates <- c(
  'poisson_spat'      = 'XGB Poisson Spat.',
  'poisson'           = 'XGB Poisson',
  'tweedie_spat'      = 'XGB Tweedie Spat.',
  'tweedie'           = 'XGB Tweedie',
  'lgbm-poisson_spat' = 'LGBM Poisson Spat.'
)
ml_results <- dplyr::bind_rows(
  purrr::imap_dfr(ml_candidates, ~ score_ml(.y, .x)),
  purrr::imap_dfr(ml_candidates, ~ score_ml_phaseB(.y, .x))
)

## 3. Baseline ladder ------------------------------------------------------

## A reader should be able to tell whether the problem warranted a fitted model
## at all. The ladder runs: constant (no information) -> single-covariate GLM ->
## spatial interpolation of observed counts -> the fitted candidates in sections
## 2 and 4. A candidate that cannot clear the rung below it is not earning its
## complexity.
##
## NOTE on `Lctn_bb`: training populations are abbreviations (XTAL, BV, CD, HB,
## ...) while the 2026 gpkg uses full site names, so any name-based
## per-population join is near-vacuous here - only 'Frigid' matches by spelling,
## and even that is 851 m from its nearest training population (HB). So
## densityModeller()'s per-population arithmetic mean collapses to the global
## training mean on this test set, as does k-NN within-population. Both are
## reported below as what they actually are. A verified name crosswalk would
## make the per-population rungs meaningful; until then they are not claimed as
## such. Separately: measured by distance, 6 of 10 sites are within 1.2 km of a
## training population, so most of this test set is RESURVEYS, not genuinely new
## populations - only Antero (86 km out) is a true extrapolation test.

train_full_sf <- wrapper(paste0(bn, 'DO:0-Seed1021-Pr.tif'), return_early = TRUE)
train_df <- sf::st_drop_geometry(train_full_sf)

## Distance from each 2026 site to the nearest training population centroid.
## This is the basis for the KNOWN/NEW regime split in section 6 and is
## PROXIMITY-BASED, not identity-based - the `Lctn_bb` name join is unusable
## (see the note above), so this stands in for a verified crosswalk. Replace
## with the crosswalk when available.
pop_centroids <- train_full_sf |>
  dplyr::group_by(Lctn_bb) |>
  dplyr::summarize(train_n = dplyr::n(), .groups = 'drop') |>
  sf::st_centroid()

site_centroids <- sf::st_as_sf(
  dplyr::distinct(df, Site, Longitude, Latitude),
  coords = c('Longitude', 'Latitude'), crs = 32613
) |>
  dplyr::group_by(Site) |>
  dplyr::summarize(.groups = 'drop') |>
  sf::st_centroid()

site_dist_m <- sf::st_distance(site_centroids, pop_centroids)
units(site_dist_m) <- NULL
site_dist <- data.frame(
  Site = site_centroids$Site,
  dist_m = round(apply(site_dist_m, 1, min)),
  nearest_train_pop = pop_centroids$Lctn_bb[apply(site_dist_m, 1, which.min)]
)
message('Site distance to nearest training population (m):')
message(paste(utils::capture.output(print(site_dist[order(site_dist$dist_m), ], row.names = FALSE)),
              collapse = '\n'))

## 3a. Constant predictors -------------------------------------------------
## The training mean is the arithmetic-mean null. The training median is the
## more important rung: MAE is minimized by the median, so "predict the training
## median" is the floor any all-plot MAE comparison must clear. Here that median
## is 0, making this rung identical to an all-zero predictor - stated outright
## rather than reported twice. For a rare species this rung SHOULD win on the
## all-plot cut; that is why the occupied-only cut exists.
train_mean   <- mean(train_df$Prsnc_All)
train_median <- stats::median(train_df$Prsnc_All)
message(sprintf('Training set: %d plots, %d populations, mean %.2f, median %.0f.',
                nrow(train_df), dplyr::n_distinct(train_df$Lctn_bb), train_mean, train_median))

constant_results <- dplyr::bind_rows(
  pred_row(sprintf('Constant: training mean (%.2f)', train_mean), rep(train_mean, nrow(df))),
  pred_row(sprintf('Constant: training median (%.0f)', train_median), rep(train_median, nrow(df)))
)

## 3b. Single-covariate Poisson GLMs ---------------------------------------
## Pr.SuitHab is the Aim 1 SDM's own output, so this rung answers a question the
## manuscript needs anyway: does the suitability surface carry count
## information, or does it only predict presence? Elevation is included because
## the phase split is an elevation split - if an elevation-only GLM reproduces
## the phase-level MAEs, the covariate-rich models are adding little.
glm_baseline <- function(formula, label){
  fit <- try(stats::glm(formula, family = stats::poisson(), data = train_df), silent = TRUE)
  if(inherits(fit, 'try-error')){
    message('GLM baseline failed, skipped: ', label)
    return(NULL)
  }
  pred_row(label, stats::predict(fit, newdata = df, type = 'response'))
}

glm_results <- dplyr::bind_rows(
  glm_baseline(Prsnc_All ~ Pr.SuitHab, 'GLM Poisson: suitability only'),
  glm_baseline(Prsnc_All ~ elevation,  'GLM Poisson: elevation only')
)

## 3c. Spatial interpolation of observed counts ----------------------------
## densityModeller() labelled its spatial null "Kriging", but gstat::krige()
## with no `model=` argument does inverse distance weighting at idp = 2.0, not
## kriging (see ?krige) - so IDW is what the pipeline reported all along and
## ordinary kriging had never been scored. functions.R was corrected on
## 2026-10-10; both are run here, labelled for what they are.
##
## Kriging is also attempted on log1p(count): ordinary kriging assumes an
## approximately stationary Gaussian field, and this response runs 0-150 with a
## median of 0, badly violating that. expm1() back-transforms to a median-ish
## rather than mean-ish prediction, appropriate where the MAE target is the
## median anyway.
## Lctn_bb carries the VERIFIED training population code for resurveyed sites,
## so knn_density_null()'s within-population variant can find same-population
## training rows. Fresh-survey sites keep their own site name, which matches no
## training population, so within-population is undefined there and correctly
## falls back to across-population.
test_lctn <- dplyr::coalesce(unname(population_crosswalk[df$Site]), df$Site)

# A fresh-survey site whose name happens to equal a training population code
# would be silently treated as a within-population lookup. That is exactly the
# 'Frigid' case, which IS a real match and so is declared in the crosswalk -
# but an undeclared collision would be a bug, so fail loudly instead.
undeclared <- setdiff(intersect(df$Site, train_df$Lctn_bb), names(population_crosswalk))
if(length(undeclared) > 0){
  stop('2026 site name(s) collide with a training population code but are not in ',
       'population_crosswalk: ', paste(undeclared, collapse = ', '),
       ' - declare the match or rename, do not leave it implicit.')
}

test_sf <- sf::st_as_sf(
  data.frame(Prsnc_All = df$Observed, Lctn_bb = test_lctn, Pr.SuitHab = df$Pr.SuitHab,
             Longitude = df$Longitude, Latitude = df$Latitude),
  coords = c('Longitude', 'Latitude'), crs = 32613, remove = FALSE
)
message(sprintf('k-NN within-population: %d of %d plots have a matched training population (%s).',
                sum(df$Site %in% names(population_crosswalk)), nrow(df),
                paste(sort(unique(stats::na.omit(unname(population_crosswalk[df$Site])))), collapse = ', ')))

# gstat needs unique training locations (duplicates make the kriging covariance
# matrix singular); co-located plots are averaged. IDW needs no inversion, but
# the same deduplicated frame is used for both so they stay comparable.
train_xy <- sf::st_coordinates(train_full_sf)
train_krige <- data.frame(Prsnc_All = train_full_sf$Prsnc_All,
                          X = train_xy[, 1], Y = train_xy[, 2]) |>
  dplyr::group_by(X, Y) |>
  dplyr::summarize(Prsnc_All = mean(Prsnc_All), .groups = 'drop') |>
  sf::st_as_sf(coords = c('X', 'Y'), crs = sf::st_crs(train_full_sf))
if(nrow(train_krige) < nrow(train_full_sf)){
  message(sprintf('Kriging/IDW: collapsed %d co-located training plots to %d unique locations.',
                  nrow(train_full_sf), nrow(train_krige)))
}

idw_pred <- gstat::krige(Prsnc_All ~ 1, train_krige, newdata = test_sf, debug.level = 0)$var1.pred

krige_fit <- function(response_fn, back_fn, label){
  d <- train_krige
  d$resp <- response_fn(d$Prsnc_All)
  emp <- try(gstat::variogram(resp ~ 1, d), silent = TRUE)
  if(inherits(emp, 'try-error')) return(NULL)
  vgm_model <- try(gstat::fit.variogram(emp, gstat::vgm(c('Exp', 'Sph', 'Gau', 'Mat'))),
                   silent = TRUE)
  if(inherits(vgm_model, 'try-error') || attr(vgm_model, 'singular')){
    # A singular fit here is a pure-nugget variogram (structured psill -> 0): no
    # spatial autocorrelation left to exploit, so ordinary kriging degenerates
    # to the global mean and would duplicate the tier-1 constant rung. That
    # degeneracy is the result worth reporting, so it is stated rather than
    # scored as its own row.
    message(sprintf('%s: variogram fit is singular / pure nugget (structured psill %.3g, nugget %.3g) - no spatial structure to krige, rung omitted as it degenerates to the training mean.',
                    label,
                    if(inherits(vgm_model, 'try-error')) NA_real_ else vgm_model$psill[nrow(vgm_model)],
                    if(inherits(vgm_model, 'try-error')) NA_real_ else vgm_model$psill[1]))
    return(NULL)
  }
  # nugget vs structured psill is the number that says whether kriging can help
  # at all; a range far larger than the between-population spacing means the
  # "spatial" fit is really just a broad regional trend.
  message(sprintf('%s: fitted %s variogram - range %.0f m, nugget %.2f, structured psill %.2f (nugget is %.0f%% of total sill).',
                  label, as.character(vgm_model$model[nrow(vgm_model)]),
                  vgm_model$range[nrow(vgm_model)], vgm_model$psill[1],
                  vgm_model$psill[nrow(vgm_model)],
                  100 * vgm_model$psill[1] / sum(vgm_model$psill)))
  k <- gstat::krige(resp ~ 1, d, newdata = test_sf, model = vgm_model, debug.level = 0)
  pred_row(label, pmax(back_fn(k$var1.pred), 0))
}

spatial_results <- dplyr::bind_rows(
  pred_row('IDW (gstat default, idp=2)', idw_pred),
  krige_fit(identity, identity, 'Ordinary kriging (raw count)'),
  krige_fit(log1p, expm1, 'Ordinary kriging (log1p count)')
)

## 3d. k-NN null models ----------------------------------------------------
## Pure lookups over *observed* counts, so unlike the gp() term they
## structurally cannot blow up when extrapolating to an unmodelled population
## (finding_density_model_extrapolation_limits.md). Within-population falls back
## to across-population for every site except the single accidental 'Frigid'
## name match (see the Lctn_bb note above), so the two are near-identical here
## by construction, not because the distinction never matters.
knn_scored <- knn_density_null(train_full_sf, test_sf, k = 5)
knn_results <- dplyr::bind_rows(
  pred_row('kNN (within pop)', knn_scored$within_population$Predicted),
  pred_row('kNN (across pop)', knn_scored$across_population$Predicted)
)

baseline_results <- dplyr::bind_rows(constant_results, glm_results,
                                     spatial_results, knn_results)

## 3e. Phase B per-phase constant and GLM baselines ------------------------
## Spatial-interpolation and k-NN rungs are not re-fit per phase: with CD
## pulled out, the Low-phase training set has one population at a single
## cluster of points, so IDW/kriging/kNN collapse to its arithmetic mean -
## already captured by the per-phase constant rung below.
train_df$Phase <- ifelse(train_df$Lctn_bb %in% LOW_ELEV_POPS, 'Low', 'High')

build_phase_baselines <- function(phase){
  td   <- train_df[train_df$Phase == phase, ]
  mask <- df$Phase == phase
  if(sum(mask) == 0 || nrow(td) == 0) return(NULL)

  mn  <- mean(td$Prsnc_All); med <- stats::median(td$Prsnc_All)
  vmn <- rep(NA_real_, nrow(df)); vmn[mask] <- mn
  vmd <- rep(NA_real_, nrow(df)); vmd[mask] <- med

  rows <- list(
    pred_row(sprintf('Constant: %s-phase training mean (%.2f)',  phase, mn),  vmn,
             scored_by = 'Phase B (per-phase)'),
    pred_row(sprintf('Constant: %s-phase training median (%.0f)', phase, med), vmd,
             scored_by = 'Phase B (per-phase)')
  )

  for(fl in list(c('Pr.SuitHab', 'suitability only'),
                 c('elevation',  'elevation only'))){
    fo <- stats::as.formula(paste('Prsnc_All ~', fl[[1]]))
    fit <- try(stats::glm(fo, family = stats::poisson(), data = td), silent = TRUE)
    if(inherits(fit, 'try-error')) next
    v <- rep(NA_real_, nrow(df))
    v[mask] <- stats::predict(fit, newdata = df[mask, ], type = 'response')
    rows <- c(rows, list(pred_row(
      sprintf('GLM Poisson (%s-phase): %s', phase, fl[[2]]), v,
      scored_by = 'Phase B (per-phase)')))
  }
  dplyr::bind_rows(rows)
}
baseline_results <- dplyr::bind_rows(baseline_results,
                                     build_phase_baselines('High'),
                                     build_phase_baselines('Low'))

## 4. Score every brms family fitted on the full training set ---------------

## All fitted families are scored, not just the Phase-1-promoted Poisson: the
## promotion was made on internal CV, so restricting the field comparison to the
## promoted family would hide whether that promotion was the right call. Only
## full-training-set fits are eligible - the many
## `-classicRepN-`/`-populationRepX-`/`-population_groupRepX-` files in
## results/count_models/models/ are per-replicate refits on data *subsets* from
## the split-strategy comparison, so scoring them here would compare models
## trained on different amounts of data.
brms_families <- c('poisson', 'negbinomial', 'hurdlenegbinomial')

score_brms <- function(family_label){
  f <- file.path(fp, 'models', paste0(bn, '-brms-', family_label, '.rds'))
  s <- file.path(fp, 'models', paste0(bn, '-brms-', family_label, '-scaling.rds'))
  if(!file.exists(f) || !file.exists(s)){
    message('Missing brms fit/scaling, skipped: ', family_label)
    return(NULL)
  }
  fit <- readRDS(f)
  scaling <- readRDS(s)
  std_newdata <- apply_standardization(dplyr::select(df, -Site, -Phase, -point_id),
                                       scaling$center, scaling$scale)
  # The approximate Hilbert-space gp() term (density_bayes.R) extrapolates via
  # sinusoidal basis functions only well-behaved inside the training domain's
  # boundary factor; plots outside it blow up (to ~1e308 per draw, sometimes
  # literal Inf after averaging). A genuine extrapolation failure, not a coding
  # bug - left in the table and flagged via MaxPredicted rather than silently
  # trimmed. On this test set the failure is localized to Antero.
  preds <- colMeans(brms::posterior_epred(fit, newdata = std_newdata, allow_new_levels = TRUE))

  label <- paste0('brms ', tools::toTitleCase(family_label))
  n_blown_up <- sum(!is.finite(preds))
  if(n_blown_up > 0){
    message(sprintf('%s: GP term non-finite at %d/%d plots (site(s): %s) - outside the training spatial domain.',
                    label, n_blown_up, length(preds),
                    paste(unique(df$Site[!is.finite(preds)]), collapse = ', ')))
  }
  pred_row(label, preds)
}

## Phase B brms: same per-phase routing as score_ml_phaseB. Each phase has
## its own scaling (fit on that phase's training rows), so standardization
## is applied per phase rather than once over df.
score_brms_phaseB <- function(family_label){
  out <- rep(NA_real_, nrow(df))
  any_hit <- FALSE
  for(ph in names(bn_phaseB)){
    f <- file.path(fp, 'models', paste0(bn_phaseB[[ph]], '-brms-', family_label, '.rds'))
    s <- file.path(fp, 'models', paste0(bn_phaseB[[ph]], '-brms-', family_label, '-scaling.rds'))
    if(!file.exists(f) || !file.exists(s)){
      message(sprintf('Phase B %s-elev brms %s missing, those plots left NA.', ph, family_label))
      next
    }
    any_hit <- TRUE
    fit <- readRDS(f); scaling <- readRDS(s)
    mask <- df$Phase == ph
    nd <- apply_standardization(dplyr::select(df[mask, ], -Site, -Phase, -point_id),
                                scaling$center, scaling$scale)
    out[mask] <- colMeans(brms::posterior_epred(fit, newdata = nd, allow_new_levels = TRUE))
  }
  if(!any_hit) return(NULL)
  label <- paste0('brms ', tools::toTitleCase(family_label))
  if(any(!is.finite(out))){
    message(sprintf('%s (Phase B): GP non-finite at %d plots (site(s): %s).',
                    label, sum(!is.finite(out)),
                    paste(unique(df$Site[!is.finite(out)]), collapse = ', ')))
  }
  pred_row(label, out, scored_by = 'Phase B (per-phase)')
}

brms_results <- dplyr::bind_rows(
  dplyr::bind_rows(lapply(brms_families, score_brms)),
  dplyr::bind_rows(lapply(brms_families, score_brms_phaseB))
)

## 5. Conformal interval coverage for the ML candidates ---------------------
## (census_uncertainty_roadmap.md Phase 1c) - the internal test/CV rows used to
## calibrate these (functions.R's densityModeller()) can't judge real coverage
## (circular), so this is where it is actually checked: real field plots, the
## same place the brms gp() term's extrapolation failure first surfaced.
## Coverage is reported overall and on occupied plots, because an interval that
## only covers the zeros is not useful for census work.

score_conformal <- function(model_suffix, label){
  f <- file.path(fp, 'models', paste0(bn, '-', model_suffix, '-conformal.rds'))
  if(!file.exists(f)) return(NULL)
  cf <- readRDS(f)
  alpha <- 1 - conformal_level

  split_pred <- as.numeric(stats::predict(cf$split$model, new_data = df)[[1]])
  split_q <- stats::quantile(cf$split$residuals, probs = c(alpha / 2, 1 - alpha / 2),
                             na.rm = TRUE, names = FALSE)

  # Full Barber/Candes/Ramdas/Tibshirani CV+ interval: every pooled
  # out-of-fold (fold model, residual) pair contributes its own order statistic
  # (mu_{-k(i)}(x0) +- R_i), not just one randomly-chosen fold - unlike
  # combine_census.R's predict_density_draw_conformal(), which only needs one
  # fold per pseudo-draw and builds the same pooled quantile via many Monte
  # Carlo draws instead. Cheap here since `df` and the fold count are small.
  fold_preds <- sapply(cf$cv_plus$fold_models, function(m) as.numeric(stats::predict(m, new_data = df)[[1]]))
  mu_minus_k <- fold_preds[, cf$cv_plus$fold_id, drop = FALSE]  # n_plots x n_pooled_residuals
  lower_mat <- sweep(mu_minus_k, 2, cf$cv_plus$residuals, '-')
  upper_mat <- sweep(mu_minus_k, 2, cf$cv_plus$residuals, '+')
  cvplus_lower <- apply(lower_mat, 1, stats::quantile, probs = alpha / 2, na.rm = TRUE)
  cvplus_upper <- apply(upper_mat, 1, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE)

  dplyr::bind_rows(
    data.frame(Model = label, Method = 'split', point_id = df$point_id, Site = df$Site,
               Phase = df$Phase, Observed = df$Observed,
               lower = pmax(split_pred + split_q[1], 0), upper = split_pred + split_q[2]),
    data.frame(Model = label, Method = 'cv_plus', point_id = df$point_id, Site = df$Site,
               Phase = df$Phase, Observed = df$Observed,
               lower = pmax(cvplus_lower, 0), upper = cvplus_upper)
  )
}

conformal_results <- dplyr::bind_rows(
  score_conformal('poisson_spat', 'XGB Poisson Spat.'),
  score_conformal('poisson', 'XGB Poisson'),
  score_conformal('tweedie_spat', 'XGB Tweedie Spat.'),
  score_conformal('tweedie', 'XGB Tweedie'),
  score_conformal('lgbm-poisson_spat', 'LGBM Poisson Spat.')
)

if(!is.null(conformal_results) && nrow(conformal_results) > 0){
  conformal_summary <- conformal_results |>
    dplyr::mutate(covered = Observed >= lower & Observed <= upper, Occupied = Observed > 0) |>
    dplyr::group_by(Model, Method) |>
    dplyr::summarize(
      nominal_level    = conformal_level,
      coverage_pct     = mean(covered) * 100,
      coverage_pct_occ = mean(covered[Occupied]) * 100,
      mean_width       = mean(upper - lower),
      mean_width_occ   = mean((upper - lower)[Occupied]),
      n                = dplyr::n(),
      n_occ            = sum(Occupied),
      .groups = 'drop'
    ) |>
    dplyr::arrange(Model, Method)

  cat('\n=== conformal interval coverage/width vs. 2026 field groundtruth (nominal ',
      conformal_level * 100, '%) ===\n', sep = '')
  print(as.data.frame(conformal_summary))
}

## 6. Summarize -------------------------------------------------------------

## Tier labels so the table reads as a ladder.
tier_of <- function(model){
  dplyr::case_when(
    grepl('^Constant', model)                    ~ '1 Constant',
    grepl('^GLM', model)                         ~ '2 Single covariate',
    grepl('^IDW|^Ordinary kriging|^kNN', model)  ~ '3 Spatial interpolation',
    TRUE                                         ~ '4 Fitted model'
  )
}

all_preds <- dplyr::bind_rows(baseline_results, ml_results, brms_results) |>
  dplyr::mutate(Tier = tier_of(Model),
                Occupied = Observed > 0,
                Error = Predicted - Observed,
                AbsError = abs(Predicted - Observed),
                # Only ever read on occupied plots (Observed >= 1), so the
                # division is safe there; Inf on true zeros by construction.
                PctError = 100 * abs(Predicted - Observed) / Observed)

# All statistics drop non-finite predictions (brms GP blow-ups) so a model with
# a few failed plots still gets a readable score on the rest; n_dropped keeps
# that visible.
fin <- function(x) x[is.finite(x)]

#' Subset a model's prediction rows, then reduce them with `f`.
#' @param phase 'High'/'Low', or NULL for both.
#' @param occupied_only TRUE restricts to plots where the plant was actually
#' found (Observed > 0).
stat_for <- function(d, f, col = 'AbsError', phase = NULL, occupied_only = FALSE){
  if(!is.null(phase)) d <- d[d$Phase == phase, ]
  if(occupied_only)  d <- d[d$Occupied, ]
  v <- fin(d[[col]])
  if(length(v) == 0) return(NA_real_)
  f(v)
}

summary_tbl <- all_preds |>
  dplyr::group_by(ScoredBy, Tier, Model) |>
  dplyr::group_modify(~ data.frame(
    # n_plots counts plots the row actually tried to score (Phase B's per-
    # phase rows are NA on the other phase by design, not by failure).
    # n_dropped counts only true blow-ups among those (non-finite preds).
    n_plots            = sum(!is.na(.x$Predicted)),
    n_dropped          = sum(!is.na(.x$Predicted) & !is.finite(.x$Predicted)),
    # all plots - reported, but see the header: for a rare species this cut
    # mostly rewards predicting zero
    MAE                = stat_for(.x, mean),
    MedianError        = stat_for(.x, stats::median, 'Error'),
    MedianAbsError     = stat_for(.x, stats::median),
    MAE_HighElev       = stat_for(.x, mean, phase = 'High'),
    MAE_LowElev        = stat_for(.x, mean, phase = 'Low'),
    # occupied plots only - the conditional "how many are there?" question
    MAE_occ            = stat_for(.x, mean, occupied_only = TRUE),
    MedianError_occ    = stat_for(.x, stats::median, 'Error', occupied_only = TRUE),
    MedianAbsError_occ = stat_for(.x, stats::median, occupied_only = TRUE),
    MAE_occ_HighElev   = stat_for(.x, mean, phase = 'High', occupied_only = TRUE),
    MAE_occ_LowElev    = stat_for(.x, mean, phase = 'Low',  occupied_only = TRUE),
    # MdAE per phase as well as MAE: the occupied low-elevation counts run to
    # 171 plants, so even within that phase MAE is pulled by its own tail.
    # MdAE is the typical-plot error and is what shows whether a model is
    # merely imprecise or systematically on the wrong scale for that phase.
    MdAE_occ_HighElev  = stat_for(.x, stats::median, phase = 'High', occupied_only = TRUE),
    MdAE_occ_LowElev   = stat_for(.x, stats::median, phase = 'Low',  occupied_only = TRUE),
    # scale-free, so the two phases are comparable despite 10x different counts
    MdAPE_occ          = stat_for(.x, stats::median, 'PctError', occupied_only = TRUE),
    MdAPE_occ_HighElev = stat_for(.x, stats::median, 'PctError', phase = 'High', occupied_only = TRUE),
    MdAPE_occ_LowElev  = stat_for(.x, stats::median, 'PctError', phase = 'Low',  occupied_only = TRUE),
    # The gp() blow-up is not always non-finite: with Upper Crossing dropped the
    # worst-extrapolating plots are gone and the rest come back huge-but-finite
    # (~1e75), which is.finite() won't catch but which still makes a mean
    # meaningless. MaxPredicted keeps that visible instead of hiding it inside
    # an aggregate. No threshold is applied - the MAE reported is the real one.
    MaxPredicted       = max(fin(.x$Predicted)),
    n_occ              = sum(.x$Occupied),
    n_HighElev         = sum(.x$Phase == 'High'),
    n_LowElev          = sum(.x$Phase == 'Low'),
    n_occ_HighElev     = sum(.x$Occupied & .x$Phase == 'High'),
    n_occ_LowElev      = sum(.x$Occupied & .x$Phase == 'Low')
  )) |>
  dplyr::ungroup() |>
  # Sort by MAE_occ_LowElev: this is the phase and cut where models have room to
  # earn their complexity (low-elev occupied counts run to 171, so there is
  # magnitude to get wrong; high-elev occupied counts are low single digits, so
  # the constants are already in-band and the ranking there is noise-dominated).
  dplyr::arrange(ScoredBy, Tier, MAE_occ_LowElev) |>
  dplyr::mutate(dplyr::across(c(MAE, MedianError, MedianAbsError, MAE_HighElev, MAE_LowElev,
                                MAE_occ, MedianError_occ, MedianAbsError_occ,
                                MAE_occ_HighElev, MAE_occ_LowElev, MdAE_occ_HighElev, MdAE_occ_LowElev,
                                MdAPE_occ, MdAPE_occ_HighElev, MdAPE_occ_LowElev),
                              ~ round(.x, 3)),
                MaxPredicted = signif(MaxPredicted, 4))

## 6b. KNOWN vs NEW population regimes -------------------------------------
## The operational question is not "which model is best" but "which model for
## which situation": estimating a total at a population you already have plots
## for is a different problem from guessing what a newly-reported population
## holds. The split uses the VERIFIED crosswalk above (identity), not distance -
## `dist_m` is retained only as context, and the two disagree materially (see
## the crosswalk note). No dead band is needed with identity.
regime_preds <- all_preds |>
  dplyr::left_join(site_dist, by = 'Site') |>
  dplyr::mutate(
    train_pop = unname(population_crosswalk[Site]),
    Regime    = ifelse(is.na(train_pop), 'NEW (fresh survey)', 'KNOWN (resurvey)')
  )

regime_tbl <- regime_preds |>
  dplyr::group_by(ScoredBy, Regime, Tier, Model) |>
  dplyr::group_modify(~ data.frame(
    n_occ              = sum(.x$Occupied & !is.na(.x$Predicted)),
    n_occ_HighElev     = sum(.x$Occupied & .x$Phase == 'High' & !is.na(.x$Predicted)),
    n_occ_LowElev      = sum(.x$Occupied & .x$Phase == 'Low'  & !is.na(.x$Predicted)),
    MAE_occ            = stat_for(.x, mean, occupied_only = TRUE),
    MAE_occ_HighElev   = stat_for(.x, mean, phase = 'High', occupied_only = TRUE),
    MAE_occ_LowElev    = stat_for(.x, mean, phase = 'Low',  occupied_only = TRUE),
    MdAE_occ_HighElev  = stat_for(.x, stats::median, phase = 'High', occupied_only = TRUE),
    MdAE_occ_LowElev   = stat_for(.x, stats::median, phase = 'Low',  occupied_only = TRUE)
  )) |>
  dplyr::ungroup() |>
  dplyr::arrange(ScoredBy, Regime, MAE_occ_LowElev) |>
  dplyr::mutate(dplyr::across(c(MAE_occ, MAE_occ_HighElev, MAE_occ_LowElev,
                                MdAE_occ_HighElev, MdAE_occ_LowElev), ~ round(.x, 3)))

cat('\n--- KNOWN vs NEW population regimes, occupied plots only ---\n')
cat('Regime membership (VERIFIED crosswalk; dist_m/nearest_train_pop are context only and disagree):\n')
print(as.data.frame(dplyr::distinct(regime_preds, Regime, Site, train_pop, dist_m, nearest_train_pop) |>
                      dplyr::arrange(Regime, dist_m)), row.names = FALSE)
for(sb in unique(regime_tbl$ScoredBy)){
  for(g in unique(regime_tbl$Regime)){
    cat('\n== ', sb, ' / ', g, ' ==\n', sep = '')
    print(as.data.frame(dplyr::select(dplyr::filter(regime_tbl, ScoredBy == sb, Regime == g),
                                      -Regime, -ScoredBy)), row.names = FALSE)
  }
}

blowups <- dplyr::filter(summary_tbl, n_dropped > 0 | MaxPredicted > 1e4)
if(nrow(blowups) > 0){
  message('GP extrapolation failure - MAE is not interpretable for: ',
          paste(blowups$Model, collapse = ', '),
          ' (max predicted ', paste(signif(blowups$MaxPredicted, 3), collapse = ', '), ').')
}

obs <- dplyr::distinct(all_preds, point_id, .keep_all = TRUE)
cat('\n=== 2026 field groundtruth (', nrow(obs), ' plots: ',
    sum(obs$Phase == 'High'), ' high-elevation, ', sum(obs$Phase == 'Low'), ' low-elevation; ',
    sum(obs$Occupied), ' occupied). Excluded: ', paste(names(excluded_sites), collapse = ', '), ' ===\n',
    'Observed counts, ALL plots:      mean ', round(mean(obs$Observed), 2), ' / median ', stats::median(obs$Observed),
    ' (high ', round(mean(obs$Observed[obs$Phase == 'High']), 2), '/', stats::median(obs$Observed[obs$Phase == 'High']),
    ', low ', round(mean(obs$Observed[obs$Phase == 'Low']), 2), '/', stats::median(obs$Observed[obs$Phase == 'Low']), ')\n',
    'Observed counts, OCCUPIED plots: mean ', round(mean(obs$Observed[obs$Occupied]), 2), ' / median ', stats::median(obs$Observed[obs$Occupied]),
    ' (high ', round(mean(obs$Observed[obs$Occupied & obs$Phase == 'High']), 2), '/', stats::median(obs$Observed[obs$Occupied & obs$Phase == 'High']),
    ', low ', round(mean(obs$Observed[obs$Occupied & obs$Phase == 'Low']), 2), '/', stats::median(obs$Observed[obs$Occupied & obs$Phase == 'Low']), ')\n', sep = '')

## The headline cut is PHASE-STRATIFIED: the two elevational phases are on
## totally different count scales (occupied median 6 high vs 60.5 low), the
## final models will be fit separately per phase, and pooling them hides which
## models earn their complexity on which phase. Pooled MAE/MAE_occ are retained
## in the CSV for context but demoted here.
for(sb in unique(summary_tbl$ScoredBy)){
  cat('\n########## ', sb, ' ##########\n', sep = '')
  sub <- dplyr::filter(summary_tbl, ScoredBy == sb)
  cat('\n--- ALL plots by phase (rare species: constants are near-optimal on the all-plot cut by construction) ---\n')
  print(as.data.frame(sub[, c('Tier', 'Model', 'MAE_HighElev', 'MAE_LowElev',
                              'MAE', 'MedianError', 'MedianAbsError')]))
  cat('\n--- OCCUPIED plots by phase: given presence, how many? (the cut that judges a count model) ---\n')
  cat('    Final models will be fit per phase, so the per-phase columns are the headline.\n')
  print(as.data.frame(sub[, c('Tier', 'Model',
                              'MAE_occ_HighElev', 'MAE_occ_LowElev',
                              'MdAE_occ_HighElev', 'MdAE_occ_LowElev',
                              'MdAPE_occ_HighElev', 'MdAPE_occ_LowElev',
                              'MAE_occ', 'MedianError_occ', 'MedianAbsError_occ', 'MdAPE_occ')]))
}

## Phase A vs Phase B direct comparison on the shared fitted stack.
## Negative delta = Phase B wins. Baselines are not paired here (per-phase
## constants/GLMs have different model names than their pooled counterparts).
has_phaseB_fits <- any(summary_tbl$ScoredBy == 'Phase B (per-phase)' &
                       summary_tbl$Tier    == '4 Fitted model')
if(has_phaseB_fits){
  ab_cmp <- summary_tbl |>
    dplyr::filter(Tier == '4 Fitted model') |>
    dplyr::select(Model, ScoredBy, MAE_occ_HighElev, MAE_occ_LowElev,
                  MdAE_occ_HighElev, MdAE_occ_LowElev) |>
    tidyr::pivot_wider(names_from = ScoredBy,
                       values_from = c(MAE_occ_HighElev, MAE_occ_LowElev,
                                       MdAE_occ_HighElev, MdAE_occ_LowElev)) |>
    dplyr::mutate(
      delta_MAE_HighElev = `MAE_occ_HighElev_Phase B (per-phase)` - `MAE_occ_HighElev_Phase A (pooled)`,
      delta_MAE_LowElev  = `MAE_occ_LowElev_Phase B (per-phase)`  - `MAE_occ_LowElev_Phase A (pooled)`
    ) |>
    dplyr::arrange(delta_MAE_LowElev)

  cat('\n--- Phase A vs Phase B on the fitted stack (negative delta = Phase B wins) ---\n')
  print(as.data.frame(ab_cmp[, c('Model',
                                 'delta_MAE_HighElev', 'delta_MAE_LowElev',
                                 'MAE_occ_HighElev_Phase A (pooled)',
                                 'MAE_occ_HighElev_Phase B (per-phase)',
                                 'MAE_occ_LowElev_Phase A (pooled)',
                                 'MAE_occ_LowElev_Phase B (per-phase)')]))
} else {
  ab_cmp <- NULL
  message('Phase B models not present yet - run First_modelling/abundance_phaseB_run.R first for the A/B delta table.')
}

## 7. Write ----------------------------------------------------------------

dir.create(file.path('..', 'results', 'tables'), showWarnings = FALSE)
out_csv <- file.path('..', 'results', 'tables', paste0(bn, '-2026groundtruth-validation.csv'))
write.csv(summary_tbl, out_csv, row.names = FALSE)

out_regime <- file.path('..', 'results', 'tables', paste0(bn, '-2026groundtruth-regime.csv'))
write.csv(regime_tbl, out_regime, row.names = FALSE)

if(!is.null(ab_cmp)){
  out_ab <- file.path('..', 'results', 'tables', paste0(bn, '-2026groundtruth-AvsB.csv'))
  write.csv(ab_cmp, out_ab, row.names = FALSE)
  message('A/B delta: ', out_ab)
}

# Plot-level predictions carry dist_m/nearest_train_pop/Regime so the regime cut
# is reproducible downstream without recomputing distances.
out_preds <- file.path('..', 'results', 'tables', paste0(bn, '-2026groundtruth-predictions.csv'))
write.csv(dplyr::left_join(all_preds, site_dist, by = 'Site') |>
            dplyr::mutate(train_pop = unname(population_crosswalk[Site]),
                          Regime = ifelse(is.na(train_pop), 'NEW (fresh survey)', 'KNOWN (resurvey)')),
          out_preds, row.names = FALSE)

if(exists('conformal_summary')){
  out_conf <- file.path('..', 'results', 'tables', paste0(bn, '-2026groundtruth-conformal.csv'))
  write.csv(conformal_summary, out_conf, row.names = FALSE)
  message('Conformal coverage: ', out_conf)
}

message('Written to: ', out_csv)
message('Regime split: ', out_regime)
message('Plot-level predictions: ', out_preds)
