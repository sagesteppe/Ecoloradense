## Phase 4 of census_uncertainty_roadmap.md - Validate before trusting it.

#' Reshape/indexing correctness + reproducibility check for `predict_density_draw()`.
#'
#' @description Two checks, both needed because of a behavior confirmed by
#' testing: for a model with a `gp()` term, `posterior_epred()` at new
#' (out-of-sample) locations draws the GP's latent value from R's global RNG
#' on every call, and that random draw is a *joint* function of the whole
#' batch of new locations in the call - predicting the same cells as part of
#' a differently-sized batch does not reproduce the same values even with
#' the same seed, because the batch itself, not just the seed, determines the
#' draw. That rules out the obvious check ("does a small subset's prediction
#' match its value inside the full-raster prediction?" - it never will, by
#' design, not because of a bug) in favor of two checks that are actually
#' meaningful:
#'
#' 1. **Reproducibility**: two full-raster calls to `predict_density_draw()`
#' with the same `draw_id` must give identical output - required for Phase
#' 3's Monte Carlo composition, since `pair_draws()` resamples density draws
#' with replacement and a reused `draw_id` must give back the identical
#' surface, not a fresh random one.
#' 2. **Reshape/indexing**: `predict_density_draw()`'s full-raster output
#' must match an independent, same-seed, same-full-`newdata`
#' `posterior_epred()` call, read back at the same cells - this is the
#' actual apples-to-apples comparison (same batch, same seed), and would
#' catch a cell-indexing/reshaping bug in `predict_density_draw()`.
#'
#' @param brms_fit,raster_template as `predict_density_draw()`.
#' @param draw_id which posterior draw to check.
#' @param n_check how many cells to spot-check (sampled from the non-`NA` cells).
#' @param cell_sample_seed RNG seed for *which* cells get spot-checked.
#' @return invisible TRUE if both checks pass; stops with an informative
#' message otherwise.
check_predict_density_draw <- function(brms_fit, raster_template, draw_id = 1,
                                        n_check = 50, cell_sample_seed = 1){

  covariate_df <- as.data.frame(raster_template, cells = TRUE)

  pred_a <- predict_density_draw(brms_fit, raster_template, draw_id, covariate_df = covariate_df)
  pred_b <- predict_density_draw(brms_fit, raster_template, draw_id, covariate_df = covariate_df)
  repro_mismatch <- max(abs(terra::values(pred_a) - terra::values(pred_b)), na.rm = TRUE)
  if(repro_mismatch > 1e-8){
    stop("check_predict_density_draw(): two calls with the same draw_id = ", draw_id,
         " did not reproduce the same surface (max diff = ", repro_mismatch, ") - ",
         "predict_density_draw()'s seeding is not making gp() prediction reproducible, ",
         "which Phase 3's with-replacement resampling of density draws depends on.")
  }

  set.seed(draw_id)
  reference <- as.vector(brms::posterior_epred(brms_fit, newdata = covariate_df, draw_ids = draw_id))

  set.seed(cell_sample_seed)
  check_cells <- sample(covariate_df$cell, min(n_check, nrow(covariate_df)))
  # index by cell number directly (r[cells]), the same convention
  # predict_density_draw() itself uses to *write* values (out[chunk$cell] <- ...) -
  # terra::extract() expects spatial geometries/coordinates, not bare cell
  # indices, and silently misinterprets them as such if given directly.
  from_raster <- pred_a[check_cells][, 1]
  from_reference <- reference[match(check_cells, covariate_df$cell)]

  mismatch <- abs(from_reference - from_raster)
  if(any(mismatch > 1e-8 * pmax(abs(from_reference), 1))){
    stop("check_predict_density_draw(): predict_density_draw()'s output does not ",
         "match a same-seed, same-full-newdata posterior_epred() call read back at ",
         "the same cells - largest mismatch = ", max(mismatch), ". This points to a ",
         "cell-indexing/reshaping bug in predict_density_draw().")
  }
  invisible(TRUE)
}

#' Sanity-check the boundary ensemble's median mask against the two
#' deterministic threshold rules.
#'
#' @description Per-cell inclusion frequency across the Phase 2 mask
#' ensemble, thresholded at 0.5 for a median mask; compares its area against
#' the two deterministic `spec_sens`/`sensitivity` masks' areas (from
#' `scripts/Second_groundtruth/make_threshold_raster_masks.R`) and flags if
#' the ensemble median falls outside the range the two deterministic rules
#' bracket.
#'
#' @param boundary_mask_paths character vector, Phase 2's mask file paths.
#' @param spec_sens_mask_path,sensitivity_mask_path paths to the deterministic masks.
#' @return list(ensemble_median_area, spec_sens_area, sensitivity_area, in_bracket = logical).
sanity_check_boundary_median <- function(boundary_mask_paths, spec_sens_mask_path,
                                          sensitivity_mask_path){

  masks <- terra::rast(boundary_mask_paths)
  inclusion_freq <- terra::app(masks, fun = function(x) mean(!is.na(x)))
  median_mask <- terra::ifel(inclusion_freq >= 0.5, 1, NA)

  cell_area <- prod(terra::res(median_mask))
  area_of <- function(m) terra::global(!is.na(m), fun = 'sum', na.rm = TRUE)[1, 1] * cell_area

  ensemble_median_area <- area_of(median_mask)
  spec_sens_area <- area_of(terra::rast(spec_sens_mask_path))
  sensitivity_area <- area_of(terra::rast(sensitivity_mask_path))

  lo <- min(spec_sens_area, sensitivity_area)
  hi <- max(spec_sens_area, sensitivity_area)

  list(ensemble_median_area = ensemble_median_area,
       spec_sens_area = spec_sens_area, sensitivity_area = sensitivity_area,
       in_bracket = ensemble_median_area >= lo && ensemble_median_area <= hi)
}

#' Re-check the suitability-as-covariate vs. spatial-random-effect spatial
#' confounding question, now that the density model's `gp()` term is
#' actually fit (roadmap Phase 4: "not just in principle").
#'
#' @description Compares fixed-effect coefficient estimates (median +
#' credible interval) between the promoted model (covariates + `gp()`) and
#' the same formula refit without the `gp()` term, on the same training
#' data. Large shifts in a covariate's estimate (especially `Pr.SuitHab`,
#' the habitat-suitability covariate most likely to be spatially confounded
#' with the `gp()` term) flag the confounding this check exists to catch -
#' this is a diagnostic only, not new fitting infrastructure: it reuses
#' `brms_count_model()`'s formula builder and refits once, cheaply, without
#' the spatial term.
#'
#' @param promoted_fit the `brms_promote_and_refit()$Model` (covariates + `gp()`).
#' @param train the same training data.frame it was fit on.
#' @param family the same brms family used for the promoted fit.
#' @param shift_tol relative-change threshold on a coefficient's median to flag.
#' @param backend,chains,iter,warmup passed to the no-`gp()` refit.
#' @return data.frame(term, with_gp_median, without_gp_median, rel_shift, flagged).
check_spatial_confounding <- function(promoted_fit, train, family, shift_tol = 0.5,
                                       backend = 'cmdstanr', chains = 2, iter = 1000, warmup = 500){

  covars <- setdiff(names(train), c('Prsnc_All', 'Longitude', 'Latitude'))
  form_no_gp <- as.formula(paste('Prsnc_All ~', paste(covars, collapse = ' + ')))

  fit_no_gp <- brms::brm(form_no_gp, data = train, family = family, backend = backend,
                          chains = chains, iter = iter, warmup = warmup, silent = 2, refresh = 0)

  with_gp <- brms::fixef(promoted_fit)[covars, 'Estimate']
  without_gp <- brms::fixef(fit_no_gp)[covars, 'Estimate']
  rel_shift <- abs(with_gp - without_gp) / pmax(abs(without_gp), 1e-8)

  data.frame(term = covars, with_gp_median = with_gp, without_gp_median = without_gp,
             rel_shift = rel_shift, flagged = rel_shift > shift_tol)
}
