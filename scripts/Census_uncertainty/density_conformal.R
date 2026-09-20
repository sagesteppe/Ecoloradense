## Phase 1c of census_uncertainty_roadmap.md - conformal intervals for the ML
## density candidates (XGB Poisson/Tweedie, spatial + non-spatial, LGBM
## Poisson Spat.). These candidates are point-estimate only in `densityModeller()`
## (../functions.R) - unlike the brms candidates (density_bayes.R), which get a
## real posterior. Source ../functions.R before this file - `conformal_calibrate_candidate()`
## is called from `densityModeller()`, and `predict_density_draw_conformal()` is
## called from `combine_census.R`'s `combined_census_montecarlo_conformal()`.

#' Calibrate split-conformal and spatial CV+ residuals for one fitted ML candidate.
#'
#' @description Two calibrations, compared rather than committing to one (same
#' spirit as Phase 1 comparing 4 brms families up front):
#' - **Split**: `fitted` predicted on `test` - the same held-out split every ML
#'   candidate (`poiss()`/`tweed()`/`gbs()`, ../functions.R) is already scored
#'   on, so this reuses an existing split rather than carving out a fresh
#'   calibration set. Residuals are kept **signed** (`Observed - Predicted`),
#'   not absolute - this project's count-data residuals are right-skewed, so
#'   the two interval tails need independent quantiles rather than a single
#'   symmetric `+-` width.
#' - **Spatial CV+** (Barber, Candes, Ramdas & Tibshirani 2021): `cv_rset` is
#'   `densityModeller()`'s existing `indx_nndm_rs` (`CAST2rsample()` on the
#'   fixed-k `CAST::knndm()` folds already built for XGB/LGBM tuning - no new
#'   spatial-fold computation here). Refits once per fold at `fitted$spec`'s
#'   already-chosen hyperparameters (`parsnip::fit(fitted$spec, ...)`, not a
#'   re-tune - same cheap-refit-at-fixed-hyperparameters idea as
#'   `brms_promote_and_refit()`), predicts the held-out fold, and pools every
#'   fold's out-of-fold residual together with which fold model produced it -
#'   the quantities the Barber et al. CV+ interval formula needs at prediction
#'   time (see `predict_density_draw_conformal()`).
#'
#' @param fitted a `parsnip` `model_fit` - `poiss()`/`tweed()`/`gbs()`'s
#' `$Model` (../functions.R). Must carry `$spec` (the finalized model spec
#' `fit()` was called with), which every plain `parsnip::fit()` result does.
#' @param train,test data.frames as `densityModeller()` uses - `test` is only
#' used for the split calibration; `train`'s own rows are only reached via
#' `cv_rset`'s analysis/assessment splits.
#' @param cv_rset an `rsample` rset (e.g. `CAST2rsample()`'s output) whose
#' splits' `analysis()`/`assessment()` subset `train`.
#' @param response column name of the observed count.
#' @return `list(split = list(model = fitted, residuals = <numeric>), cv_plus =
#' list(fold_models = <list of K model_fits>, fold_id = <integer, one per
#' pooled residual, indexing into fold_models>, residuals = <numeric, pooled
#' across folds>))`.
conformal_calibrate_candidate <- function(fitted, train, test, cv_rset, response = 'Prsnc_All'){

  split_pred <- as.numeric(stats::predict(fitted, new_data = test)[[1]])
  split_residuals <- test[[response]] - split_pred

  form <- stats::as.formula(paste(response, '~ .'))

  fold_results <- lapply(cv_rset$splits, function(sp){
    analysis_set   <- rsample::analysis(sp)
    assessment_set <- rsample::assessment(sp)
    fold_fit  <- parsnip::fit(fitted$spec, form, data = analysis_set)
    fold_pred <- as.numeric(stats::predict(fold_fit, new_data = assessment_set)[[1]])
    list(model = fold_fit, residuals = assessment_set[[response]] - fold_pred)
  })

  fold_models <- lapply(fold_results, `[[`, 'model')
  fold_id <- rep(seq_along(fold_results),
                  vapply(fold_results, function(r) length(r$residuals), integer(1)))
  pooled_residuals <- unlist(lapply(fold_results, `[[`, 'residuals'), use.names = FALSE)

  list(
    split   = list(model = fitted, residuals = split_residuals),
    cv_plus = list(fold_models = fold_models, fold_id = fold_id, residuals = pooled_residuals)
  )
}

#' Width of a residual-quantile interval, for the `densityModeller()` comparison table.
#'
#' @description Informational only, **not** a coverage estimate - computing
#' coverage against the same `test`/`cv_rset` rows the residuals were
#' calibrated on would be circular. Real coverage is checked where this
#' project already checks real generalization:
#' `validate_against_2026_groundtruth.R`, against genuinely new field plots.
#'
#' @param residuals numeric vector of signed residuals (`conformal_calibrate_candidate()`'s
#' `$split$residuals` or `$cv_plus$residuals`).
#' @param level nominal interval level (default 0.9, i.e. a 5%/95% residual split).
#' @return numeric scalar, `upper_quantile - lower_quantile`.
conformal_interval_width <- function(residuals, level = 0.9){
  alpha <- 1 - level
  q <- stats::quantile(residuals, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE, names = FALSE)
  unname(diff(q))
}

#' Predict one conformal pseudo-draw's density surface onto a raster template.
#'
#' @description The Phase 3-facing counterpart to `combine_census.R`'s
#' `predict_density_draw()` (brms posterior draws), for the ML candidates'
#' conformal calibration instead. Returns the same shape (single-layer
#' `terra::rast`) so `combined_census_montecarlo_conformal()` can loop it the
#' same "predict -> use -> discard" way Phase 3 already does for brms.
#'
#' **One resampled residual for the whole raster, not one per cell.** An
#' independently-resampled residual at every cell would mostly cancel out
#' under `terra::global(fun='sum')` (same reasoning Phase 2 already applied to
#' the boundary field: per-cell iid noise vs. a spatially coherent
#' realization) - it would understate this candidate's real census-size
#' uncertainty. Each pseudo-draw instead applies one residual as a flat shift
#' across the entire surface, which is also the more faithful reading of what
#' a conformal interval actually licenses (a bound on one new point's error,
#' not independent per-cell error).
#'
#' `method = 'cv_plus'` samples one `(fold, residual)` pair from the pooled
#' out-of-fold pool - using that fold's own held-out model to predict the
#' whole raster (so different pseudo-draws use genuinely different fitted
#' models, not just different residual shifts) plus that same fold's own
#' residual, matching Barber et al.'s CV+ pairing (`mu_{-k(i)}(x0) +- R_i`).
#' Unlike a single evaluation of the CV+ formula (all K fold models' order
#' statistics at once - see `validate_against_2026_groundtruth.R`'s
#' `score_conformal()`), this only needs one fold per draw; enough draws
#' approximate the same pooled quantile via Monte Carlo instead.
#'
#' @param conformal one candidate's `conformal_calibrate_candidate()` output.
#' @param method 'split' or 'cv_plus'.
#' @param raster_template a `terra::rast` covering the prediction domain, one
#' layer per covariate the model needs, named to match (same contract as
#' `predict_density_draw()`).
#' @param draw_id integer, this pseudo-draw's id (only used as the default seed).
#' @param covariate_df optional pre-extracted `as.data.frame(raster_template, cells=TRUE)` -
#' pass this in when calling repeatedly so it's built once, not per draw.
#' @param seed RNG seed for this draw's residual/fold resampling; defaults to
#' `draw_id` so the same draw is always reproducible.
#' @return a `terra::rast` (single layer) of predicted density for this pseudo-draw.
predict_density_draw_conformal <- function(conformal, method = c('split', 'cv_plus'),
                                            raster_template, draw_id,
                                            covariate_df = NULL, seed = draw_id){
  method <- match.arg(method)
  if(is.null(covariate_df)){
    covariate_df <- as.data.frame(raster_template, cells = TRUE)
  }

  set.seed(seed)
  if(method == 'split'){
    point_pred <- as.numeric(stats::predict(conformal$split$model, new_data = covariate_df)[[1]])
    resid_draw <- sample(conformal$split$residuals, 1)
  } else {
    i <- sample.int(length(conformal$cv_plus$residuals), 1)
    fold_model <- conformal$cv_plus$fold_models[[conformal$cv_plus$fold_id[i]]]
    point_pred <- as.numeric(stats::predict(fold_model, new_data = covariate_df)[[1]])
    resid_draw <- conformal$cv_plus$residuals[i]
  }

  out <- terra::rast(raster_template, nlyrs = 1)
  out[covariate_df$cell] <- pmax(point_pred + resid_draw, 0)  # counts can't go negative
  out
}
