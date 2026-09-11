## Phase 1 of census_uncertainty_roadmap.md - Bayesian density candidates.
##
## Fits brms count-model families (poisson, negbinomial, hurdle_poisson,
## hurdle_negbinomial) as additional candidates in the same comparison table
## `densityModeller()` (../functions.R) already builds for the five ML
## candidates (XGBoost Poisson/Tweedie spatial+non-spatial, LightGBM Poisson
## spatial). Source ../functions.R before this file - `brms_cv_compare()` and
## `brms_promote_and_refit()` are called from `densityModeller()`.

#' Z-score covariates before fitting a brms model.
#'
#' @description Stan's default random inits are drawn uniformly on (-2, 2) on
#' the unconstrained parameter scale; combined with this project's raw
#' environmental covariates spanning wildly different natural units -
#' `elevation` in the thousands, `NDVI` near 0, raw UTM coordinates in the
#' millions - that blows up the initial linear predictor enough to send the
#' count-family log-likelihood to `log(0)`, and Stan can't start sampling at
#' all ("Initialization between (-2, 2) failed after 100 attempts", confirmed
#' via `Census_uncertainty/diagnose_brms_failure.R`). Same numerical failure
#' mode `sagesteppe/safeHavens`'s `bayesianSDM()` already hit and fixed -
#' matches its approach here: environmental covariates are z-scored, and the
#' `gp()` coordinates are centred and rescaled to **kilometres** rather than
#' left in raw UTM metres (`safeHavens::add_planar_coords()`'s
#' `scale_km = TRUE`: "avoids the numerical issues that arise when GP
#' length-scales must span millions of metres") - full z-scoring isn't used
#' for the coordinates themselves so the GP lengthscale stays interpretable in
#' (scaled) spatial units rather than an arbitrary unitless one.
#'
#' @param train data.frame to compute scaling parameters from.
#' @param newdata optional data.frame to apply the identical transform to (a
#' held-out test set at CV time, or a raster-derived prediction frame at
#' Phase 3 time) - must have the same covariate columns as `train`.
#' @param response,gp_coords columns excluded from scaling (`gp_coords` are
#' centred and divided by `gp_km_divisor`, not by their own sd).
#' @param gp_km_divisor divisor applied to centred `gp_coords` (1000 for
#' metre-unit CRS's, matching `safeHavens::add_planar_coords()`).
#' @return list(train, newdata, center, scale) - `center`/`scale` are named
#' by column, and are what a caller must reapply to standardize any further
#' new data (e.g. Phase 3's raster covariates) with this exact transform.
standardize_covariates <- function(train, newdata = NULL, response = 'Prsnc_All',
                                    gp_coords = c('Longitude', 'Latitude'),
                                    gp_km_divisor = 1000){

  covars <- setdiff(names(train), response)
  center <- vapply(train[covars], mean, numeric(1), na.rm = TRUE)
  scale  <- vapply(train[covars], stats::sd, numeric(1), na.rm = TRUE)
  scale[gp_coords[gp_coords %in% covars]] <- gp_km_divisor
  scale[!is.finite(scale) | scale == 0] <- 1  # guard constant/degenerate columns

  rescale <- function(d){
    for(v in covars) d[[v]] <- (d[[v]] - center[[v]]) / scale[[v]]
    d
  }

  list(train = rescale(train), newdata = if(is.null(newdata)) NULL else rescale(newdata),
       center = center, scale = scale)
}

#' Apply a previously-fit standardize_covariates() transform to new data.
#'
#' @description The Phase 3/validation counterpart to `standardize_covariates()`:
#' that function both fits *and* applies a transform, but raster-scale
#' prediction (`combine_census.R`) and the Phase 4 checks (`validate_pipeline.R`)
#' need to apply an *already-fit* model's saved `center`/`scale` (from
#' `brms_promote_and_refit()`'s `-scaling.rds`) to brand-new data instead.
#'
#' @param newdata data.frame (or a raster-derived data.frame) with the same
#' covariate columns `center`/`scale` were computed on.
#' @param center,scale named numeric vectors from `standardize_covariates()`.
#' @return `newdata` with every named column in `center`/`scale` rescaled.
apply_standardization <- function(newdata, center, scale){
  for(v in names(center)) if(v %in% names(newdata)) newdata[[v]] <- (newdata[[v]] - center[[v]]) / scale[[v]]
  newdata
}

#' Build the shared brms formula: response ~ covariates + approximate gp(x, y).
#'
#' @description Uses brms's approximate (Hilbert-space) `gp()` term
#' (`gp_k`/`gp_c`) rather than an exact GP, so the fitted posterior stays a
#' plain coefficient+basis-weight matrix (see `brms_promote_and_refit()`) and
#' stays fittable in Stan at plot-level n.
#' @param train data.frame the formula's covariates are read off of (every
#' column except `Prsnc_All` and `gp_coords`).
#' @param gp_coords character(2), coordinate column names for the `gp()` term.
#' @param gp_k,gp_c approximate-GP basis dimension / boundary factor.
#' **`gp_k` must be a number, not `NA`**: per `?brms::gp`, `k = NA` (which
#' looks like "let brms pick a default") actually means "compute an *exact*
#' GP" - confirmed by testing, this was silently forcing a full n x n latent
#' GP (n = ~538 plot-level rows here) instead of the intended cheap
#' Hilbert-space approximation, and is what was producing chronic 100%
#' max-treedepth in `brms_cv_compare()`'s fold refits regardless of
#' `adapt_delta`. `20` is a moderate, computationally cheap starting basis
#' size for a 2D spatial term at this domain size (~40km x ~110km, `c = 5/4`);
#' increase it only if posterior spatial predictions look under-resolved.
#' @return a formula.
brms_density_formula <- function(train, gp_coords = c('Longitude', 'Latitude'),
                                  gp_k = 20, gp_c = 5/4){

  covars <- setdiff(names(train), c('Prsnc_All', gp_coords))
  gp_term <- sprintf("gp(%s, %s, k = %s, c = %s)",
                      gp_coords[1], gp_coords[2],
                      if(is.na(gp_k)) 'NA' else gp_k, gp_c)
  as.formula(paste('Prsnc_All ~', paste(covars, collapse = ' + '), '+', gp_term))
}

#' Fit one brms count-model candidate and predict onto a held-out test set.
#'
#' @description Mirrors the `list(Model = , Predictions = )` return shape of
#' `poiss()`/`tweed()`/`gbs()` (../functions.R) so a brms candidate slots into
#' `densityModeller()`'s `mods`/`namev`/`mets()` pipeline without changes to
#' that pipeline.
#'
#' `control = list(adapt_delta = 0.99)` (brms default is 0.8) on every `brm()`
#' call in this file: without it, `brms_cv_compare()`'s spatial `kfold()` arm
#' was observed hitting max-treedepth on the large majority of its 10 fold
#' refits, and then crashing entirely ("Gaussian process covariance matrix is
#' not positive definite") when `.predictor_gp_new()` tried to predict the
#' `gp()` term at held-out fold locations from a poorly-explored posterior.
#' Kept as a general safety margin even after the actual root cause of that
#' pathology was found and fixed (`gp_k` was `NA`, silently requesting an
#' *exact* GP rather than the intended Hilbert-space approximation - see
#' `brms_density_formula()`'s `gp_k` docs). Same `adapt_delta` fix
#' `sagesteppe/safeHavens::bayesianSDM()` already applies to its own GP-term
#' brms model, for the same reason.
#'
#' @param train,test data.frames (not sf) with `Prsnc_All`, `Pr.SuitHab`,
#' `Longitude`, `Latitude` and covariate columns - the same shape
#' `densityModeller()` passes to `poiss()`/`tweed()`/`gbs()`, plus the two
#' coordinate columns for the `gp()` term.
#' @param family a brms family object: `poisson()`, `negbinomial()`,
#' `brms::hurdle_poisson()`, or `brms::hurdle_negbinomial()`.
#' @param gp_coords,gp_k,gp_c as `brms_density_formula()`.
#' @param backend,chains,iter,warmup,cores,seed passed to `brms::brm()`.
#' @return `list(Model = <brmsfit>, Predictions = data.frame(Observed, Predicted, Pr.suit))`.
brms_count_model <- function(train, test, family,
                              gp_coords = c('Longitude', 'Latitude'),
                              gp_k = 20, gp_c = 5/4,
                              backend = 'cmdstanr', chains = 4, iter = 2000,
                              warmup = 1000, cores = chains, seed = 1){

  std <- standardize_covariates(train, test, gp_coords = gp_coords)
  form <- brms_density_formula(std$train, gp_coords, gp_k, gp_c)

  fit <- brms::brm(
    form, data = std$train, family = family, backend = backend,
    prior = brms_default_priors(), init = 0.1, control = list(adapt_delta = 0.99),
    chains = chains, iter = iter, warmup = warmup, cores = cores, seed = seed,
    silent = 2, refresh = 0
  )

  preds <- data.frame(
    Observed  = test$Prsnc_All,
    Predicted = colMeans(brms::posterior_epred(fit, newdata = std$newdata)),
    Pr.suit   = test$Pr.SuitHab
  )

  list(Model = fit, Predictions = preds, center = std$center, scale = std$scale)
}

#' Weakly-informative priors for the brms density candidates.
#'
#' @description brms's own default is a flat (improper) prior on `class = 'b'`
#' (the environmental-covariate slopes) - the actual gap that mattered for the
#' initialization failure (../Census_uncertainty/diagnose_brms_failure.R),
#' since a flat prior does nothing to keep slope draws away from extreme
#' values. Combined with standardized covariates (`standardize_covariates()`),
#' `normal(0, 1)` here is weakly-informative on the (standardized) predictor
#' scale - same family of choice as `sagesteppe/safeHavens::bayesianSDM()`'s
#' `prior_type = "normal"` option (`build_priors()`), just without that
#' function's horseshoe/user-configurable machinery, since this project's
#' covariates are already Boruta/RFE-selected rather than a large
#' possibly-irrelevant set that would call for shrinkage.
#'
#' `normal(0, 1)`, not the originally-tried `normal(0, 5)`: confirmed by
#' testing (`Census_uncertainty/isolate_prior_geometry.R`-style A/B, both
#' with and without the `gp()` term) that `normal(0, 5)` was the actual cause
#' of Phase 1's chronic 100% max-treedepth hurdle-negbinomial fits, not
#' `gp()` itself (removing `gp()` alone left the pathology at 97.5%) and not
#' any single missing covariate (an aridity-PCA composite and a low/high
#' elevation flag were each tried as additions first - both left the
#' pathology untouched or worse). At ~500 plot-level rows across 15 uneven
#' populations and 15 covariates, several slopes (NDVI/SAVI/MSP/MAP/
#' mean_curv/maximal_curv) are only weakly identified by the data - under
#' `normal(0, 5)` their posteriors sprawled out toward the width of the prior
#' itself (`Est.Error` up to ~3.6, `Q2.5`/`Q97.5` spanning roughly -7 to +7),
#' creating exactly the diffuse, hard-to-traverse geometry that produces long
#' NUTS trajectories. `normal(0, 1)` regularizes those weak slopes toward
#' zero instead of leaving them to sprawl, and treedepth saturation dropped
#' to 0% (with `gp()`) with only 3/1000 residual divergences, from 100%
#' max-treedepth under `normal(0, 5)`.
#'
#' Deliberately does **not** touch `class = 'Intercept'` (or `shape`/`hu`):
#' confirmed via `brms::get_prior()` on this exact formula/family that brms's
#' own default Intercept prior is already response-scale-informed
#' (`student_t(3, <data-derived location>, 2.5)`, i.e. it already knows this
#' project's counts run roughly 0-200) - overriding it with a generic normal
#' prior would replace a better default with a worse one. `shape`
#' (`inv_gamma(0.4, 0.3)`) and `hu` (`beta(1, 1)`) defaults are likewise
#' already reasonable, not flat.
#' @return a `brmsprior` object.
brms_default_priors <- function(){
  brms::set_prior('normal(0, 1)', class = 'b')
}

#' Fixed-k spatial fold ids for brms::kfold(), for the density training set.
#'
#' @description `splitData()`'s `nndm_indices` (../functions.R:868-954, built
#' via plain `CAST::nndm()`) is near-leave-one-out - confirmed by testing on
#' synthetic data, its fold count equals `nrow(train)`. That's fine for the
#' XGBoost/LightGBM candidates (`caret`/`tune_race_anova` fits are cheap
#' enough to repeat that many times) but is not feasible for `brms::kfold()`,
#' which does a full Stan refit per fold. This function instead runs
#' `CAST::knndm(train.sf, modeldomain, k = k)` - the same fixed-k spatial CV
#' CAST call the RF suitability model already uses (../functions.R:103,
#' `CAST::knndm(Train.sf, rast_dat, k = 10)`), applied to the density
#' training points instead - and returns one fold id per row (no `NA`s,
#' since `knndm()`'s folds partition every training point, unlike `nndm()`'s
#' buffered near-LOO folds).
#'
#' @param train data.frame with `Longitude`/`Latitude` columns (as passed to
#' `brms_count_model()`).
#' @param k number of spatial folds.
#' @param coords character(2), coordinate column names.
#' @param crs_utm CRS `train`'s coordinates are in (`splitData()` uses UTM 13N, EPSG:32613).
#' @param buffer_m metres to buffer the training points' union by when
#' building the modeldomain surrogate (`splitData()` uses 5000m for the same purpose).
#' @return list(ids = integer vector length `nrow(train)`, knndm = the raw `CAST::knndm()` object).
brms_spatial_fold_ids <- function(train, k = 10, coords = c('Longitude', 'Latitude'),
                                   crs_utm = 32613, buffer_m = 5000){

  train_sf <- sf::st_as_sf(train, coords = coords, crs = crs_utm)
  modeldomain <- sf::st_union(train_sf) |> sf::st_buffer(buffer_m)

  kn <- CAST::knndm(train_sf, modeldomain, k = k)

  ids <- rep(NA_integer_, nrow(train))
  for(i in seq_along(kn$indx_test)){
    ids[kn$indx_test[[i]]] <- i
  }

  list(ids = ids, knndm = kn)
}

#' Cheap-search comparison of the 4 brms families under spatial + non-spatial CV.
#'
#' @description Fits each family once per CV mode at reduced chains/iter (the
#' same cheap-search-then-expensive-final-fit split `adaptive_PAratio_search()`
#' already uses, ../functions.R:755-835), then scores every fit with the same
#' `Observed`/`Predicted` -> `mets()` convention (../functions.R:1050-1059) the
#' rest of `densityModeller()`'s comparison table uses, so results append
#' directly onto the existing `metrrs` table.
#'
#' Non-spatial CV here is single-fit PSIS-LOO (`brms::loo(fit, cores = )`),
#' not a `K`-fold refit loop - matches `sagesteppe/safeHavens::bayesianSDM()`'s
#' evaluation strategy (confirmed by reading its source: no `brms::kfold()`
#' call anywhere in that package, evaluation is `brms::loo()`/`loo::loo()`
#' throughout). One Stan fit plus parallelized (`cores`) importance sampling
#' replaces what would otherwise be `k_nonspatial` serial full refits per
#' family - that serial loop was confirmed to be the slower half of this
#' comparison (the spatial arm still needs actual held-out spatial folds,
#' which LOO doesn't respect, so it keeps `brms::kfold(fit, folds = )`).
#' `moment_match`/`reloo` are deliberately left off here (unlike the full
#' final-fit LOO one might run post-promotion): at `cheap_iter`/`cheap_chains`
#' precision this is a coarse ranking pass, not a publishable estimate, and
#' `reloo`'s exact per-observation refits would reintroduce the same serial
#' cost this change exists to remove.
#'
#' @param train,test as `brms_count_model()`.
#' @param nndm_indices `splitData()$nndm_indices`, for the spatial fold ids.
#' @param families named list of brms family objects to compare.
#' @param cheap_chains,cheap_iter,cheap_warmup reduced-precision fit settings.
#' @param k_spatial number of folds for the spatial (`brms_spatial_fold_ids()`)
#' comparator; the non-spatial comparator is single-fit PSIS-LOO and has no
#' fold count.
#' @param backend,cores,seed passed through to fitting/kfold/loo.
#' @return list(table = data.frame(Model, Metric, Value) matching `metrrs`'s shape,
#' fits = named list of every family/CV-mode fit, best = list(family_name, cv_mode)).
brms_cv_compare <- function(train, test,
                             families = list(
                               Poisson          = poisson(),
                               NegBinomial      = brms::negbinomial(),
                               HurdlePoisson    = brms::hurdle_poisson(),
                               HurdleNegBinomial = brms::hurdle_negbinomial()
                             ),
                             cheap_chains = 1, cheap_iter = 500, cheap_warmup = 250,
                             k_spatial = 10, backend = 'cmdstanr',
                             cores = cheap_chains, seed = 1){

  spat_ids <- brms_spatial_fold_ids(train, k = k_spatial)$ids

  score_one <- function(fam, cv_mode){
    if(cv_mode == 'spatial'){
      keep <- !is.na(spat_ids)
      dat <- train[keep, , drop = FALSE]
    } else {
      dat <- train
    }

    # standardize once per fold-subset (not via brms_count_model(), which would
    # standardize a second time on top of this and double-divide the already-
    # rescaled gp() coordinates) - self-consistent between this fit and the
    # manual posterior_epred() call further down, which is all this needs.
    std <- standardize_covariates(dat)

    fit <- brms::brm(
      brms_density_formula(std$train), data = std$train, family = families[[fam]],
      backend = backend, prior = brms_default_priors(), init = 0.1,
      control = list(adapt_delta = 0.99),
      chains = cheap_chains, iter = cheap_iter, warmup = cheap_warmup,
      cores = cores, seed = seed, silent = 2, refresh = 0
    )

    if(cv_mode == 'spatial'){
      cv <- brms::kfold(fit, folds = spat_ids[keep], chains = cheap_chains,
                         iter = cheap_iter, warmup = cheap_warmup, cores = cores)
      elpd <- cv$estimates['elpd_kfold', 'Estimate']
    } else {
      cv <- brms::loo(fit, cores = cores)
      elpd <- cv$estimates['elpd_loo', 'Estimate']
    }

    # per-fold held-out prediction, scored the same way as the rest of the table.
    preds <- data.frame(Observed = dat$Prsnc_All,
                         Predicted = colMeans(brms::posterior_epred(fit, newdata = std$train)))

    list(fit = fit, cv = cv, elpd = elpd, metrics = mets(preds))
  }

  combos <- expand.grid(family = names(families), cv_mode = c('spatial', 'non_spatial'),
                         stringsAsFactors = FALSE)

  results <- Map(score_one, combos$family, combos$cv_mode)
  names(results) <- paste(combos$family, combos$cv_mode, sep = '_')

  namev <- paste0(
    ifelse(combos$family == 'NegBinomial', 'NegBin',
           ifelse(combos$family == 'HurdlePoisson', 'Hurdle Poisson',
                  ifelse(combos$family == 'HurdleNegBinomial', 'Hurdle NegBin', combos$family))),
    ifelse(combos$cv_mode == 'spatial', ' Spat.', '')
  )

  table <- dplyr::bind_rows(lapply(results, `[[`, 'metrics')) |>
    dplyr::mutate(Model = rep(namev, each = 3), .before = 1)

  elpd <- sapply(results, `[[`, 'elpd')
  best_combo <- combos[which.max(elpd), ]

  list(table = table, fits = results,
       best = list(family_name = best_combo$family, cv_mode = best_combo$cv_mode))
}

#' Refit the winning brms family at full precision on the full training data.
#'
#' @description The expensive-final-fit half of Phase 1: no folds, full
#' chains/iterations, `backend = 'cmdstanr'`, caches to disk with the same
#' `if(file.exists(f)) readRDS() else fit + saveRDS()` idiom `modeller()`/
#' `densityModeller()` already use (../functions.R). Also extracts and saves
#' the posterior draws matrix Phase 3 consumes.
#'
#' @param best_family a brms family object (the winner from `brms_cv_compare()$best`).
#' @param train,seed as `brms_count_model()`; no `test` arg since this refit
#' uses every training row.
#' @param fp,bn as `densityModeller()` - used to build the cache path.
#' @param family_label short label used in the cached filenames (e.g. "hurdle_negbinomial").
#' @param full_chains,full_iter,full_warmup,cores full-precision fit settings.
#' @return list(Model = <brmsfit>, draws = <draws_matrix>, model_path =,
#' draws_path =, scaling_path = <path to the saved center/scale used to fit
#' this model - Phase 3/`combine_census.R` and Phase 4/`validate_pipeline.R`
#' need to `apply_standardization()` any new prediction data with this exact
#' transform before calling `posterior_epred()` on this model>).
brms_promote_and_refit <- function(best_family, train, seed, fp, bn, family_label,
                                    full_chains = 4, full_iter = 2000, full_warmup = 1000,
                                    cores = full_chains, backend = 'cmdstanr'){

  model_path   <- file.path(fp, 'models', paste0(bn, '-brms-', family_label, '.rds'))
  draws_path   <- file.path(fp, 'models', paste0(bn, '-brms-', family_label, '-draws.rds'))
  scaling_path <- file.path(fp, 'models', paste0(bn, '-brms-', family_label, '-scaling.rds'))

  if(!file.exists(model_path)){
    std <- standardize_covariates(train)
    saveRDS(list(center = std$center, scale = std$scale), scaling_path)

    form <- brms_density_formula(std$train)
    fitted <- brms::brm(
      form, data = std$train, family = best_family, backend = backend,
      prior = brms_default_priors(), init = 0.1, control = list(adapt_delta = 0.99),
      chains = full_chains, iter = full_iter, warmup = full_warmup,
      cores = cores, seed = seed, silent = 2, refresh = 0
    )
    saveRDS(fitted, model_path)
  } else {
    fitted <- readRDS(model_path)
  }

  if(!file.exists(draws_path)){
    draws <- posterior::as_draws_matrix(fitted)
    saveRDS(draws, draws_path)
  } else {
    draws <- readRDS(draws_path)
  }

  list(Model = fitted, draws = draws, model_path = model_path,
       draws_path = draws_path, scaling_path = scaling_path)
}
