## Phase 1 of census_uncertainty_roadmap.md - Bayesian density candidates.
##
## Fits brms count-model families (poisson, negbinomial, hurdle_poisson,
## hurdle_negbinomial) as additional candidates in the same comparison table
## `densityModeller()` (../functions.R) already builds for the five ML
## candidates (XGBoost Poisson/Tweedie spatial+non-spatial, LightGBM Poisson
## spatial). Source ../functions.R before this file - `brms_cv_compare()` and
## `brms_promote_and_refit()` are called from `densityModeller()`.

#' Build the shared brms formula: response ~ covariates + approximate gp(x, y).
#'
#' @description Uses brms's approximate (Hilbert-space) `gp()` term
#' (`gp_k`/`gp_c`) rather than an exact GP, so the fitted posterior stays a
#' plain coefficient+basis-weight matrix (see `brms_promote_and_refit()`) and
#' stays fittable in Stan at plot-level n.
#' @param train data.frame the formula's covariates are read off of (every
#' column except `Prsnc_All` and `gp_coords`).
#' @param gp_coords character(2), coordinate column names for the `gp()` term.
#' @param gp_k,gp_c approximate-GP basis dimension / boundary factor. `gp_k = NA`
#' lets brms pick a default basis size.
#' @return a formula.
brms_density_formula <- function(train, gp_coords = c('Longitude', 'Latitude'),
                                  gp_k = NA, gp_c = 5/4){

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
                              gp_k = NA, gp_c = 5/4,
                              backend = 'cmdstanr', chains = 4, iter = 2000,
                              warmup = 1000, cores = chains, seed = 1){

  form <- brms_density_formula(train, gp_coords, gp_k, gp_c)

  fit <- brms::brm(
    form, data = train, family = family, backend = backend,
    chains = chains, iter = iter, warmup = warmup, cores = cores, seed = seed,
    silent = 2, refresh = 0
  )

  preds <- data.frame(
    Observed  = test$Prsnc_All,
    Predicted = colMeans(brms::posterior_epred(fit, newdata = test)),
    Pr.suit   = test$Pr.SuitHab
  )

  list(Model = fit, Predictions = preds)
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
#' Non-spatial CV here is an ordinary random `K`-fold partition
#' (`brms::kfold(fit, K = )`), not `rsample::bootstraps(train, 10)`
#' (../functions.R:1123) - `brms::kfold()`'s `folds` argument needs a disjoint
#' partition, which a bootstrap resample (with replacement) isn't. Still
#' "non-spatial CV" in spirit, just not the identical resampling object the
#' XGBoost/LightGBM candidates use.
#'
#' @param train,test as `brms_count_model()`.
#' @param nndm_indices `splitData()$nndm_indices`, for the spatial fold ids.
#' @param families named list of brms family objects to compare.
#' @param cheap_chains,cheap_iter,cheap_warmup reduced-precision fit settings.
#' @param k_spatial,k_nonspatial number of folds for the spatial
#' (`brms_spatial_fold_ids()`) and non-spatial `brms::kfold()` comparators.
#' @param backend,cores,seed passed through to fitting/kfold.
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
                             k_spatial = 10, k_nonspatial = 10, backend = 'cmdstanr',
                             cores = cheap_chains, seed = 1){

  spat_ids <- brms_spatial_fold_ids(train, k = k_spatial)$ids

  score_one <- function(fam, cv_mode){
    if(cv_mode == 'spatial'){
      keep <- !is.na(spat_ids)
      dat <- train[keep, , drop = FALSE]
    } else {
      dat <- train
    }

    fit <- brms_count_model(dat, test, family = families[[fam]],
                             backend = backend, chains = cheap_chains,
                             iter = cheap_iter, warmup = cheap_warmup,
                             cores = cores, seed = seed)$Model

    kf <- if(cv_mode == 'spatial'){
      brms::kfold(fit, folds = spat_ids[keep], chains = cheap_chains,
                  iter = cheap_iter, warmup = cheap_warmup, cores = cores)
    } else {
      brms::kfold(fit, K = k_nonspatial, chains = cheap_chains,
                  iter = cheap_iter, warmup = cheap_warmup, cores = cores)
    }

    # per-fold held-out prediction, scored the same way as the rest of the table.
    preds <- data.frame(Observed = dat$Prsnc_All,
                         Predicted = colMeans(brms::posterior_epred(fit, newdata = dat)))

    list(fit = fit, kfold = kf, metrics = mets(preds))
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

  elpd <- sapply(results, function(r) r$kfold$estimates['elpd_kfold', 'Estimate'])
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
#' @return list(Model = <brmsfit>, draws = <draws_matrix>, model_path =, draws_path =).
brms_promote_and_refit <- function(best_family, train, seed, fp, bn, family_label,
                                    full_chains = 4, full_iter = 2000, full_warmup = 1000,
                                    cores = full_chains, backend = 'cmdstanr'){

  model_path <- file.path(fp, 'models', paste0(bn, '-brms-', family_label, '.rds'))
  draws_path <- file.path(fp, 'models', paste0(bn, '-brms-', family_label, '-draws.rds'))

  if(!file.exists(model_path)){
    form <- brms_density_formula(train)
    fitted <- brms::brm(
      form, data = train, family = best_family, backend = backend,
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

  list(Model = fitted, draws = draws, model_path = model_path, draws_path = draws_path)
}
