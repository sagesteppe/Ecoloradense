## Phase 3 of census_uncertainty_roadmap.md - Combine.
##
## The boundary ensemble (Phase 2, boundary_simulation.R) and the density
## posterior (Phase 1, density_bayes.R) are fit independently with no shared
## MCMC chain, so their draws are combined by Monte Carlo composition rather
## than a single joint model - see the roadmap's Phase 3 notes for why.

#' Assemble the covariate raster template both Phase 3 draw sources need.
#'
#' @description Neither `predict_density_draw()` (brms) nor
#' `predict_density_draw_conformal()` (ML candidates, density_conformal.R)
#' can be handed just the suitability probability raster alone - both need
#' every covariate column the promoted model's formula/feature set actually
#' uses, named to match (confirmed by testing: `stats::predict()` on a
#' `parsnip` fit errors loudly, "object '<covariate>' not found", the moment
#' a required column is missing - a safe failure mode, but a real gap
#' `run_census_uncertainty()` had left unaddressed, passing `terra::rast(pr_path)`
#' - a single layer - directly as `raster_template`). This assembles the
#' real thing: `rastReader()`'s full covariate stack (../functions.R) plus
#' the suitability raster itself, attached as `Pr.SuitHab` (resampled onto
#' the covariate stack's grid first if it doesn't already match exactly -
#' confirmed necessary in practice, not just defensive).
#'
#' @param p2proc path to processed raster data (`rastReader()`'s `p2proc` arg).
#' @param pr_path path to the promoted product's suitability probability
#' raster (`results/suitability_maps/<bn>DO:0-Seed....-Pr.tif`).
#' @param region_ext optional `terra::ext`/SpatExtent to crop to before
#' returning - strongly recommended for anything but a final full-domain
#' run: the full domain is ~330M cells (confirmed against the existing
#' suitability raster), a many-hour, tens-of-GB job to predict across
#' unclipped - see census_uncertainty_roadmap.md Phase 3's "mask tightly to
#' regional bounding boxes... first". `combined_census_montecarlo()`/
#' `combined_census_montecarlo_conformal()` crop to `region_bbox` internally
#' too, so passing the same extent here is redundant-but-harmless with those
#' - it mainly matters for a raster handed to Phase 4's
#' `check_predict_density_draw()`, which does not crop internally.
#' @return a `terra::rast` with `rastReader()`'s full covariate stack plus a
#' `Pr.SuitHab` layer, ready to pass as `raster_template` to either
#' `predict_density_draw()` or `predict_density_draw_conformal()`.
build_density_raster_template <- function(p2proc, pr_path, region_ext = NULL){

  rast_dat <- rastReader('dem_1-3arc', p2proc)
  suit_r <- terra::rast(pr_path)

  if(!is.null(region_ext)){
    rast_dat <- terra::crop(rast_dat, region_ext)
    suit_r <- terra::crop(suit_r, region_ext)
  }

  if(!terra::compareGeom(suit_r, rast_dat, stopOnError = FALSE)){
    suit_r <- terra::resample(suit_r, rast_dat, method = 'bilinear')
  }
  names(suit_r) <- 'Pr.SuitHab'
  rast_dat$Pr.SuitHab <- suit_r

  rast_dat
}

#' Pair boundary draws against density draws for the combination loop.
#'
#' @description Resamples the smaller boundary stack with replacement,
#' consuming as many *distinct* density draws as available before repeating
#' any - per the roadmap's "don't throw away density draws you already paid
#' for."
#'
#' @param n_boundary,n_density number of available boundary/density draws.
#' @param n_combined number of combined pairs to produce.
#' @param seed RNG seed.
#' @return data.frame(boundary_id, density_id), nrow = n_combined.
pair_draws <- function(n_boundary, n_density, n_combined, seed = 1){

  set.seed(seed)
  boundary_id <- sample.int(n_boundary, n_combined, replace = n_combined > n_boundary)
  density_id  <- sample.int(n_density,  n_combined, replace = n_combined > n_density)
  data.frame(boundary_id = boundary_id, density_id = density_id)
}

#' Predict one posterior draw's density surface onto a raster template.
#'
#' @description Builds `newdata` from `raster_template`'s cell values and
#' calls `brms::posterior_epred(..., draw_ids = draw_id)` for a single draw -
#' the roadmap's "predict -> use -> discard per draw" (Phase 3), using brms's
#' own posterior_epred rather than a hand-rolled GP-basis matmul (a design
#' fork discussed and settled earlier: brms's built-in prediction, not a
#' custom GPU path).
#'
#' **Important, confirmed by testing**: for a model with a `gp()` term,
#' `posterior_epred()` at genuinely new (out-of-sample) locations is *not*
#' deterministic given only `draw_ids` - it draws the GP's latent value at
#' each new location from R's global RNG on every call (fixing `set.seed()`
#' first *does* make it reproducible). Two consequences: (1) this function
#' seeds on `draw_id` itself so the same draw always reproduces the same
#' surface - required for Phase 3's Monte Carlo composition, since
#' `pair_draws()` resamples density draws with replacement and a reused
#' `draw_id` must give back the identical realization, not a fresh random
#' one; (2) **do not chunk `newdata` across multiple `posterior_epred()`
#' calls for a `gp()` model** - each call independently resamples the GP at
#' its own subset of new locations, so adjacent chunks would not be
#' consistent pieces of one spatially-coherent realization, only
#' independently-conditioned fragments. `chunk_cells` therefore defaults to
#' "no chunking"; only lower it if `newdata` truly won't fit in memory for
#' one call, and treat that raster as an approximation, not an exact draw.
#'
#' @param brms_fit the promoted `brms_promote_and_refit()$Model` - fit on
#' standardized covariates (`density_bayes.R`'s `standardize_covariates()`), so
#' `center`/`scale` (its `-scaling.rds`) must be supplied here to transform
#' `raster_template`'s raw covariate values into the same scale before
#' `posterior_epred()` - otherwise predictions are silently wrong (the model's
#' coefficients are in standardized-covariate units, not the raster's raw units).
#' @param raster_template a `terra::rast` covering the full prediction domain,
#' one layer per covariate the model's formula needs, named to match.
#' @param draw_id integer, which posterior draw to predict.
#' @param center,scale named numeric vectors from the promoted model's
#' `-scaling.rds` (`brms_promote_and_refit()$scaling_path`).
#' @param covariate_df optional, pre-extracted `as.data.frame(raster_template, cells=TRUE)`
#' - pass this in when calling repeatedly (e.g. per draw in
#' `combined_census_montecarlo()`) so it's built once, not on every call. Must
#' still be in raw (unstandardized) units - the `center`/`scale` transform is
#' applied inside this function, once per chunk.
#' @param chunk_cells max rows of `covariate_df` per `posterior_epred()` call;
#' `Inf` (default) means one call for the whole raster - see caveat above.
#' @param seed RNG seed for this draw's `gp()` resampling at new locations;
#' defaults to `draw_id` so the same draw is always reproducible.
#' @return a `terra::rast` (single layer) of predicted density for this draw.
predict_density_draw <- function(brms_fit, raster_template, draw_id, center, scale,
                                  covariate_df = NULL, chunk_cells = Inf, seed = draw_id){

  if(is.null(covariate_df)){
    covariate_df <- as.data.frame(raster_template, cells = TRUE)
  }

  out <- terra::rast(raster_template, nlyrs = 1)

  chunk_cells <- min(chunk_cells, nrow(covariate_df))  # seq(..., by = Inf) gives NaN, not one chunk
  row_starts <- seq(1, nrow(covariate_df), by = chunk_cells)
  for(s in row_starts){
    e <- min(s + chunk_cells - 1, nrow(covariate_df))
    chunk <- apply_standardization(covariate_df[s:e, , drop = FALSE], center, scale)

    set.seed(seed)
    pred <- brms::posterior_epred(brms_fit, newdata = chunk, draw_ids = draw_id)
    out[chunk$cell] <- as.vector(pred)
  }
  out
}

#' Align one Phase 2 boundary mask onto a target raster's exact geometry.
#'
#' @description Phase 2's masks are built on the RF suitability surface's
#' full extent; the density side works on `raster_template` cropped to
#' `region_bbox`. `terra::mask()` requires matching geometry (same extent/
#' resolution/CRS) between its two arguments, so every mask has to be
#' aligned to the (already-cropped) density raster before use - resampling
#' rather than just cropping guards against any grid-origin mismatch between
#' the suitability and covariate rasters, and `method = 'near'` is used since
#' masks are binary/categorical, not continuous.
#'
#' @param mask_path path to one Phase 2 mask file.
#' @param target a `terra::rast` whose geometry the mask should match.
#' @return the aligned mask, a `terra::rast`.
align_mask_to_template <- function(mask_path, target){
  terra::resample(terra::rast(mask_path), target, method = 'near')
}

#' Combine one boundary/density draw pair into a single census-size value.
#'
#' @param density_raster one `predict_density_draw()` output.
#' @param aligned_mask one `align_mask_to_template()` output - already sharing
#' `density_raster`'s exact geometry, so `terra::mask()` is safe here.
#' @param cell_area numeric, area of one raster cell (e.g. `prod(terra::res(density_raster))`).
#' @return numeric scalar, the census-size estimate for this draw pair.
census_from_pair <- function(density_raster, aligned_mask, cell_area){

  masked <- terra::mask(density_raster, aligned_mask)
  total <- terra::global(masked, fun = 'sum', na.rm = TRUE)[1, 1]
  # terra::global(..., na.rm=TRUE) returns NA (not 0, unlike base R's sum())
  # when every cell is NA - a draw whose mask excludes everything is a
  # genuine, meaningful "0 individuals" outcome, not a missing value, and
  # must not silently drop out of the Monte Carlo draws vector.
  if(is.na(total)) total <- 0
  as.numeric(total) * cell_area
}

#' Running 2.5/97.5 percentile of a growing vector.
#' @param x numeric vector.
#' @return c(lower, upper).
running_ci <- function(x) stats::quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)

#' Monte Carlo combination of the density posterior and boundary ensemble
#' into a census-size credible interval.
#'
#' @description Masks `raster_template` to `region_bbox` first (roadmap's
#' "mask tightly to regional bounding boxes used in manuscript"), then loops
#' paired draws up to `n_max`, predicting -> masking -> summing -> discarding
#' each draw's raster (never holds more than one density raster in memory at
#' once). Once `length(draws) >= min_n` (respecting the `ess_tail` ~400 floor
#' the roadmap cites from Stan's own tail-ESS diagnostic), recomputes the
#' running 2.5/97.5 percentile each iteration and stops once it stabilizes
#' within `stop_tol` (relative change) over a trailing window of `window`
#' draws, or at `n_max`.
#'
#' @param brms_fit,raster_template,center,scale as `predict_density_draw()`.
#' @param boundary_mask_paths character vector, Phase 2's mask file paths.
#' @param region_bbox a `terra::ext`/SpatExtent (or object `terra::crop()`
#' accepts) to mask `raster_template` to before looping.
#' @param n_max maximum combined draws.
#' @param min_n minimum draws before checking for stabilization (default 400).
#' @param window trailing-window size (in draws) used to judge stabilization.
#' @param stop_tol relative-change tolerance on the CI bounds to declare
#' stabilization.
#' @param out_csv path to write the draws + running-quantile trace to.
#' @return list(draws = numeric vector, ci = c(lower, upper), n = length(draws)).
combined_census_montecarlo <- function(brms_fit, raster_template, center, scale,
                                        boundary_mask_paths,
                                        region_bbox, n_max = 2000, min_n = 400,
                                        window = 100, stop_tol = 0.01, seed = 1,
                                        out_csv = NULL){

  raster_template <- terra::crop(raster_template, region_bbox)
  cell_area <- prod(terra::res(raster_template))

  # built/aligned once, reused across all draws - covariate extraction and mask
  # alignment are both comparatively cheap (small relative to a posterior_epred()
  # call), and masks especially get reused often since they're resampled with
  # replacement against the larger density stack.
  covariate_df <- as.data.frame(raster_template, cells = TRUE)
  aligned_masks <- lapply(boundary_mask_paths, align_mask_to_template, target = raster_template)

  n_boundary <- length(boundary_mask_paths)
  n_density <- nrow(posterior::as_draws_matrix(brms_fit))

  pairs <- pair_draws(n_boundary, n_density, n_max, seed = seed)

  draws <- numeric(n_max)
  trace <- matrix(NA_real_, nrow = n_max, ncol = 2)
  n_done <- 0

  for(i in seq_len(n_max)){
    dens_r <- predict_density_draw(brms_fit, raster_template, pairs$density_id[i], center, scale,
                                    covariate_df = covariate_df)
    draws[i] <- census_from_pair(dens_r, aligned_masks[[pairs$boundary_id[i]]], cell_area)
    n_done <- i

    if(i >= min_n){
      trace[i, ] <- running_ci(draws[1:i])
      if(i >= min_n + window){
        prev <- trace[i - window, ]
        rel_change <- abs(trace[i, ] - prev) / pmax(abs(prev), .Machine$double.eps)
        if(all(rel_change < stop_tol)) break
      }
    }
  }

  draws <- draws[1:n_done]
  trace <- trace[1:n_done, , drop = FALSE]

  if(!is.null(out_csv)){
    write.csv(
      data.frame(draw = seq_len(n_done), census_size = draws,
                 running_lower = trace[, 1], running_upper = trace[, 2]),
      out_csv, row.names = FALSE
    )
  }

  list(draws = draws, ci = running_ci(draws), n = n_done)
}

#' Monte Carlo combination of an ML candidate's conformal pseudo-draws and the
#' boundary ensemble into a census-size credible interval.
#'
#' @description The conformal-pseudo-draw sibling of `combined_census_montecarlo()`
#' (census_uncertainty_roadmap.md Phase 1c/3), for whichever ML candidate is in
#' play instead of the brms posterior. Structurally identical - mask to
#' `region_bbox`, build `covariate_df`/aligned masks once, loop draws via
#' `predict_density_draw_conformal()` (`density_conformal.R`) instead of
#' `predict_density_draw()`, `census_from_pair()`, running-quantile
#' stabilization via `running_ci()` - but deliberately does **not** share code
#' with `combined_census_montecarlo()`: that path is already Phase 4-validated
#' for reproducibility, and duplicating this ~15-line loop is lower risk than
#' refactoring it into a shared helper both paths would depend on.
#'
#' Unlike `pair_draws()`, there's no finite `n_density` to sample an id from -
#' conformal pseudo-draws are generated fresh by resampling `conformal`'s
#' residual/fold pool on every call to `predict_density_draw_conformal()`, not
#' drawn from a fixed pre-existing posterior - so only the boundary side needs
#' an explicit resampled index here.
#'
#' @param conformal one ML candidate's `conformal_calibrate_candidate()`
#' output (`densityModeller()`'s `$conformal[[<candidate name>]]`).
#' @param method 'split' or 'cv_plus' - which conformal calibration to draw from.
#' @param raster_template,boundary_mask_paths,region_bbox,n_max,min_n,window,stop_tol,seed,out_csv
#' as `combined_census_montecarlo()`.
#' @return list(draws = numeric vector, ci = c(lower, upper), n = length(draws)).
combined_census_montecarlo_conformal <- function(conformal, method, raster_template,
                                                  boundary_mask_paths, region_bbox,
                                                  n_max = 2000, min_n = 400, window = 100,
                                                  stop_tol = 0.01, seed = 1, out_csv = NULL){

  raster_template <- terra::crop(raster_template, region_bbox)
  cell_area <- prod(terra::res(raster_template))

  covariate_df <- as.data.frame(raster_template, cells = TRUE)
  aligned_masks <- lapply(boundary_mask_paths, align_mask_to_template, target = raster_template)

  n_boundary <- length(boundary_mask_paths)
  set.seed(seed)
  boundary_ids <- sample.int(n_boundary, n_max, replace = n_max > n_boundary)

  draws <- numeric(n_max)
  trace <- matrix(NA_real_, nrow = n_max, ncol = 2)
  n_done <- 0

  for(i in seq_len(n_max)){
    dens_r <- predict_density_draw_conformal(conformal, method, raster_template, draw_id = i,
                                              covariate_df = covariate_df)
    draws[i] <- census_from_pair(dens_r, aligned_masks[[boundary_ids[i]]], cell_area)
    n_done <- i

    if(i >= min_n){
      trace[i, ] <- running_ci(draws[1:i])
      if(i >= min_n + window){
        prev <- trace[i - window, ]
        rel_change <- abs(trace[i, ] - prev) / pmax(abs(prev), .Machine$double.eps)
        if(all(rel_change < stop_tol)) break
      }
    }
  }

  draws <- draws[1:n_done]
  trace <- trace[1:n_done, , drop = FALSE]

  if(!is.null(out_csv)){
    write.csv(
      data.frame(draw = seq_len(n_done), census_size = draws,
                 running_lower = trace[, 1], running_upper = trace[, 2]),
      out_csv, row.names = FALSE
    )
  }

  list(draws = draws, ci = running_ci(draws), n = n_done)
}
