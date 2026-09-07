## Phase 3 of census_uncertainty_roadmap.md - Combine.
##
## The boundary ensemble (Phase 2, boundary_simulation.R) and the density
## posterior (Phase 1, density_bayes.R) are fit independently with no shared
## MCMC chain, so their draws are combined by Monte Carlo composition rather
## than a single joint model - see the roadmap's Phase 3 notes for why.

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
#' @param brms_fit the promoted `brms_promote_and_refit()$Model`.
#' @param raster_template a `terra::rast` covering the full prediction domain,
#' one layer per covariate the model's formula needs, named to match.
#' @param draw_id integer, which posterior draw to predict.
#' @param covariate_df optional, pre-extracted `as.data.frame(raster_template, cells=TRUE)`
#' - pass this in when calling repeatedly (e.g. per draw in
#' `combined_census_montecarlo()`) so it's built once, not on every call.
#' @param chunk_cells max rows of `covariate_df` per `posterior_epred()` call;
#' `Inf` (default) means one call for the whole raster - see caveat above.
#' @param seed RNG seed for this draw's `gp()` resampling at new locations;
#' defaults to `draw_id` so the same draw is always reproducible.
#' @return a `terra::rast` (single layer) of predicted density for this draw.
predict_density_draw <- function(brms_fit, raster_template, draw_id,
                                  covariate_df = NULL, chunk_cells = Inf, seed = draw_id){

  if(is.null(covariate_df)){
    covariate_df <- as.data.frame(raster_template, cells = TRUE)
  }

  out <- terra::rast(raster_template, nlyrs = 1)

  chunk_cells <- min(chunk_cells, nrow(covariate_df))  # seq(..., by = Inf) gives NaN, not one chunk
  row_starts <- seq(1, nrow(covariate_df), by = chunk_cells)
  for(s in row_starts){
    e <- min(s + chunk_cells - 1, nrow(covariate_df))
    chunk <- covariate_df[s:e, , drop = FALSE]

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
#' @param brms_fit,raster_template as `predict_density_draw()`.
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
combined_census_montecarlo <- function(brms_fit, raster_template, boundary_mask_paths,
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
    dens_r <- predict_density_draw(brms_fit, raster_template, pairs$density_id[i], covariate_df = covariate_df)
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
