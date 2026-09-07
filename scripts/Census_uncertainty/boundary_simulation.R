## Phase 2 of census_uncertainty_roadmap.md - Boundary simulation.
##
## No Bayesian refit of suitability: reuses the existing RF mean (`-Pr.tif`)
## and infinitesimal-jackknife SE (`-SE.tif`) surfaces `modeller()` already
## produces (../functions.R:216-390) as-is. Estimates the RF residuals'
## spatial covariance, simulates spatially-coherent realizations of the
## suitability surface via FFT/circulant embedding (`fields::circulantEmbedding`),
## and wombles each realization into a boundary/inclusion mask.

#' Estimate the spatial covariance structure of the RF suitability residuals.
#'
#' @description Fits an empirical + theoretical variogram to
#' (observed occurrence - RF mean prediction) at held-out points, the same
#' `gstat::variogram()`/`gstat::fit.variogram()` machinery `densityModeller()`'s
#' kriging null already uses (../functions.R:1090-1094). The result
#' parameterizes the stationary/isotropic covariance `simulate_boundary_fields()`
#' asks `fields::circulantEmbeddingSetup()` to simulate.
#'
#' @param pr_raster the RF mean-prediction raster (`terra::rast`), a `modeller()`
#' `-Pr.tif` output.
#' @param holdout_sf an sf POINT object with an `Occurrence` column (0/1 or
#' continuous suitability) at held-out locations - the same test set `modeller()`
#' evaluates against.
#' @param model_types gstat short model codes to try (`fit.variogram()` picks
#' the best-fitting one via SSE).
#' @return list(vgm_model = the fitted gstat variogramModel, fields_cov =
#' character, one of "Exponential"/"Gaussian"/"Matern" (gstat model codes not
#' in this set fall back to "Exponential", noted via a message - the
#' Monte Carlo boundary ensemble is a simulation envelope, not a
#' likelihood-critical fit, so this approximation is acceptable), theta =
#' numeric (gstat's fitted range, passed directly as fields' `theta`; the two
#' packages' range parameterizations aren't identical, another acceptable
#' approximation for the same reason).
estimate_residual_covariance <- function(pr_raster, holdout_sf,
                                          model_types = c('Exp', 'Gau', 'Sph', 'Mat')){

  holdout_sf$resid <- holdout_sf$Occurrence - terra::extract(pr_raster, holdout_sf)[, 2]

  emp_vgm <- gstat::variogram(resid ~ 1, holdout_sf)
  vgm_model <- gstat::fit.variogram(emp_vgm, gstat::vgm(model_types))

  spatial_rows <- vgm_model$model != 'Nug'
  if(!any(spatial_rows) || max(vgm_model$range[spatial_rows]) <= 0){
    stop("estimate_residual_covariance(): the best-fitting variogram is pure ",
         "nugget (or has a zero range) - these holdout residuals show no ",
         "detectable spatial structure at this sample density, so there is ",
         "nothing for circulant embedding to simulate. Check holdout_sf's ",
         "coverage/density before proceeding.")
  }
  best <- which(spatial_rows)[which.max(vgm_model$psill[spatial_rows])]
  code <- as.character(vgm_model$model[best])
  # fields::circulantEmbeddingSetup()/stationary.cov() resolve `Covariance` by
  # do.call(Covariance, list(d = <distance matrix>, ...)) - it has to be an
  # actual function taking a `d` argument this way. Exponential() and Matern()
  # both fit that convention; fields' Gaussian covariance (gauss.cov(), a
  # thin Exp.cov(..., p=2) wrapper) uses a different, coordinate-based
  # convention (x1/x2, not d) and isn't a drop-in fit here, so it - and any
  # other gstat model code without a direct fit - falls back to Exponential.
  fields_cov <- switch(code, Exp = 'Exponential', Mat = 'Matern',
                        { message(sprintf("gstat model '%s' has no direct fields:: equivalent under this API, using Exponential", code))
                          'Exponential' })

  list(vgm_model = vgm_model, fields_cov = fields_cov, theta = vgm_model$range[best])
}

#' Simulate a spatially-coherent boundary-field ensemble via circulant embedding.
#'
#' @description Draws `n_draws` realizations of the suitability surface as
#' `pr_raster + unit_field * se_raster`, where `unit_field` is a stationary,
#' unit-variance Gaussian random field simulated with
#' `fields::circulantEmbedding()` (O(n log n), FFT-based - the roadmap's
#' explicit alternative to a dense covariance matrix + Cholesky, which would
#' be a ~1TB memory wall at full 1/3-arc resolution). `fields::circulantEmbeddingSetup()`
#' handles the periodicity-avoiding grid padding internally (its `bigGrid`) -
#' no manual domain-doubling needed. Each realization is written to disk
#' immediately (`out_dir/draw_<i>.tif`) rather than held in memory, the same
#' "predict -> write -> discard" idiom `modeller()`'s tile loop already uses
#' (../functions.R:223-264).
#'
#' @param pr_raster,se_raster `terra::rast` mean and SE surfaces, same extent/resolution.
#' @param cov_params `estimate_residual_covariance()`'s return value.
#' @param n_draws number of realizations.
#' @param out_dir directory realizations are written to.
#' @param seed RNG seed.
#' @return character vector of written file paths.
simulate_boundary_fields <- function(pr_raster, se_raster, cov_params, n_draws = 500,
                                      out_dir, seed = 1){

  # fields::stationary.cov() resolves cov.args$Covariance (e.g. "Exponential")
  # by name lookup on the search path, not by namespace - it isn't found via
  # fields:: qualification alone, so the namespace needs to actually be attached.
  if(!'package:fields' %in% search()) library(fields)

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  grid <- list(x = terra::xFromCol(pr_raster, 1:terra::ncol(pr_raster)),
               y = terra::yFromRow(pr_raster, terra::nrow(pr_raster):1))

  # circulantEmbeddingSetup()'s default padding (M ~= 2x each grid dimension)
  # can leave small negative spectral weights from floating-point error rather
  # than a genuinely invalid covariance - confirmed by testing: for a
  # 60x60 grid, the default M=128 gives min(Re(wght)) ~ -9e-6, and bumping to
  # M=486+ (~8x) brings it to a small *positive* number with no other change.
  # Retry with growing padding before giving up.
  fft_sizes <- sort(unique(2^(0:15) %o% 3^(0:15)))
  fft_sizes <- fft_sizes[fft_sizes <= 1e6]
  m_max <- max(length(grid$x), length(grid$y))

  setup <- NULL
  for(pad_mult in c(2, 4, 8, 16, 32)){
    M <- rep(min(fft_sizes[fft_sizes >= pad_mult * m_max]), 2)
    setup <- fields::circulantEmbeddingSetup(
      grid, M = M, cov.args = list(Covariance = cov_params$fields_cov, theta = cov_params$theta),
      giveWarnings = FALSE
    )
    if(!any(Re(setup$wght) < 0)) break
  }
  if(any(Re(setup$wght) < 0)){
    stop("circulantEmbeddingSetup still produced negative spectral weights for ",
         "theta = ", cov_params$theta, " even at ", pad_mult, "x padding - this ",
         "looks like a genuinely invalid covariance/grid combination, not just ",
         "under-padding. Try a coarser check of cov_params (a smaller theta, or ",
         "refit the variogram) before increasing n_draws.")
  }

  set.seed(seed)
  paths <- character(n_draws)
  for(i in seq_len(n_draws)){
    unit_field <- fields::circulantEmbedding(setup)
    unit_field <- unit_field[nrow(unit_field):1, ]  # fields' y runs bottom-to-top; terra's rows run top-to-bottom
    realization <- pr_raster + terra::setValues(se_raster, as.vector(t(unit_field))) * se_raster

    f <- file.path(out_dir, paste0('draw_', i, '.tif'))
    terra::writeRaster(realization, f, overwrite = TRUE)
    paths[i] <- f
  }
  paths
}

#' Womble one realization into a boundary/inclusion mask.
#'
#' @description Not the full model-based boundary-likelihood curve of
#' Womble (1951) / Banerjee & Gelfand's Bayesian wombling (evaluated as an
#' option and set aside - see census_uncertainty_roadmap.md's Phase 2 notes:
#' the reference implementation, github.com/arh926/spWombling, had several
#' bugs, two of which got fixed and PR'd upstream - github.com/arh926/spWombling/pull/3
#' - but a further unresolved numerical issue and an inactive maintainer made
#' it impractical to depend on here). Instead: compute gradient magnitude via
#' the Sobel-Feldman operator (Sobel & Feldman 1969), `spatialEco::sobal(r,
#' method = 'intensity')` - a tested, citable, terra-native implementation,
#' not hand-rolled - then find the realization's own value level where mean
#' gradient magnitude peaks (its steepest-change "ridge") and threshold at
#' that level, a data-driven per-draw analogue of the deterministic
#' `spec_sens`/`sensitivity` threshold rules (../functions.R:158-164,
#' `scripts/Second_groundtruth/make_threshold_raster_masks.R`), except the
#' level comes from the realization's own gradient structure instead of a
#' fixed ROC-based cutoff. Applies the same >=2-pixel patch filter as the
#' deterministic masks (`scripts/Second_groundtruth/transect_placement.R:71-78`)
#' via `terra::patches()` + `terra::freq()`, so output masks are structurally
#' comparable to those.
#'
#' @param realization_path path to one `simulate_boundary_fields()` output.
#' @param n_levels number of value-bins used to locate the gradient ridge.
#' @param min_patch_px minimum patch size, in pixels, to retain (roadmap/manuscript
#' precedent: 2).
#' @return the binary inclusion mask, a `terra::rast` (1 = inside the boundary, NA = outside).
womble_boundary <- function(realization_path, n_levels = 20, min_patch_px = 2){

  r <- terra::rast(realization_path)
  grad <- spatialEco::sobal(r, method = 'intensity')

  vals <- terra::values(r)[, 1]
  grad_vals <- terra::values(grad)[, 1]
  ok <- !is.na(vals) & !is.na(grad_vals)

  bins <- cut(vals[ok], breaks = n_levels)
  mean_grad_by_bin <- tapply(grad_vals[ok], bins, mean, na.rm = TRUE)
  ridge_bin <- names(which.max(mean_grad_by_bin))
  boundary_level <- as.numeric(sub('^\\(([^,]+),.*$', '\\1', ridge_bin))

  mask <- terra::ifel(r >= boundary_level, 1, NA)

  patch_ids <- terra::patches(mask, directions = 4)
  sizes <- terra::freq(patch_ids)
  small_patches <- sizes$value[sizes$count < min_patch_px]
  if(length(small_patches) > 0){
    mask <- terra::ifel(patch_ids %in% small_patches, NA, mask)
  }
  mask
}

#' Drive the simulate + womble loop for the full boundary ensemble.
#'
#' @param pr_path,se_path paths to the RF mean/SE surfaces (`modeller()`'s `-Pr.tif`/`-SE.tif`).
#' @param holdout_sf as `estimate_residual_covariance()`.
#' @param out_dir base output directory; realizations go to `out_dir/draws`, masks to `out_dir/masks`.
#' @param n_draws,seed as `simulate_boundary_fields()`.
#' @return character vector of written mask file paths.
boundary_ensemble_run <- function(pr_path, se_path, holdout_sf, out_dir, n_draws = 500, seed = 1){

  pr_raster <- terra::rast(pr_path)
  se_raster <- terra::rast(se_path)

  cov_params <- estimate_residual_covariance(pr_raster, holdout_sf)

  draw_paths <- simulate_boundary_fields(
    pr_raster, se_raster, cov_params, n_draws = n_draws,
    out_dir = file.path(out_dir, 'draws'), seed = seed
  )

  mask_dir <- file.path(out_dir, 'masks')
  dir.create(mask_dir, recursive = TRUE, showWarnings = FALSE)

  mask_paths <- character(length(draw_paths))
  for(i in seq_along(draw_paths)){
    mask <- womble_boundary(draw_paths[i])
    f <- file.path(mask_dir, paste0('mask_', i, '.tif'))
    terra::writeRaster(mask, f, overwrite = TRUE)
    mask_paths[i] <- f
  }
  mask_paths
}
