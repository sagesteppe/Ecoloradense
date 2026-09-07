## Driver for census_uncertainty_roadmap.md's Phases 1-4. Sources the
## phase-specific scripts in order (each documents which functions it
## contributes) and runs them for one resolution/PAratio product at a time -
## the same "one bn per run" granularity `densityModeller()`/`wrapper()`
## already use elsewhere in this repo.
##
## This is a template, not a push-button script: the paths below are the
## repo's existing conventions, but region_bbox, the promoted family/cov_type,
## and the draw budgets (n_draws/n_max) are per-product decisions the roadmap
## leaves to be set at run time, not hardcoded here.
##
## Run with the working directory set to scripts/ - the same convention every
## other script/Rmd in this repo already assumes (relative paths like
## '../results/...' throughout functions.R).

source('functions.R')
source('Census_uncertainty/density_bayes.R')
source('Census_uncertainty/boundary_simulation.R')
source('Census_uncertainty/combine_census.R')
source('Census_uncertainty/validate_pipeline.R')

#' Run Phases 1-4 for one resolution/PAratio product.
#'
#' @param x the same `-Pr.tif` filename `wrapper()` takes elsewhere in this
#' repo (e.g. `'1-3arc-Iteration1-PA1_3.3DO_0-Pr.tif'`).
#' @param holdout_sf sf POINT object with an `Occurrence` column, for Phase 2's
#' residual-covariance estimate.
#' @param region_bbox a `terra::ext` to mask the combined draws to (Phase 3).
#' @param spec_sens_mask_path,sensitivity_mask_path Section 2.7's deterministic masks.
#' @param n_boundary_draws,n_combined_max as `boundary_ensemble_run()`/`combined_census_montecarlo()`.
#' @param fp,out_dir base paths for density-model and boundary-ensemble outputs.
#' @return list(density = densityModeller()'s return, boundary_masks =
#' character vector, census = combined_census_montecarlo()'s return,
#' validation = list of the three Phase 4 checks' results).
run_census_uncertainty <- function(x, holdout_sf, region_bbox,
                                    spec_sens_mask_path, sensitivity_mask_path,
                                    n_boundary_draws = 500, n_combined_max = 2000,
                                    fp = file.path('..', 'results', 'CountModels'),
                                    out_dir = file.path('..', 'results', 'boundary_ensemble')){

  bn <- gsub('DO.*$', '', x)

  # Phase 1 - density (fits the 5 existing ML candidates + 4 new Bayesian
  # ones, promotes the winner; wrapper() adds Longitude/Latitude first).
  df <- wrapper(x, return_early = FALSE)
  density_result <- densityModeller(df, bn = bn, fp = fp)

  # Phase 2 - boundary ensemble, from the RF suitability model's own outputs.
  pr_path <- file.path('..', 'results', 'suitability_maps', paste0(bn, '-Pr.tif'))
  se_path <- file.path('..', 'results', 'suitability_maps', paste0(bn, '-SE.tif'))
  boundary_masks <- boundary_ensemble_run(pr_path, se_path, holdout_sf,
                                           out_dir = out_dir, n_draws = n_boundary_draws)

  # Phase 4a/4b - validate before combining: reshape/reproducibility check on
  # the promoted density model, and the boundary ensemble's plausibility.
  check_predict_density_draw(density_result$brms_promoted$Model,
                              terra::rast(pr_path), draw_id = 1)
  boundary_check <- sanity_check_boundary_median(boundary_masks, spec_sens_mask_path,
                                                  sensitivity_mask_path)

  # Phase 3 - combine.
  census_out_csv <- file.path('..', 'results', 'tables', paste0(bn, '-census_size_posterior.csv'))
  census <- combined_census_montecarlo(
    density_result$brms_promoted$Model, terra::rast(pr_path), boundary_masks,
    region_bbox, n_max = n_combined_max, out_csv = census_out_csv
  )

  # Phase 4c - spatial confounding re-check, now that gp() is actually fit.
  confounding_check <- check_spatial_confounding(
    density_result$brms_promoted$Model, df, family = poisson()  # match the promoted family
  )

  list(density = density_result, boundary_masks = boundary_masks, census = census,
       validation = list(boundary = boundary_check, confounding = confounding_check))
}
