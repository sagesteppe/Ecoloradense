# Roadmap: credible-interval census size estimation (H4a/H4b)

Plan for propagating real uncertainty through to a final census-size number,
instead of reporting a point estimate with no interval around it. Two
independent sources of uncertainty - population boundary and plant density -
get modelled separately, then combined by Monte Carlo composition rather than
a single joint model (see Phase 3 for why).

## Phase 0 - Infrastructure

- Install `cmdstanr` + CmdStan toolchain (not present in this environment yet).
  `torch`/GPU is **no longer needed** - superseded below: Phase 1 uses brms's
  own `posterior_epred()` for raster-scale prediction (not a hand-rolled
  matmul) and Phase 2 dropped the Bayesian-wombling route that would have
  needed it, so nothing in the implemented pipeline runs on GPU.
- Also install `spatialEco` (Sobel-Feldman gradient operator, Phase 2) - CRAN,
  already installed in the working environment this was built in.
- **Decision:** keep XGBoost (Poisson + Tweedie) and LightGBM as retained
  candidates alongside the new Bayesian candidates, rather than dropping them
  now that Tweedie has no clean Bayesian path (no native `brms` family, and a
  hurdle wrapped around Tweedie would double-model the zero process) - worth
  seeing whether the point-estimate ML models actually beat the Bayesian
  candidates rather than assuming they won't. The Phase 1 model-comparison
  table is therefore the existing five candidates plus the four Bayesian
  additions.

## Phase 1 - Density: Bayesian candidate models

The field data is integer counts per 3x3m quadrat (individuals by maturity
stage), so the natural Bayesian candidates mirror what a count-data GLM
family offers, run the same way the existing XGBoost step already compares
multiple distributions rather than committing to one up front:

- `poisson()` - direct Bayesian analogue of the existing XGBoost-Poisson
  candidate. Native, trivial.
- `negbinomial()` - handles overdispersion relative to Poisson, which
  clustered/mat-forming species like this one often show. Native.
- `hurdle_poisson()` - splits "quadrat has any plants at all" from "how many,
  given present" - useful if zeros are more common than a plain Poisson/NB
  would predict. Native.
- `hurdle_negbinomial()` - both excess zeros and overdispersion in one
  family; likely the single most defensible choice for this data, and the
  natural default if only one Bayesian candidate is fit. Native.

Implementation notes:
- Fit at the plot level, same covariates/n as the existing XGBoost/LightGBM
  candidates - not raster-cell level.
- `gp(x, y)` spatial term (approximate/Hilbert-space, `k`/`c` args - keeps the
  posterior a plain coefficient+basis-weight matrix), `backend = "cmdstanr"`.
- Spatial CV: **not** a literal reuse of `splitData()`'s `nndm_indices` -
  confirmed by testing that `CAST::nndm()` returns near-leave-one-out folds
  (fold count = nrow(train)), fine for the cheap XGBoost/ranger fits but not
  for repeated Stan refits. Uses a fresh fixed-`k` `CAST::knndm(k=10)` call
  instead (the same mechanism the RF suitability model already uses,
  `functions.R:103`), applied to the density training points. Non-spatial
  mode is an ordinary random `K=10` partition via `brms::kfold()`, not
  literally `rsample::bootstraps()` (`brms::kfold()`'s disjoint-fold
  semantics don't fit a with-replacement resample). Each of the 4 families is
  compared under both modes (per an earlier decision to match the XGBoost
  Poisson/Tweedie precedent) - 8 comparison fits.
- Apply the same cheap-search/expensive-final-fit split already used in
  `adaptive_PAratio_search()`: fewer iterations during CV comparison, full
  precision only on the promoted model.
- Output: the promoted model's posterior draws (a plain coefficient matrix) -
  a standard 4-chain fit gives ~4000 for free.
- Implemented in `scripts/Census_uncertainty/density_bayes.R`
  (`brms_count_model()`, `brms_spatial_fold_ids()`, `brms_cv_compare()`,
  `brms_promote_and_refit()`), wired into `densityModeller()`
  (`functions.R:1184-1198`) after the existing LGBM block. Two pre-existing
  bugs had to be fixed first for `densityModeller()` to be reachable at all:
  `wrapper()`'s unconditional early return (`functions.R:1253-1290`, now
  gated behind `return_early`) and `splitData()`'s expectation of
  `Longitude`/`Latitude` columns that `wrapper()`'s sf-geometry-only output
  didn't have (now added in `wrapper()` before the `densityModeller()` call).

## Phase 2 - Boundary: simulate, don't refit

- No Bayesian refit of suitability - reuse the existing random forest mean
  surface + jackknife SE surface as-is.
- Estimate the spatial correlation structure of the RF residuals (needed to
  make the simulation spatially coherent, not just per-cell noise).
- Generate boundary-field realizations via conditional Gaussian simulation,
  at true 1/3-arc resolution - no coarsening, since resampling the
  uncertainty ensemble down would undercut the resolution-comparison thesis
  the rest of the manuscript is testing (H2).
  **Do not do this via a dense covariance matrix + Cholesky, even on GPU** -
  a dense covariance matrix at full 1/3-arc scale is on the order of 1 TB
  even for an aggressively-masked domain (~500k cells -> n^2 entries), which
  is a memory wall no GPU fixes. Use FFT/circulant-embedding simulation
  instead (`fields::circulantEmbedding`, the package's standard tool for
  simulating a stationary/isotropic random field on a big regular grid,
  O(n log n) instead of O(n^2) memory) - this runs on CPU, not the GPU
  stack; the GPU budget in this roadmap is for the density-prediction
  matmuls (Phase 3), not boundary simulation.
- Womble each realization - i.e. detect the boundary as the ridge of highest
  gradient magnitude in the simulated field. **Not** Bayesian wombling in the
  end, despite real effort spent trying: `github.com/arh926/spWombling`
  (Halder, Banerjee & Dey 2023, *JASA*, extending Banerjee & Gelfand 2006's
  first-order wombling with covariates and a spatial random effect) turned
  out to have real bugs on the pinned commit
  (`05f43e92d97f811bea2c52081ea13dbe6c2aa3b0`) - `spatial_gradient()`
  referenced undeclared arguments (100% failure), and `bayes_cwomb()` had an
  `ntimes` off-by-vector bug plus an `apply()`/`t()` orientation bug in its
  matern1 branch. Two of these got diagnosed, fixed, validated (by
  reproducing the paper's own Section 5 coverage check), and PR'd upstream -
  **github.com/arh926/spWombling/pull/3** - but a further, unresolved
  numerical issue (`eigen()` hitting non-finite values a few steps later)
  combined with an inactive-looking maintainer made depending on it
  impractical for this pipeline. Fell back to a non-Bayesian but tested,
  citable implementation instead: `spatialEco::sobal(r, method='intensity')`
  (Sobel & Feldman 1969) for the gradient-magnitude step, thresholded at the
  realization's own gradient-ridge value (a data-driven per-draw analogue of
  the deterministic `spec_sens`/`sensitivity` rules), with the same >=2-pixel
  patch filter the deterministic masks already use.
  The trade-off this accepts: draw-to-draw *variability* in the boundary
  (from the circulant-embedding ensemble) still flows through to Phase 3's
  final census-size interval, but there's no direct posterior statement
  about the boundary curves themselves (no "this edge is credibly real, here,
  at X% probability") the way `bayes_cwomb()` on the two candidate curves
  would have given - that would need revisiting if the maintainer merges/
  fixes the remaining issue, or if it's worth patching further ourselves.
- Start the budget at ~500 draws; this is the expensive stack, so it sets the
  pace rather than matching the density stack's free ~4000.
- Implemented in `scripts/Census_uncertainty/boundary_simulation.R`
  (`estimate_residual_covariance()`, `simulate_boundary_fields()`,
  `womble_boundary()`, `boundary_ensemble_run()`). Two real numerical
  gotchas surfaced and got fixed during testing, not just theoretical risks:
  `fields::circulantEmbeddingSetup()`'s default padding (~2x grid dimension)
  can leave small negative spectral weights from floating-point error rather
  than genuine non-positive-definiteness - `simulate_boundary_fields()` now
  retries with growing padding (2x/4x/8x/16x/32x) before giving up; and
  `fields::stationary.cov()` resolves its `Covariance` argument by name
  lookup on the search path (not via `fields::`-qualification), so the
  namespace has to actually be attached, and there is no function literally
  named `"Gaussian"` in `fields` (its squared-exponential covariance,
  `gauss.cov()`, uses an incompatible coordinate-based calling convention) -
  gstat `"Gau"` model fits fall back to Exponential.

## Phase 3 - Combine

The boundary model and the density model are fit independently with no
shared MCMC chain, so their draws can't be paired index-for-index as if
correlated - but as long as the two uncertainty sources are genuinely
independent (reasonable here: different data, different models), they can be
validly combined by Monte Carlo composition instead of one joint multivariate
model. A joint model would be "more correct" in principle but reopens the
spatial-confounding problem from early planning discussion, for a likely
marginal gain in interval accuracy.

- Mask tightly first to regional bounding boxes used in manuscript. 
- For each of ~1000-2000 combined draws: pair one boundary mask with one
  density posterior-predictive raster (resample the smaller boundary stack
  with replacement against the larger density stack, not the other way
  around - don't throw away density draws you already paid for), multiply by
  cell area, sum to one census-size value.
- Predict -> use -> discard per draw: `brms::posterior_epred()` for the
  density-raster prediction step (not a hand-rolled GPU matmul - settled
  earlier as a design fork, since brms's own machinery is correctness-first
  and the raster-scale cost turned out manageable without it), write into
  `terra`, discard, move to the next draw. Never materialize all N layers
  simultaneously.
- Track the running 2.5/97.5 percentile of the combined census-size draws as
  N grows; stop once it stops moving instead of committing to a fixed N up
  front. Stan's own `ess_tail` diagnostic uses ~400 as an absolute floor for
  trusting a tail quantile; 1000-2000 is the comfortable zone.
- **Confirmed by testing, easy to get wrong**: `posterior_epred()` for a
  `gp()` model at new (out-of-sample) locations draws the GP's latent value
  from R's global RNG on *every call* - not reproducible given only
  `draw_ids`, and the draw is a joint function of the whole batch of new
  locations, so predicting the same cells inside a differently-sized batch
  does not reproduce the same values either. Two implications, both handled
  in `predict_density_draw()`: it seeds on `draw_id` (so a reused density
  draw - `pair_draws()` resamples with replacement - reproduces identically,
  as the Monte Carlo composition requires), and it does **not** chunk
  `newdata` across multiple `posterior_epred()` calls for a `gp()` model
  (chunking would make adjacent chunks independently-resampled fragments,
  not consistent pieces of one spatially-coherent draw).
- Implemented in `scripts/Census_uncertainty/combine_census.R`
  (`pair_draws()`, `predict_density_draw()`, `align_mask_to_template()`,
  `census_from_pair()`, `combined_census_montecarlo()`). Also fixed in
  testing: Phase 2's masks are built on the RF suitability surface's full
  extent while the density side works on a region-cropped raster, so every
  mask needs `terra::resample(..., method='near')` onto the density raster's
  exact (post-crop) geometry before `terra::mask()` - `align_mask_to_template()`
  does this once per mask, not once per draw, since masks get reused often
  under with-replacement resampling; and `terra::global(x, fun='sum',
  na.rm=TRUE)` returns `NA` (not `0`, unlike base R's `sum()`) when every
  cell is `NA` - a draw whose mask excludes everything is a genuine "0
  individuals" outcome and must not silently drop out of the draws vector.

## Phase 4 - Validate before trusting it

- Small-scale numerical check before scaling up - not CPU-vs-GPU anymore
  (moot once Phase 3 settled on brms's own prediction, not a hand-rolled
  matmul), but two checks that are actually meaningful given the gp()-resampling
  behavior above: two full-raster `predict_density_draw()` calls with the same
  `draw_id` must reproduce identically (required for Phase 3's with-replacement
  resampling), and its output must match an independent same-seed,
  same-full-`newdata` `posterior_epred()` call read back at the same cells
  (the actual reshaping/indexing check).
- Sanity-check the boundary ensemble's median against the two existing
  deterministic threshold rules (`spec_sens`/`sensitivity` from Section 2.7)
  - they should roughly bracket it, not disagree wildly.
- Re-check the spatial confounding question (suitability-as-covariate vs. a
  spatial random effect over the same coordinates) once the density model's
  `gp()` term is actually fit, not just in principle.
- Implemented in `scripts/Census_uncertainty/validate_pipeline.R`
  (`check_predict_density_draw()`, `sanity_check_boundary_median()`,
  `check_spatial_confounding()`).

## Phase 5 - Write-up

- Sections 2.6/2.8 methods text for the Bayesian density candidates and the
  simulate-womble-combine boundary pipeline.
- Swap point-estimate density/boundary figures for versions with credible
  bands.
