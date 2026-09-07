# Roadmap: population-blocked split + sample-size rarefaction (H1/H2 robustness)

Motivation: two related confounds surfaced while reviewing the ground-truth
evaluation and planning the H2 resolution comparison. (1) `spThin`'s thinning
distance is the cell hypotenuse, so it's resolution-dependent - finer
resolutions mechanically retain more training points from the same raw
occurrence pool (already flagged in-line in `SDManuscript.Rmd`: "this must be
wrong... at least ~3.8x were generated for finest resolution"). Any H2
comparison right now conflates "resolution helps" with "resolution N is
bigger." (2) Separately, heavily-documented populations may dominate the
niche estimate purely from historical survey effort, independent of
resolution or model choice. This roadmap adds a population-identity field to
the occurrence data and four companion analyses to partition these effects
out. Nothing below is implemented yet - this is the design handoff.

## Phase 0 - Population identity [implemented]

- `data/collections/occurrences_coloradense/occurrences.shp` (the iNaturalist
  + herbaria records used for ENM fitting) carries only `eventDate` +
  `geometry` - no population/site ID.
- `data/GroundTruthPts/ground_truth_median-VISITED.gpkg`'s `areas` layer has
  27 named population/site polygons (Cochetopa, Antero, Bellview, Gothic,
  etc.) - the natural join key via point-in-polygon.
- **Coverage check (done)**: point-in-polygon against `areas`, checked
  against the actual `Data4modelling/*-presence-iterN.gpkg` pools (post-
  thinning, so this is coverage of what's really used for fitting, not the
  raw occurrence pool): iter1 ranges 70.2% (90m) - 75.0% (3m) inside a
  polygon; iter0 ranges 66.7% (90m) - 74.5% (3m). Coverage is incomplete but
  substantial, and improves slightly at finer resolutions (more retained
  points overall, same polygon footprint).
- **Fallback (decided)**: `assign_population()` in `functions.R`. Primary
  assignment is point-in-polygon. For points outside every polygon, checked
  the distance-to-nearest-polygon distribution (iter1 3m pool): ~55% of
  unassigned points sit within 500m, ~71% within 1km, but the tail runs to
  ~27km - those are real, distinct, unsurveyed populations, not GPS/boundary
  slop, so annexing them to whichever named site is nearest would corrupt
  the population-blocked split's "held-out population" framing. Landed on a
  conservative 500m nearest-polygon fallback for the strict `population`
  field (used wherever true membership matters - Phase 2/3's per-population
  rarefaction); points beyond that are `NA` there and should be dropped from
  those population-level analyses only, per the original plan, not guessed
  at. A second field, `population_zone` (always the single nearest `Site`,
  no `NA`, no distance cutoff), is a separate full-landscape partition used
  only as a spatial blocking key by Phase 1's split - see below - not a
  membership claim.

## Phase 1 - Leave-one-population-out CV (new iter1 variant) [implemented, not yet run]

- The existing "classic" iter1 split partitions individual *points*
  (twin/nndm-based, i.e. `spatial_class_split()` in `functions.R`). Added a
  second variant that blocks by *population*: whole populations assigned
  entirely to train or entirely to test. This tests transferability to a
  population the model has never partially seen, not just interpolation
  within populations it has some data from - closer in spirit to Borokini et
  al.'s spatial-block CV, but blocked by real population identity rather
  than an arbitrary geographic grid.
- **Design turn**: originally implemented as a single random population-
  block draw (whole populations shuffled and greedily packed into a test
  bin until the target test proportion of presence records was hit),
  reusing the same "compute the split once, reuse across every ratio/seed
  replicate" pattern as `spatial_class_split()`. Dropped in favor of full
  leave-one-population-out: a single random draw only characterizes model/
  seed stochasticity for whichever populations happened to land in test on
  that one draw, not variability from a different, equally arbitrary draw -
  which is the thing this robustness check is actually supposed to speak
  to. Full LOPO removes that source of variance instead of estimating the
  risk of a bad one, and - relative to what `adaptive_PAratio_search()`
  already spends per resolution (Boruta feature selection plus a knndm-
  tuned hyperparameter grid, repeated across dozens of ratio/seed
  replicates) - fitting ~15-20 more single models, one per population, is a
  drop in the bucket, not a meaningful multiplier on the expensive part of
  the pipeline. So the earlier "check compute budget before committing to
  LOPO" concern doesn't actually bind here.
- **Implementation**: `population_lopo_search()` in `functions.R`. Does NOT
  re-run the PA-ratio bracketing search per population - it takes a single
  fixed `PAratio` (in practice, whichever `adaptive_PAratio_search()` already
  chose as `best_ratio` for that resolution/iteration) and fits it once per
  population, with that population's zone (`population_zone` from
  `assign_population()` - the nearest named population regardless of
  distance, so every point, presence or background, lands in exactly one
  zone) held out entirely as test, every other point as train. Populations
  with fewer than `min_presence` (default 3) presence records don't get a
  fold; populations whose zone happens to catch only one class in test get
  `NA` metrics rather than a fit, since PR-AUC/ROC-AUC need both classes -
  checked against the real iter1 90m pool, this affects a small edge zone
  ("13043 Ridge-Long-Cassi Valley", 0 absences in its test fold) out of 15
  populations that otherwise qualify. Each fold's `modeller()` fit gets its
  own seed (`base_seed + fold index`) since `modeller()`'s cache file stem
  doesn't include the held-out population - sharing a seed across folds
  would silently reload/reuse the first fold's fit for every later one.
  `Modelling_iter0.Rmd` and `Modelling.Rmd` now each add an LOPO sweep per
  resolution (`lopo_90`/`lopo_30`/`lopo_10`, using `search_XX$best_ratio`),
  writing to a new `results_populationsplit/` tree (mirrors
  `results_classicsplit/`'s models/modelsTune/tables/test_data/evaluations/
  suitability_maps layout) and a per-resolution
  `<res>-Iteration<iteration>-LOPO-PA<ratio>DO:<distOrder>.csv` per-population
  table, plus a `-summary.csv` (see next point).
- **Population size imbalance (on eval)**: populations range from n=1
  presence record up to n=50 (Cochetopa) in the real iter1 90m pool - an
  unweighted mean PR-AUC/ROC-AUC across populations would let a 3-point
  fold (which, with that few positives, can only take a handful of discrete
  PR-AUC values - essentially a coin flip) count the same as a genuinely
  well-estimated 50-point fold. `population_lopo_search()` returns
  `list(per_population, summary)`; `summary` comes from the new
  `summarize_lopo()` helper, which weights each population's contribution
  by its own `n_presence_test` (dropping any `NA`/zero-class rows first).
  Sanity-checked against synthetic fold sizes matching the real population
  distribution: unweighted mean 0.756 vs. presence-weighted 0.866 for the
  same folds - the unweighted version was letting one noisy 3-point fold
  pull the average down far more than its actual evidentiary weight
  justifies. Any manuscript figure/number built on top of this table should
  use `summary`, not a raw `mean()` of the per-population column - or, per
  discussion, at minimum plot metric-vs-n per population rather than
  collapsing to a single number, so the reliability gradient is visible
  rather than hidden.
- **Not yet run**: the actual model fits (Boruta + knndm CV + ranger tuning)
  weren't executed this session - this machine's R install is missing
  `ranger`, `Boruta`, and `gdalUtilities` (this looks like a lighter/dev R
  environment than wherever `PROJ_ROOT` in `functions.R` points, which is
  presumably where the real fits run). The fold-construction logic itself
  (population join, per-population train/test id assignment, the zero-class
  guard) was smoke-tested against the real iter1 90m pool and looks correct;
  run `Modelling_iter0.Rmd` and `Modelling.Rmd`'s new `lopo_XX` chunks on
  the modeling machine for the actual fits.
- **Manuscript note**: state explicitly that this population-blocked check
  was not part of the original sampling design - it's a robustness check
  added after the fact, not pre-registered.
- The 2026 ground-truth field evaluation
  (`Second_groundtruth/digitize_evaluate_2026_groundtruth.R`) rolls over
  unchanged onto this variant: score both the classic-split model and the
  population-split model against the same held-out field points, so the
  manuscript can report whether field-validated performance differs between
  a model validated by interpolation vs. one validated by
  extrapolation-to-new-population.

## Phase 2 - Within-population rarefaction (both components, per decision)

1. **Environmental coverage (cheap, no refitting)** - for each population,
   add points one at a time (random order, repeated draws for a
   distribution) and track how much of that population's own covariate range
   gets captured - e.g. convex-hull or PCA-space volume across the 50
   environmental variables (§2.2), normalized against the volume using all
   of that population's available points. Answers "how many points does it
   take to represent this population's variability" without touching the
   model. Run this first - it's fast and may show some populations are
   already saturated at low N, narrowing which populations need step 2.
2. **Downstream model performance (expensive, real refits)** - for
   populations still under-sampled after step 1, refit with N points from
   that population (all other populations' data held fixed) and score
   against the held-out 2026 field points at increasing N, to find where the
   field-validated AUC-PR curve for that population's contribution flattens.

## Phase 3 - Species-level rarefaction (learning curve)

- Same idea as Phase 2.2 but pooled across the whole species: vary total
  occurrence N (and/or number of populations represented) in the training
  set, refit, score against held-out field points, plot the learning curve.
  Answers "how many points/populations are needed to model the species" as
  opposed to any one population.
- Natural to extend across resolutions once the 1 and 1/3 arc-second
  `-Pr.tif` surfaces exist on this machine (only 3 arc-second is built here
  currently).

## Phase 4 - Resolution N-confound (H2 robustness - first priority)

- Rarefy all three resolutions down to a common matched N (the smallest,
  i.e. whatever the 3 arc-second thinned pool provides) before comparing
  performance, since `spThin`'s radius is resolution-dependent and
  mechanically inflates finer-resolution N.
- Repeat the subsample-and-refit many times (bootstrap-style), report the
  performance *distribution* per resolution at matched N - not a single
  point estimate from one draw.
- Supplementary to the main full-data H2 comparison, not a replacement for
  it - don't cripple the primary models just to answer this robustness
  question.

## Why these fixes are needed (reviewer-facing framing)

Both Zimmer et al. 2023 and Borokini et al. 2023 fold newly-collected field
points back into training and report performance changes without a
sample-size-matched or population-blocked control, which leaves open whether
reported gains reflect genuine generalization or just more data / more
spatially autocorrelated neighbors of existing training points. Phases 1 and
4 here are the direct fix: population-blocked evaluation tests
extrapolation, not interpolation, and matched-N rarefaction isolates
resolution/data effects from raw sample-size effects.
