# split_tau exploration — autonomous brief

You have the data this repo's working copy lacks (DEM rasters under `../data/spatial/processed`,
the `Data4modelling` gpkgs). Run this end to end; only stop to report at the explicit gates or
guardrails below. Working directory is `scripts/`.

## Goal (one line)
Find whether `split_tau` buys **genuine holdout replicates** for the structured CV designs
(spatial_mapping, LMPO) without wrecking spatial stringency — and if so, **measure how much a
structured split's diagnostics actually vary** when resampled. That variance is the one number
the manuscript is currently missing.

## Why (established already)
- spatial_mapping/LMPO have **no CV replicates** today: seeds 5001–5003 vary only the
  pseudo-absence draw + ranger RNG on **one cached partition**; the test set is frozen.
- We showed (manuscript §3.6, recomputed from `results_summary/all_fit_metrics.csv`) that the
  pseudo-absence **draw** is a negligible variance source (median within-cell SD ≈ 0.005 nAP),
  while the **train/test partition** dominates (random splits ≈ 0.05 nAP, 6–9× more, n-matched).
  So the variance that matters lives in the test partition — which the structured designs freeze.
- `split_tau` (in `final_sdm/Partitioner.R`, default 0 = unchanged) loosens the outer-split
  medoid picker (`kmed_pick` in `partition_split`) from the exact argmin to a seeded soft draw
  toward the k-means centroid, so each seed lands a **different but still env-central** test set.
- Open question: how much does a **structured** split vary when resampled? Can't tell yet
  (frozen). Phase 2 answers it.

## Phase 1 — probe the split (fast, no model fitting)
```
# wd = scripts/
Rscript final_sdm/prototype_split_tau.R 30 spatial_mapping
Rscript final_sdm/prototype_split_tau.R 30 LMPO
# repeat for resolutions 90 and 10 if each run is quick
```
Each writes a CSV and prints, per `tau ∈ {0, 0.15, 0.25, 0.5, 1}`: mean/min pairwise **Jaccard**
of the presence TEST sets across 8 seeds, mean **frac_beyond_p95** (env-stringency), mean n_test.

**Interpret, per (flavor, resolution):**
- If `tau = 0` already gives mean Jaccard ≲ 0.7 → kmeans-init noise alone yields replicates;
  note it — reseeding via `n_splits > 1` may suffice with no looseness. Set `tau* = 0`.
- Else pick **`tau*` = the smallest tau with mean Jaccard ≲ 0.7 AND mean frac_beyond_p95 within
  ~0.05 of its tau=0 value** (stringency preserved). Record `tau*` per (flavor, resolution).
- If no tau drops Jaccard without collapsing frac_beyond_p95 → report that; the structured test
  is too rigid to resample without losing stringency (itself a finding — stop, report).

Proceed to Phase 2 only when Phase 1 is clearly green for at least spatial_mapping @ 30.

## Phase 2 — bounded confirmatory re-fit (only if Phase 1 green)
Measure the structured-split diagnostic SD and compare it to the random (0.05 nAP) and
absence-only (0.005 nAP) SDs.
- **One resolution** (30 = 1 arc) and a **subset of ratios** spanning the range (e.g. 1:4, 1:10,
  1:20) — NOT the full 1–30 grid.
- Vary the **outer split** per replicate: build partitions at `tau*` with a **different seed per
  `k`**, `n_splits ≈ 5–8`. Today `get_partition(k)` in `First_modelling/highPA_grid.R` only
  changes the seed; with `tau* > 0` that now yields genuinely distinct splits. Thread `tau*`
  into the `partition_params(seed = split_seed_of(k))` call there (add `split_tau = tau*`).
- Write Phase-2 outputs to a **scratch tree / new cache paths** (e.g.
  `Results_spatialmapping_splittau_probe/`) so existing production results and cached partitions
  are untouched.
- Fit that small grid × the `n_splits` partitions (the only model fitting here). Then compute
  within-cell SD across the split replicates for nAP / ROC-AUC / Cox slope.

**Decision this yields:** does the spatially stringent structured split have **low** diagnostic
variance (stringency buys stability — strong result) or **high** variance like random splits
(then single structured fits are just as much of a gamble)? Either answer is publishable —
report the three-way SD comparison (structured vs random vs absence-only).

## Guardrails
- **Do NOT run the full 1:1–1:30 × all-resolutions × multi-seed sweep (~4 days).** Not needed;
  not without explicit go-ahead. Phase 2 is the bounded subset above.
- Keep `split_tau`'s default at 0. Do not overwrite production trees or their cached
  `*-partition.rds` — Phase 2 writes to a scratch tree only.
- If any step looks like it will take more than a few hours, **stop and report with a runtime
  estimate** rather than proceeding.
- No prod partition rebuilds, no git pushes, no edits to the manuscript.

## Report back
1. Phase 1: the Jaccard / frac_beyond_p95 tables (attach CSVs), chosen `tau*` per
   (flavor, resolution), and whether `tau = 0` already suffices.
2. Phase 2 (if run): structured-split SD vs random-split SD vs absence-only SD, with the
   stringency-vs-stability interpretation.
3. Any blockers (missing data, package, runtime overruns).

## Pointers
- `final_sdm/Partitioner.R` — `partition_params()` (has `split_tau`), `partition_split()`
  (`kmed_pick`), `build_partition()`.
- `final_sdm/prototype_split_tau.R` — the Phase-1 probe (header documents it).
- `First_modelling/highPA_grid.R` — production fit pipeline; `get_partition()` is the Phase-2
  wiring point.
- Context memory: `spatial-lmpo-replicates-split-tau`. Manuscript variance decomposition:
  SDManuscript.Rmd §2.8 (methods) and §3.6 (results, the `Variance Decomposition` chunk).
