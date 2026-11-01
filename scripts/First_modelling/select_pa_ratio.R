## Choose the presence:absence ratio per resolution and CV design, from HOLDOUT (test) data
## only - the 2026 field ground truth is kept out of selection entirely, so it stays an
## independent test of the chosen models.
##
## Reworked copy of select_far_ratio.R. Key differences:
##   * reads the already-scored per-replicate metrics (all_fit_metrics.csv, the same rows the
##     manuscript PA-ratio grid uses) instead of re-predicting from the stored models, so no
##     sf/terra/ranger and no field-absence bookkeeping.
##   * three CV designs: Classic split, Spatial mapping, LMPO (Spatial kNN dropped).
##   * candidate ratios capped at 1:30.
##   * evaluation metrics: nAP, ROC-AUC, Cox slope, Cox intercept, WEIGHTED EQUALLY in a
##     composite. The Cox metrics are scored TWO-SIDED - distance to the target (slope -> 1,
##     intercept -> 0) regardless of side. A holistic composite is used because the gap
##     between test and ground-truth behaviour means no single metric should drive the choice.
##   * prevalence is carried inside nAP via the ratio-implied baseline pi = 1/(1+ratio); set
##     fixed_prevalence <- 0.05 to score nAP at a fixed landscape prevalence instead.
##
## Composite, per resolution and design (fixed before looking at any results):
##   distance    - each metric's per-replicate distance to its target: nAP -> |nAP - 1|,
##                 ROC-AUC -> |ROC - 1|, Cox slope -> |slope - 1|, Cox intercept -> |intercept|
##   normalize   - the four metrics live on different scales, so each metric's distance is
##                 min-max scaled to [0, 1] ACROSS the candidate ratios within a
##                 (design, resolution) group, so one metric's wider spread can't dominate
##   composite   - the equal-weight mean of the four normalized distances (lower is better)
##   one-SE rule - the SMALLEST ratio whose mean composite is within one standard error
##                 (across replicates) of the best (lowest) mean composite
##   model       - the replicate at that ratio whose composite is closest to the median
##
## Usage (working directory scripts/): Rscript First_modelling/select_pa_ratio.R

suppressPackageStartupMessages(library(tidyverse))
source('functions.R')  # sets PROJ_ROOT

## --- parameters -------------------------------------------------------------------------
max_ratio          <- 30      # candidate ratios: PAratio <= max_ratio
nAP_target         <- 1
roc_target         <- 1
cox_slope_target   <- 1
cox_intercept_target <- 0
weights            <- c(nAP = 1, roc = 1, slope = 1, intercept = 1)  # equal; edit to reweight
fixed_prevalence   <- NA      # NA -> nAP uses the ratio-implied prevalence 1/(1+ratio);
                              # set 0.05 to score nAP at a fixed landscape prevalence instead

design_cols <- c(
  Results_random_iter1a_high_pa         = 'classic',
  Results_spatialmapping_iter1a_high_pa = 'spatial_mapping',
  Results_LMPO_iter1a_high_pa           = 'lmpo')

out_dir <- file.path(PROJ_ROOT, 'results_summary')

## nAP: prevalence-normalized average precision. AP = pr_auc; a no-skill classifier scores
## AP = pi, so nAP rescales onto [baseline -> 0, perfect -> 1]. pi defaults to the fit's own
## holdout prevalence 1/(1+ratio); fixed_prevalence overrides it with a landscape value.
nAP <- function(ap, ratio, prevalence = NA) {
  pi <- if (is.na(prevalence)) 1 / (1 + ratio) else prevalence
  (ap - pi) / (1 - pi)
}

## --- 1. per-replicate metrics on test data ----------------------------------------------
raw <- read.csv(file.path(out_dir, 'all_fit_metrics.csv'))
raw <- raw[raw$tree %in% names(design_cols) & !is.na(raw$PAratio) & raw$PAratio <= max_ratio, ]

reps <- tibble(
  design      = unname(design_cols[raw$tree]),
  resolution  = raw$resolution,
  ratio       = raw$PAratio,
  replicate   = raw$seed,
  nAP         = nAP(raw$pr_auc, raw$PAratio, fixed_prevalence),
  roc_auc     = raw$roc_auc,
  cox_slope   = raw$cox_slope,
  cox_intercept = raw$cox_intercept) |>
  ## per-replicate distance of each metric to its target (Cox distances are two-sided)
  mutate(d_nAP   = abs(nAP       - nAP_target),
         d_roc   = abs(roc_auc   - roc_target),
         d_slope = abs(cox_slope - cox_slope_target),
         d_int   = abs(cox_intercept - cox_intercept_target)) |>
  filter(is.finite(d_nAP), is.finite(d_roc), is.finite(d_slope), is.finite(d_int))

reps$design <- factor(reps$design, levels = unname(design_cols))

## --- 2. equal-weight composite ----------------------------------------------------------
## Per-ratio mean distance for each metric; these set the min-max scaling range within each
## (design, resolution) group. min-max is linear, so normalizing per replicate and averaging
## gives the same per-ratio mean as normalizing the per-ratio means - but also a valid SE.
per_ratio_d <- reps |>
  group_by(design, resolution, ratio) |>
  summarise(across(c(d_nAP, d_roc, d_slope, d_int), ~ mean(.x, na.rm = TRUE)),
            .groups = 'drop')

# group min/max of each metric's per-ratio mean distance -> the scaling constants
scaling <- per_ratio_d |>
  group_by(design, resolution) |>
  summarise(across(c(d_nAP, d_roc, d_slope, d_int),
                   list(mn = ~ min(.x, na.rm = TRUE), mx = ~ max(.x, na.rm = TRUE))),
            .groups = 'drop')

norm01 <- function(x, mn, mx) ifelse(mx > mn, (x - mn) / (mx - mn), 0)  # vectorized; flat metric -> 0
w <- weights / sum(weights)

reps <- reps |>
  left_join(scaling, by = c('design', 'resolution')) |>
  mutate(composite = w['nAP']       * norm01(d_nAP,   d_nAP_mn,   d_nAP_mx) +
                     w['roc']       * norm01(d_roc,   d_roc_mn,   d_roc_mx) +
                     w['slope']     * norm01(d_slope, d_slope_mn, d_slope_mx) +
                     w['intercept'] * norm01(d_int,   d_int_mn,   d_int_mx))

## per-ratio composite (lower is better) + the raw metric means for context
per_ratio <- reps |>
  group_by(design, resolution, ratio) |>
  summarise(n = n(),
            composite_mean = mean(composite, na.rm = TRUE),
            composite_se   = sd(composite, na.rm = TRUE) / sqrt(sum(!is.na(composite))),
            nAP_mean       = mean(nAP, na.rm = TRUE),
            roc_mean       = mean(roc_auc, na.rm = TRUE),
            cox_slope_med  = median(cox_slope, na.rm = TRUE),
            cox_int_med    = median(cox_intercept, na.rm = TRUE),
            .groups = 'drop')

## --- 3. selection: one-SE rule on the composite -----------------------------------------
## best = lowest mean composite; the one-SE window is best + one SE, and the SMALLEST ratio
## inside it is chosen (favour a leaner PA ratio when the composite is statistically tied).
select_ratio <- function(d){
  best <- d[which.min(d$composite_mean), ]
  best_se0 <- ifelse(is.na(best$composite_se), 0, best$composite_se)  # single replicate -> 0
  within_1se <- filter(d, composite_mean <= best$composite_mean + best_se0)
  tibble(best_ratio = best$ratio, best_composite = best$composite_mean, best_se = best$composite_se,
         chosen_ratio = min(within_1se$ratio))
}

selection <- per_ratio |>
  group_by(design, resolution) |>
  group_modify(~ select_ratio(.x)) |>
  ungroup()

## median-performing replicate at the chosen ratio (closest composite to the median)
median_model <- selection |>
  filter(!is.na(chosen_ratio)) |>
  pmap_dfr(function(design, resolution, chosen_ratio, ...){
    r <- filter(reps, design == !!design, resolution == !!resolution, ratio == !!chosen_ratio)
    r[which.min(abs(r$composite - median(r$composite, na.rm = TRUE))), ] |>
      transmute(design, resolution, chosen_ratio = ratio, median_replicate = replicate,
                composite, nAP, roc_auc, cox_slope, cox_intercept)
  })

## --- 4. write + report ------------------------------------------------------------------
write.csv(per_ratio,    file.path(out_dir, 'pa-ratio-per-ratio-scores.csv'), row.names = FALSE)
write.csv(selection,    file.path(out_dir, 'pa-ratio-selection.csv'),        row.names = FALSE)
write.csv(median_model, file.path(out_dir, 'pa-ratio-median-models.csv'),    row.names = FALSE)

options(width = 200)
message('\n=== chosen PA ratio (one-SE rule on the equal-weight composite of 4 metrics) ===')
selection |>
  left_join(select(per_ratio, design, resolution, ratio, nAP_mean, roc_mean,
                   cox_slope_med, cox_int_med), by = c('design', 'resolution', 'chosen_ratio' = 'ratio')) |>
  select(design, resolution, chosen_ratio, best_ratio,
         nAP_mean, roc_mean, cox_slope_med, cox_int_med) |>
  arrange(design, resolution) |>
  as.data.frame() |> print(digits = 3, row.names = FALSE)
message('\n=== median models at the chosen ratio ===')
print(as.data.frame(median_model), digits = 3, row.names = FALSE)
