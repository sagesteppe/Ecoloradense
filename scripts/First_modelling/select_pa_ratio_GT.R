## Counterfactual PA-ratio selection run on the 2026 FIELD GROUND TRUTH instead of the holdout
## test split, SEPARATELY per interpolation domain (<=270 m and >270 m). Mirrors
## select_pa_ratio.R exactly - same equal-weight composite of 4 metrics (nAP, ROC-AUC, Cox
## slope, Cox intercept, two-sided Cox distances), same min-max normalization within group, and
## the same 5% rule - but the "group" now also includes `domain`, so each design x resolution
## picks a ratio per domain. This asks: if we had tuned the PA ratio on the field data itself,
## which ratio would each design choose, and does it match the holdout-chosen ratio?
##
## Differences from select_pa_ratio.R:
##   * reads all_fit_2026gt.csv (field) not all_fit_metrics.csv (holdout).
##   * nAP baseline is the domain's field prevalence (stored per row), not 1/(1+ratio).
##   * Extrapolation is dropped (selection is an interpolation question); calibration = raw to
##     stay comparable with the raw holdout selection (set CALIB <- 'calibrated' to switch).
##   * spatial/LMPO use the tau split-replicate trees (phase2, 1arc only); classic = random.
##
## Usage: Rscript First_modelling/select_pa_ratio_GT.R

suppressPackageStartupMessages(library(tidyverse))
PROJ_ROOT <- '/home/sagesteppe/Documents/Ecoloradense'
out_dir   <- file.path(PROJ_ROOT, 'results_summary')

## --- parameters (match select_pa_ratio.R) -----------------------------------------------
max_ratio <- 30
rel_tol   <- 0.05                       # smallest ratio within 5% of the best mean composite
CALIB     <- 'calibrated'               # 'raw' (comparable to holdout) or 'calibrated'
weights   <- c(nAP = 1, roc = 1, slope = 1, intercept = 1)
targets   <- c(nAP = 1, roc = 1, slope = 1, intercept = 0)

design_cols <- c(
  Results_random_iter1a_high_pa          = 'classic',
  Results_spatialmapping_splittau_phase2 = 'spatial_mapping',
  Results_LMPO_splittau_phase2           = 'lmpo')

## --- 1. per-replicate field metrics, per domain -----------------------------------------
g <- read.csv(file.path(out_dir, 'all_fit_2026gt.csv'))
g <- g[g$tree %in% names(design_cols) & g$calibration == CALIB &
       g$domain != 'Extrapolation' & !is.na(g$PAratio) & g$PAratio <= max_ratio, ]
g$nAP <- (g$pr_auc - g$prevalence) / (1 - g$prevalence)   # field-prevalence baseline

reps <- tibble(
  design = unname(design_cols[g$tree]), resolution = g$resolution, domain = g$domain,
  ratio = g$PAratio, replicate = g$seed,
  nAP = g$nAP, roc_auc = g$roc_auc, cox_slope = g$cox_slope, cox_intercept = g$cox_intercept) |>
  mutate(d_nAP   = abs(nAP           - targets['nAP']),
         d_roc   = abs(roc_auc       - targets['roc']),
         d_slope = abs(cox_slope     - targets['slope']),
         d_int   = abs(cox_intercept - targets['intercept'])) |>
  filter(is.finite(d_nAP), is.finite(d_roc), is.finite(d_slope), is.finite(d_int))

## --- 2. equal-weight composite (normalized within design x resolution x domain) ----------
grp <- c('design', 'resolution', 'domain')
per_ratio_d <- reps |> group_by(across(all_of(c(grp, 'ratio')))) |>
  summarise(across(c(d_nAP, d_roc, d_slope, d_int), ~ mean(.x, na.rm = TRUE)), .groups = 'drop')
scaling <- per_ratio_d |> group_by(across(all_of(grp))) |>
  summarise(across(c(d_nAP, d_roc, d_slope, d_int),
                   list(mn = ~ min(.x, na.rm = TRUE), mx = ~ max(.x, na.rm = TRUE))), .groups = 'drop')
norm01 <- function(x, mn, mx) ifelse(mx > mn, (x - mn) / (mx - mn), 0)
w <- weights / sum(weights)

reps <- reps |> left_join(scaling, by = grp) |>
  mutate(composite = w['nAP']   * norm01(d_nAP,   d_nAP_mn,   d_nAP_mx) +
                     w['roc']   * norm01(d_roc,   d_roc_mn,   d_roc_mx) +
                     w['slope'] * norm01(d_slope, d_slope_mn, d_slope_mx) +
                     w['intercept'] * norm01(d_int, d_int_mn, d_int_mx))

per_ratio <- reps |> group_by(across(all_of(c(grp, 'ratio')))) |>
  summarise(n = n(), composite_mean = mean(composite, na.rm = TRUE),
            nAP_mean = mean(nAP, na.rm = TRUE), roc_mean = mean(roc_auc, na.rm = TRUE),
            cox_slope_med = median(cox_slope, na.rm = TRUE),
            cox_int_med = median(cox_intercept, na.rm = TRUE), .groups = 'drop')

## --- 3. selection: 5% rule (smallest ratio within rel_tol of the best mean composite) ----
select_ratio <- function(d){
  best <- d[which.min(d$composite_mean), ]
  within_tol <- filter(d, composite_mean <= best$composite_mean * (1 + rel_tol))
  tibble(best_ratio = best$ratio, best_composite = best$composite_mean,
         chosen_ratio = min(within_tol$ratio))
}
selection <- per_ratio |> group_by(across(all_of(grp))) |>
  group_modify(~ select_ratio(.x)) |> ungroup()

## median-performing replicate at the chosen ratio
median_model <- selection |> filter(!is.na(chosen_ratio)) |>
  pmap_dfr(function(design, resolution, domain, chosen_ratio, ...){
    r <- filter(reps, design == !!design, resolution == !!resolution,
                domain == !!domain, ratio == !!chosen_ratio)
    r[which.min(abs(r$composite - median(r$composite, na.rm = TRUE))), ] |>
      transmute(design, resolution, domain, chosen_ratio = ratio, median_replicate = replicate,
                composite, nAP, roc_auc, cox_slope, cox_intercept)
  })

## --- 4. write + report (and compare to the holdout-chosen ratio) -------------------------
write.csv(selection,    file.path(out_dir, 'pa-ratio-selection-GT.csv'),     row.names = FALSE)
write.csv(median_model, file.path(out_dir, 'pa-ratio-median-models-GT.csv'), row.names = FALSE)

holdout <- read.csv(file.path(out_dir, 'pa-ratio-selection.csv'))[, c('design','resolution','chosen_ratio')]
names(holdout)[3] <- 'holdout_ratio'
cmp <- selection |>
  select(design, resolution, domain, GT_ratio = chosen_ratio) |>
  left_join(holdout, by = c('design','resolution')) |>
  left_join(select(per_ratio, design, resolution, domain, ratio, nAP_mean, roc_mean,
                   cox_slope_med, cox_int_med),
            by = c('design','resolution','domain','GT_ratio' = 'ratio')) |>
  arrange(domain, design, resolution)

options(width = 200)
message(sprintf('\n=== PA ratio chosen on FIELD GT (%s, 5%% rule), per domain, vs holdout ===', CALIB))
cmp |> as.data.frame() |> print(digits = 3, row.names = FALSE)
message('\n=== GT median models at the GT-chosen ratio ===')
print(as.data.frame(median_model), digits = 3, row.names = FALSE)
