## Consolidate every per-fit 2026 FIELD-GROUNDTRUTH metric CSV under [Rr]esults_* into one
## wide table. These '*-2026gt.csv' files score each fit against the common external 2026
## field dataset (n = 283 presence/absence points), as opposed to the fit's own train/test
## holdout (which aggregate_all_metrics.R consolidates). The '-2026gt-270m.csv' variant,
## restricted to points within 270 m of known occurrences, is deliberately NOT included.
##
## Unlike the own-holdout files, these carry an actual `prevalence` column (the field
## prevalence of the 283-point set), so nAP is computed against that real prevalence rather
## than 1/(1+ratio). Boyce is absent (it needs presence/background, not the field P/A set).
##
## One row per fit, columns:
##   tree, resolution, iteration, PAratio, distOrder, seed, fname, n,
##   pr_auc, roc_auc, brier_class, brier_scaled, prevalence,
##   cox_intercept, cox_intercept_se, cox_slope, cox_slope_se
##
## Two variants, selected by the sole argument:
##   (none)  full 283-point field set  -> *-2026gt.csv       -> all_fit_2026gt.csv
##   270m    near-population subset     -> *-2026gt-270m.csv  -> all_fit_2026gt_270m.csv
##          (points within 270 m of a known occurrence; smaller n, higher prevalence)
##
## Run with wd = scripts/:
##   Rscript First_modelling/aggregate_2026gt_metrics.R
##   Rscript First_modelling/aggregate_2026gt_metrics.R 270m

suppressPackageStartupMessages(library(dplyr))
source('functions.R')

args     <- commandArgs(trailingOnly = TRUE)
suffix   <- if (length(args) && args[1] == '270m') '-2026gt-270m' else '-2026gt'
out_stem <- paste0('all_fit', gsub('-', '_', suffix))   # all_fit_2026gt | all_fit_2026gt_270m

roots <- list.dirs(PROJ_ROOT, recursive = FALSE)
roots <- roots[grepl('^[Rr]esults(_|$)', basename(roots))]

pull_metric <- function(met, name, col = 'estimate'){
  v <- met[[col]][met$metric == name]
  if(length(v) == 0) NA else v[1]
}

rows <- list()
for(root in roots){
  tables_dir <- file.path(root, 'tables')
  if(!dir.exists(tables_dir)) next
  csvs <- list.files(tables_dir, pattern = paste0('-Seed[0-9]+', suffix, '\\.csv$'), full.names = TRUE)
  for(csv in csvs){
    fname <- sub(paste0(suffix, '$'), '', sub('\\.csv$', '', basename(csv)))
    pf <- fname_parse(fname); if(is.na(pf$resolution)) next
    met <- tryCatch(read.csv(csv, stringsAsFactors = FALSE), error = function(e) NULL)
    if(is.null(met) || !'metric' %in% names(met)) next
    rows[[length(rows) + 1]] <- data.frame(
      tree             = basename(root),
      resolution       = pf$resolution,
      iteration        = pf$iteration,
      PAratio          = pf$PAratio,
      distOrder        = pf$distOrder,
      seed             = pf$seed,
      fname            = fname,
      n                = pull_metric(met, 'pr_auc', 'n'),
      pr_auc           = pull_metric(met, 'pr_auc'),
      roc_auc          = pull_metric(met, 'roc_auc'),
      brier_class      = pull_metric(met, 'brier_class'),
      brier_scaled     = pull_metric(met, 'brier_scaled'),
      prevalence       = pull_metric(met, 'prevalence'),
      cox_intercept    = pull_metric(met, 'cox_intercept'),
      cox_intercept_se = pull_metric(met, 'cox_intercept', 'std_error'),
      cox_slope        = pull_metric(met, 'cox_slope'),
      cox_slope_se     = pull_metric(met, 'cox_slope', 'std_error'),
      stringsAsFactors = FALSE)
  }
}
out <- do.call(rbind, rows)

out_dir <- file.path(PROJ_ROOT, 'results_summary')
dir.create(out_dir, showWarnings = FALSE)
out_path <- file.path(out_dir, paste0(out_stem, '.csv'))
write.csv(out, out_path, row.names = FALSE)

message(sprintf('wrote %d rows (fits) across %d trees to %s',
                nrow(out), dplyr::n_distinct(out$tree), out_path))
print(dplyr::count(out, tree, resolution, name = 'fits') |> as.data.frame())
