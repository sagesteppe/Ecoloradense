## Consolidate every per-fit metric CSV under results_* into one wide table.
##
## One row per fit, columns:
##   tree, resolution, iteration, PAratio, distOrder, seed, fname,
##   pr_auc, roc_auc, brier_class, boyce,
##   cox_intercept, cox_intercept_se, cox_slope, cox_slope_se
##
## Run with wd = scripts/:
##   Rscript First_modelling/aggregate_all_metrics.R
## Writes results_summary/all_fit_metrics.csv.

suppressPackageStartupMessages(library(dplyr))
source('functions.R')

roots <- list.dirs(PROJ_ROOT, recursive = FALSE)
roots <- roots[grepl('^[Rr]esults(_|$)', basename(roots))]

## Pull one metric's value. `prefer` selects a `calibration` variant when the table carries
## that column (backfill_metrics.R writes both 'raw' and, where a -platt.rds exists,
## 'calibrated' rows for the scale-dependent metrics). When the preferred variant is absent
## for this fit (e.g. legacy tables with no calibration column, or no Platt map), it falls
## back to the first matching row (raw).
pull_metric <- function(met, name, col = 'estimate', prefer = NULL){
  sel <- met$metric == name
  if(!is.null(prefer) && 'calibration' %in% names(met)){
    cal <- sel & met$calibration == prefer
    if(any(cal)) sel <- cal
  }
  v <- met[[col]][sel]
  if(length(v) == 0) NA else v[1]
}

rows <- list()
for(root in roots){
  tables_dir <- file.path(root, 'tables')
  if(!dir.exists(tables_dir)) next
  csvs <- list.files(tables_dir, pattern = '-Seed[0-9]+\\.csv$', full.names = TRUE)
  for(csv in csvs){
    fname <- sub('\\.csv$', '', basename(csv))
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
      pr_auc           = pull_metric(met, 'pr_auc'),
      roc_auc          = pull_metric(met, 'roc_auc'),
      brier_class      = pull_metric(met, 'brier_class'),
      boyce            = pull_metric(met, 'boyce'),
      ## cox calibration: RAW on the holdout. The Platt map is fit on the training inner-CV
      ## OOF, so it transfers badly to extrapolative holdouts (spatial/LMPO) and over-sharpens
      ## the slope there - see memory 'platt-overcorrects-extrapolative-holdout'. Platt's real
      ## benefit is a within-AOA deployment property, shown on the 2026 GT path, not here.
      ## (Calibrated rows still live in the per-fit tables; pass prefer='calibrated' to use them.)
      cox_intercept    = pull_metric(met, 'cox_intercept'),
      cox_intercept_se = pull_metric(met, 'cox_intercept', 'std_error'),
      cox_slope        = pull_metric(met, 'cox_slope'),
      cox_slope_se     = pull_metric(met, 'cox_slope', 'std_error'),
      stringsAsFactors = FALSE)
  }
}
out <- do.call(rbind, rows)

# Supplement with 2026-ground-truth point-prediction scores (predict_gridsearch_to_points.R
# scores every model at the 283 2026 field points - an EXTERNAL holdout, independent of the
# fit-time train/test split that drove the pr_auc/roc_auc/brier/boyce/cox columns above).
# Only rows for fits whose model is in the pointpredict table will carry gt2026_* values;
# the rest stay NA.
gt_path <- file.path(PROJ_ROOT, 'results', 'tables',
                     '2026-groundtruth-gridsearch-pointpredict-all.csv')
if(file.exists(gt_path)){
  gt <- read.csv(gt_path, stringsAsFactors = FALSE)
  keep <- c('tree', 'model', 'n', 'pr_auc', 'roc_auc', 'brier_class', 'brier_scaled',
            'prevalence', 'cox_intercept', 'cox_intercept_se',
            'cox_slope', 'cox_slope_se',
            'dca_prop_useful', 'dca_mean_excess_nb',
            'dca_nb_pt25', 'dca_nb_pt50', 'dca_nb_pt75')
  gt <- gt[, intersect(keep, names(gt)), drop = FALSE]
  # fname is unique within a tree but can collide across trees - join on both.
  idcols <- c('tree', 'model')
  names(gt)[!names(gt) %in% idcols] <- paste0('gt2026_', names(gt)[!names(gt) %in% idcols])
  names(gt)[names(gt) == 'model'] <- 'fname'
  out <- merge(out, gt, by = c('tree', 'fname'), all.x = TRUE, sort = FALSE)
  # bring fname-identifier columns back to the front
  front <- c('tree', 'resolution', 'iteration', 'PAratio', 'distOrder', 'seed', 'fname')
  out <- out[, c(front, setdiff(names(out), front))]
  message(sprintf('merged %d gt2026 rows (external 2026 holdout) by fname.', nrow(gt)))
} else {
  message('no 2026-groundtruth-gridsearch-pointpredict-all.csv found; skipping join.')
}

out_dir <- file.path(PROJ_ROOT, 'results_summary')
dir.create(out_dir, showWarnings = FALSE)
out_path <- file.path(out_dir, 'all_fit_metrics.csv')
write.csv(out, out_path, row.names = FALSE)

message(sprintf('wrote %d rows (fits) across %d trees to %s',
                nrow(out), dplyr::n_distinct(out$tree), out_path))
print(dplyr::count(out, tree, resolution, name = 'fits') |> as.data.frame())
