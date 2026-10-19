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
roots <- roots[grepl('^results(_|$)', basename(roots))]

# {res}-Iteration{it}-PA1:{ratio}DO:{distOrder}-Seed{seed}
# PAratio is stored as "1:X" in the filename (X may be fractional, e.g. 1.079).
re_fname <- '^(.+)-Iteration([0-9]+)-PA1:([0-9.]+)DO:([0-9]+)-Seed([0-9]+)$'

parse_fname <- function(f){
  m <- regmatches(f, regexec(re_fname, f))[[1]]
  if(length(m) != 6) return(NULL)
  list(resolution = m[2], iteration = as.integer(m[3]),
       PAratio = as.numeric(m[4]), distOrder = as.integer(m[5]),
       seed = as.integer(m[6]))
}

pull_metric <- function(met, name, col = 'estimate'){
  v <- met[[col]][met$metric == name]
  if(length(v) == 0) NA else v[1]
}

rows <- list()
for(root in roots){
  tables_dir <- file.path(root, 'tables')
  if(!dir.exists(tables_dir)) next
  csvs <- list.files(tables_dir, pattern = '-Seed[0-9]+\\.csv$', full.names = TRUE)
  for(csv in csvs){
    fname <- sub('\\.csv$', '', basename(csv))
    pf <- parse_fname(fname); if(is.null(pf)) next
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
