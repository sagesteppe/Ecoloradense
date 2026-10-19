## Backfill any missing per-fit metrics across every results_* tree.
##
## Metrics checked and (re)computed from the saved model + test data:
##   brier_class    direct mean squared error
##   boyce          ecospat::ecospat.boyce (test absences as background)
##   cox_intercept, cox_slope  via coxCalibration() in functions.R
##
## For each tables/{fname}.csv, the script finds which of those metrics are missing,
## reloads models/{fname}.rds and test_data/{fname}.csv, predicts probabilities ONCE,
## computes only the missing metrics, and appends those rows to the CSV. Fits whose
## CSVs are already complete are skipped, so the script is safe to rerun.
##
## Run with wd = scripts/:
##   Rscript First_modelling/backfill_metrics.R                 # all result trees
##   Rscript First_modelling/backfill_metrics.R results_LMPO    # specific trees
##
## Supersedes backfill_boyce.R.

suppressPackageStartupMessages({
  library(ranger)
  library(ecospat)
})
source('functions.R')

args <- commandArgs(trailingOnly = TRUE)
if(length(args) == 0){
  roots <- list.dirs(PROJ_ROOT, recursive = FALSE)
  roots <- roots[grepl('^results(_|$)', basename(roots))]
} else {
  roots <- normalizePath(args, mustWork = TRUE)
}
if(length(roots) == 0) stop('no results_* trees found under ', PROJ_ROOT)

KNOWN <- c('brier_class', 'boyce', 'cox_intercept', 'cox_slope')

predict_probs <- function(rf_model, Test){
  # type = 'response' + probability = TRUE gives the P(class=1) column at a fraction
  # of the cost of modeller()'s se.method='infjack' path (which is only needed for
  # thresholding at fit time).
  pr <- predict(rf_model, Test, type = 'response', probability = TRUE,
                num.threads = max(1, parallel::detectCores() - 1))$predictions
  if(is.matrix(pr) && ncol(pr) >= 2) pr[, 2] else as.numeric(pr)
}

compute_all <- function(probs, truth){
  pres <- probs[as.character(truth) == '1']
  abs  <- probs[as.character(truth) == '0']

  brier <- data.frame(metric = 'brier_class', estimator = 'binary',
    estimate = mean((as.numeric(as.character(truth)) - probs)^2),
    std_error = NA_real_, stringsAsFactors = FALSE)

  boyce_est <- NA_real_
  if(length(pres) >= 2 && length(abs) >= 5){
    r <- tryCatch(ecospat::ecospat.boyce(fit = abs, obs = pres,
                                          PEplot = FALSE, rm.duplicate = TRUE)$cor,
                  error = function(e) NA_real_)
    if(length(r) == 1 && is.numeric(r)) boyce_est <- as.numeric(r)
  }
  boyce <- data.frame(metric = 'boyce', estimator = 'continuous',
    estimate = boyce_est, std_error = NA_real_, stringsAsFactors = FALSE)

  cox <- tryCatch(coxCalibration(truth = truth, prob = probs),
                  error = function(e) NULL)
  if(is.null(cox) || !nrow(cox))
    cox <- data.frame(metric = c('cox_intercept', 'cox_slope'),
                      estimator = 'glm',
                      estimate = NA_real_, std_error = NA_real_,
                      stringsAsFactors = FALSE)
  # keep only the two rows this script owns, in the expected order/columns
  cox <- cox[cox$metric %in% c('cox_intercept', 'cox_slope'), ]
  cox <- cox[, intersect(c('metric', 'estimator', 'estimate', 'std_error'), names(cox))]

  rbind(brier, boyce, cox)
}

tot <- c(scanned = 0L, no_missing = 0L, updated = 0L, failed = 0L, no_model = 0L,
         added_rows = 0L)
for(root in roots){
  tables_dir <- file.path(root, 'tables')
  if(!dir.exists(tables_dir)) next
  csvs <- list.files(tables_dir, pattern = '-Seed[0-9]+\\.csv$', full.names = TRUE)
  if(length(csvs) == 0) next
  message('\n=== ', basename(root), ' (', length(csvs), ' fits) ===')
  for(csv in csvs){
    tot['scanned'] <- tot['scanned'] + 1L
    fname <- sub('\\.csv$', '', basename(csv))

    met <- tryCatch(read.csv(csv, stringsAsFactors = FALSE), error = function(e) NULL)
    if(is.null(met) || !'metric' %in% names(met)){
      tot['failed'] <- tot['failed'] + 1L; next
    }
    missing <- setdiff(KNOWN, unique(met$metric))
    if(length(missing) == 0){
      tot['no_missing'] <- tot['no_missing'] + 1L; next
    }

    model_path <- file.path(root, 'models',    paste0(fname, '.rds'))
    test_path  <- file.path(root, 'test_data', paste0(fname, '.csv'))
    if(!file.exists(model_path) || !file.exists(test_path)){
      tot['no_model'] <- tot['no_model'] + 1L; next
    }

    new_rows <- tryCatch({
      rf <- readRDS(model_path)
      te <- read.csv(test_path, stringsAsFactors = FALSE)
      if(names(te)[1] %in% c('X', '')) te <- te[, -1, drop = FALSE]
      te$Occurrence <- as.factor(te$Occurrence)
      p <- predict_probs(rf, te)
      compute_all(p, te$Occurrence)
    }, error = function(e){
      message('  ', fname, ': ', conditionMessage(e)); NULL
    })
    if(is.null(new_rows)){
      tot['failed'] <- tot['failed'] + 1L; next
    }

    new_rows <- new_rows[new_rows$metric %in% missing, , drop = FALSE]
    if(!nrow(new_rows)){
      tot['failed'] <- tot['failed'] + 1L; next
    }

    # align the appended rows to the existing CSV's column order (res/iter from `met`)
    new_rows$resolution <- unique(met$resolution)[1]
    new_rows$iteration  <- unique(met$iteration)[1]
    for(c in setdiff(names(met), names(new_rows))) new_rows[[c]] <- NA
    new_rows <- new_rows[, names(met), drop = FALSE]

    write.csv(rbind(met, new_rows), csv, row.names = FALSE)
    tot['updated']    <- tot['updated']    + 1L
    tot['added_rows'] <- tot['added_rows'] + nrow(new_rows)
    cat('.')
  }
  cat('\n')
}

message('\n--- backfill_metrics.R summary ---')
print(tot)
message(sprintf('updated %d of %d fits (%d rows added; %d already complete; %d failed; %d missing model/test).',
                tot['updated'], tot['scanned'], tot['added_rows'],
                tot['no_missing'], tot['failed'], tot['no_model']))
