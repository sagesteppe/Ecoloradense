## One-off summary of Brier score, Cox calibration (intercept/slope) and decision
## curve analysis across every finished iteration-0 and iteration-1 SDM production
## fit - the calibration/clinical-utility companion to summarize_iter_eval_metrics.R's
## PR-AUC / ROC-AUC tables. Every run in results/tables/ predates modeller() writing
## brier_class/cox_* rows and the evaluations/*-dca.csv curves, so rather than read
## those rows this re-predicts each saved ranger model (results/models/) onto its
## saved holdout set (results/test_data/) and recomputes everything with the same
## functions.R helpers modeller() now uses. pr_auc/roc_auc are recomputed here too
## as a check that the re-predictions match what modeller() scored (they agree to
## within ~0.002 PR-AUC).
##
## Brier score depends on holdout prevalence, which the PAratio sets artificially,
## so a scaled Brier (1 - brier / brier of predicting the holdout prevalence for
## every point) is reported alongside the raw score. Decision curves are likewise
## prevalence-dependent; each curve is reduced to a few scalars for the tables
## (see dca_scalars()), the full per-run curves go to results/evaluations/
## {fname}-dca.csv (the file modeller() now writes), and the mean curve per
## resolution x iteration x PAratio goes to its own table.
##
## Run with the working directory set to scripts/ (this project's usual
## convention) - e.g. `Rscript First_modelling/summarize_iter_calibration_dca.R`.

suppressPackageStartupMessages({
  library(dplyr)
  library(ranger)
})
source('functions.R') # presenceScores(), dca_scalars()

results_root <- file.path('..', 'results')
results_dir <- file.path(results_root, 'tables')

# fname format (see functions.R's modeller()): {resolution}-Iteration{0,1}-PA{PAratio}DO:{distOrder}-Seed{seed}.csv
fname_re <- '^(1-3arc|1arc|3arc)-Iteration([01])-PA(1:[0-9.]+)DO:([0-9.]+)-Seed([0-9]+)\\.csv$'

files <- list.files(results_dir, pattern = fname_re, full.names = TRUE)
message(sprintf('Found %d finished iteration-0/1 production eval files.', length(files)))

parse_fname <- function(f){
  bn <- basename(f)
  m <- regmatches(bn, regexec(fname_re, bn))[[1]]
  data.frame(
    resolution = m[2],
    iteration  = m[3],
    PAratio    = as.numeric(sub('^1:', '', m[4])),
    distOrder  = as.numeric(m[5]),
    seed       = as.numeric(m[6])
  )
}

eval_one <- function(f){
  fname <- sub('\\.csv$', '', basename(f))
  model_path <- file.path(results_root, 'models', paste0(fname, '.rds'))
  test_path  <- file.path(results_root, 'test_data', paste0(fname, '.csv'))
  if(!file.exists(model_path) || !file.exists(test_path)){
    message('Skipping ', fname, ': missing saved model or test data.')
    return(NULL)
  }

  rf_model <- readRDS(model_path)
  Test <- read.csv(test_path)
  prob <- predict(rf_model, Test)$predictions[,2]
  truth <- Test$Occurrence

  scores <- presenceScores(truth = truth, prob = prob)
  dca_path <- file.path(results_root, 'evaluations', paste0(fname, '-dca.csv'))
  if(!file.exists(dca_path)) write.csv(scores$dca, dca_path, row.names = FALSE)

  list(
    metrics = cbind(scores$metrics, parse_fname(f), n_test = length(truth)),
    dca = cbind(scores$dca, parse_fname(f))
  )
}

res <- lapply(seq_along(files), function(i){
  if(i %% 25 == 0) message(sprintf('  %d / %d', i, length(files)))
  eval_one(files[i])
})
res <- Filter(Negate(is.null), res)

all_metrics <- do.call(rbind, lapply(res, `[[`, 'metrics'))
all_dca     <- do.call(rbind, lapply(res, `[[`, 'dca'))

# sanity check: recomputed pr_auc should match what modeller() stored
stored <- do.call(rbind, lapply(files, function(f){
  d <- read.csv(f); d <- d[d$metric == 'pr_auc', ]
  if(nrow(d) == 0) return(NULL)
  cbind(stored = d$estimate, parse_fname(f))
}))
chk <- dplyr::inner_join(
  dplyr::filter(all_metrics, metric == 'pr_auc'), stored,
  by = c('resolution', 'iteration', 'PAratio', 'distOrder', 'seed'))
message(sprintf('pr_auc recomputed vs stored: max abs difference %.2g over %d runs.',
                max(abs(chk$estimate - chk$stored)), nrow(chk)))

summarise_metrics <- function(d, ...){
  d |>
    dplyr::filter(metric != 'pr_auc', metric != 'roc_auc') |>
    dplyr::group_by(..., metric) |>
    dplyr::summarize(mean = mean(estimate), median = median(estimate), sd = sd(estimate),
                     min = min(estimate), max = max(estimate), n = dplyr::n(), .groups = 'drop')
}

summary_by_res_iter <- summarise_metrics(all_metrics, resolution, iteration) |>
  dplyr::arrange(iteration, resolution, metric)

summary_by_res_iter_pa <- summarise_metrics(all_metrics, resolution, iteration, PAratio) |>
  dplyr::arrange(iteration, resolution, PAratio, metric)

dca_curve <- all_dca |>
  dplyr::group_by(resolution, iteration, PAratio, strategy, threshold) |>
  dplyr::summarize(mean_net_benefit = mean(net_benefit), sd_net_benefit = sd(net_benefit),
                   n = dplyr::n(), .groups = 'drop') |>
  dplyr::arrange(iteration, resolution, PAratio, strategy, threshold)

out_per_run        <- file.path(results_dir, 'iter0_iter1_calibration_dca_per_run.csv')
out_by_res_iter    <- file.path(results_dir, 'iter0_iter1_calibration_dca_summary_by_resolution.csv')
out_by_res_iter_pa <- file.path(results_dir, 'iter0_iter1_calibration_dca_summary_by_resolution_PAratio.csv')
out_dca_curve      <- file.path(results_dir, 'iter0_iter1_dca_curve_by_resolution_PAratio.csv')
write.csv(all_metrics, out_per_run, row.names = FALSE)
write.csv(summary_by_res_iter, out_by_res_iter, row.names = FALSE)
write.csv(summary_by_res_iter_pa, out_by_res_iter_pa, row.names = FALSE)
write.csv(dca_curve, out_dca_curve, row.names = FALSE)

message('\n=== Brier / Cox / DCA by resolution x iteration (across every PAratio/seed run) ===')
print(as.data.frame(summary_by_res_iter))
message('\nWritten to: ', out_per_run)
message('Written to: ', out_by_res_iter)
message('Written to: ', out_by_res_iter_pa)
message('Written to: ', out_dca_curve)
