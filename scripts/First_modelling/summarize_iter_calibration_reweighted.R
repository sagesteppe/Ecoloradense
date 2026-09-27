## Holdout Brier / Cox calibration / decision-curve scores for every finished
## iteration-0 and iteration-1 SDM production fit, reweighted to a common set of
## target prevalences - the companion to summarize_iter_calibration_dca.R.
##
## Why: each run's holdout inherits the presence:absence ratio its model was built
## with (PAratio), so an unweighted holdout grades every model against the very
## prevalence it was trained to predict - a 1:1 model is scored on a ~50% presence
## holdout, a 1:4 model on a ~20% one - and both look calibrated. Brier, Cox and DCA
## are prevalence-dependent, so those scores can't separate PAratios. Reweighting
## each holdout (functions.R's prevalenceWeights(): presences x target/observed,
## absences x (1 - target)/(1 - observed)) scores every model against the SAME
## prevalence, answering "if presences were `target` of the landscape, whose
## probabilities come out right?" without dropping any holdout points. ROC-AUC is
## unchanged by the weights (reported as a check); PR-AUC, Brier, Cox and DCA move.
##
## The Cox standard errors under weighting are sandwich (robust) errors; seed-to-seed
## spread (the sd column of the summaries) is the more honest uncertainty here.
##
## Run with the working directory set to scripts/ - e.g.
## `Rscript First_modelling/summarize_iter_calibration_reweighted.R`.

suppressPackageStartupMessages({
  library(dplyr)
  library(ranger)
})
source('functions.R') # presenceScores(), prevalenceWeights()

# Target prevalences to score every holdout against, 1% to 40% by 1%. The low end is
# the realistic one - at the pixel scale most rare species occupy ~1% of a landscape;
# 0.34 is the 2026 field ground truth's observed prevalence (its sampling was partly
# targeted near known populations, so it overstates landscape prevalence).
target_prevalences <- round(seq(0.01, 0.40, by = 0.01), 2)

results_root <- file.path('..', 'results')
results_dir <- file.path(results_root, 'tables')

# fname format (see functions.R's modeller()): {resolution}-Iteration{0,1}-PA{PAratio}DO:{distOrder}-Seed{seed}.csv
fname_re <- '^(1-3arc|1arc|3arc)-Iteration([01])-PA(1:[0-9.]+)DO:([0-9.]+)-Seed([0-9]+)\\.csv$'

files <- list.files(results_dir, pattern = fname_re, full.names = TRUE)
message(sprintf('Found %d finished iteration-0/1 production eval files; scoring each at %d target prevalences.',
                length(files), length(target_prevalences)))

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

  do.call(rbind, lapply(target_prevalences, function(tp){
    scores <- presenceScores(truth = truth, prob = prob, target_prevalence = tp)
    cbind(scores$metrics, target_prevalence = tp, parse_fname(f), n_test = length(truth))
  }))
}

res <- lapply(seq_along(files), function(i){
  if(i %% 25 == 0) message(sprintf('  %d / %d', i, length(files)))
  eval_one(files[i])
})
all_metrics <- do.call(rbind, Filter(Negate(is.null), res))

# sanity check: ROC-AUC is invariant to class reweighting, so it must not move with target
roc_spread <- all_metrics |>
  dplyr::filter(metric == 'roc_auc') |>
  dplyr::group_by(resolution, iteration, PAratio, distOrder, seed) |>
  dplyr::summarize(spread = max(estimate) - min(estimate), .groups = 'drop')
message(sprintf('ROC-AUC spread across target prevalences: max %.2g (should be ~0).', max(roc_spread$spread)))

summarise_metrics <- function(d, ...){
  d |>
    dplyr::filter(metric != 'roc_auc', metric != 'prevalence') |>
    dplyr::group_by(..., metric) |>
    dplyr::summarize(mean = mean(estimate), median = median(estimate), sd = sd(estimate),
                     min = min(estimate), max = max(estimate), n = dplyr::n(), .groups = 'drop')
}

summary_by_res_iter_pa <- summarise_metrics(all_metrics, target_prevalence, resolution, iteration, PAratio) |>
  dplyr::arrange(target_prevalence, iteration, resolution, PAratio, metric)

out_per_run        <- file.path(results_dir, 'iter0_iter1_calibration_reweighted_per_run.csv')
out_by_res_iter_pa <- file.path(results_dir, 'iter0_iter1_calibration_reweighted_summary_by_resolution_PAratio.csv')
write.csv(all_metrics, out_per_run, row.names = FALSE)
write.csv(summary_by_res_iter_pa, out_by_res_iter_pa, row.names = FALSE)

message('\n=== mean scaled Brier by PAratio (rows) x selected target prevalences (cols), iteration x resolution ===')
summary_by_res_iter_pa |>
  dplyr::filter(metric == 'brier_scaled', target_prevalence %in% c(0.01, 0.05, 0.10, 0.20, 0.34, 0.40)) |>
  dplyr::select(target_prevalence, iteration, resolution, PAratio, mean) |>
  tidyr::pivot_wider(names_from = target_prevalence, values_from = mean, names_prefix = 'pi=') |>
  as.data.frame() |>
  print(digits = 2)
message('\nWritten to: ', out_per_run)
message('Written to: ', out_by_res_iter_pa)
