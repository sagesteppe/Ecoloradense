## Holdout Brier / Cox calibration / decision-curve scores for every finished
## iteration-0 and iteration-1 SDM production fit.
##
## Two modes, controlled by --prevalences:
##
##   Unweighted (default): one unweighted pass per fit, using the holdout's
##   own prevalence. Also writes DCA curves to results/evaluations/ and the
##   mean curve per resolution x iteration x PAratio. Sanity-checks that the
##   recomputed pr_auc matches what modeller() stored. Output stems:
##   iter0_iter1_calibration_dca_*.
##
##   Reweighted (--prevalences=0.01,0.05,...): re-scores every fit against
##   each target prevalence via functions.R's prevalenceWeights(), answering
##   "if presences were `target` of the landscape, whose probabilities come
##   out right?". Sanity-checks that ROC-AUC is invariant to the reweighting.
##   For highPA_grid.R trees, also writes a per-split summary. Output stems:
##   iter0_iter1_calibration_reweighted_*.
##
## Run with the working directory set to scripts/ - e.g.
##   Rscript First_modelling/summarize_iter_calibration.R                        # unweighted on results/
##   Rscript First_modelling/summarize_iter_calibration.R results_highPA         # unweighted on an alt tree
##   Rscript First_modelling/summarize_iter_calibration.R results_highPA --prevalences=0.01,0.05,0.10,0.20,0.34,0.40

suppressPackageStartupMessages({
  library(dplyr)
  library(ranger)
})
source('functions.R') # presenceScores(), dca_scalars(), prevalenceWeights()

args <- commandArgs(trailingOnly = TRUE)
prev_arg <- grep('^--prevalences=', args, value = TRUE)
target_prevalences <- if(length(prev_arg)){
  as.numeric(strsplit(sub('^--prevalences=', '', prev_arg[1]), ',')[[1]])
} else NULL
reweight <- !is.null(target_prevalences)

positional <- args[!grepl('^--', args)]
results_root <- file.path('..', if(length(positional) >= 1) positional[1] else 'results')
results_dir <- file.path(results_root, 'tables')
mode_label <- if(reweight) 'reweighted' else 'dca'
out_stem <- paste0('iter0_iter1_calibration_', mode_label)

# fname format (see functions.R's modeller()): {resolution}-Iteration{0,1}-PA{PAratio}DO:{distOrder}-Seed{seed}.csv
fname_re <- '^(1-3arc|1arc|3arc)-Iteration([01])-PA(1:[0-9.]+)DO:([0-9.]+)-Seed([0-9]+)\\.csv$'

files <- list.files(results_dir, pattern = fname_re, full.names = TRUE)
message(sprintf('Found %d finished iteration-0/1 production eval files%s.',
                length(files),
                if(reweight) sprintf('; scoring each at %d target prevalences', length(target_prevalences)) else ''))

parse_fname <- function(f){
  bn <- basename(f)
  m <- regmatches(bn, regexec(fname_re, bn))[[1]]
  seed <- as.numeric(m[6])
  data.frame(
    resolution = m[2],
    iteration  = m[3],
    PAratio    = as.numeric(sub('^1:', '', m[4])),
    distOrder  = as.numeric(m[5]),
    seed       = seed,
    # highPA_grid.R's spatial splits: seeds 5000 + 100 * (split - 1) + 1..n
    split      = if(grepl('results_highPA(_fieldabs)?$', results_root) && seed > 5000) (seed - 5001) %/% 100 + 1 else NA
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

  if(reweight){
    list(metrics = do.call(rbind, lapply(target_prevalences, function(tp){
      scores <- presenceScores(truth = truth, prob = prob, target_prevalence = tp)
      cbind(scores$metrics, target_prevalence = tp, parse_fname(f), n_test = length(truth))
    })))
  } else {
    scores <- presenceScores(truth = truth, prob = prob)
    dca_path <- file.path(results_root, 'evaluations', paste0(fname, '-dca.csv'))
    if(!file.exists(dca_path)) write.csv(scores$dca, dca_path, row.names = FALSE)
    list(
      metrics = cbind(scores$metrics, parse_fname(f), n_test = length(truth)),
      dca     = cbind(scores$dca, parse_fname(f))
    )
  }
}

res <- lapply(seq_along(files), function(i){
  if(i %% 25 == 0) message(sprintf('  %d / %d', i, length(files)))
  eval_one(files[i])
})
res <- Filter(Negate(is.null), res)
all_metrics <- do.call(rbind, lapply(res, `[[`, 'metrics'))

summarise_metrics <- function(d, ..., drop_metrics){
  d |>
    dplyr::filter(!metric %in% drop_metrics) |>
    dplyr::group_by(..., metric) |>
    dplyr::summarize(mean = mean(estimate), median = median(estimate), sd = sd(estimate),
                     min = min(estimate), max = max(estimate), n = dplyr::n(), .groups = 'drop')
}

out_per_run <- file.path(results_dir, paste0(out_stem, '_per_run.csv'))
write.csv(all_metrics, out_per_run, row.names = FALSE)
message('Written to: ', out_per_run)

if(reweight){
  # ROC-AUC is invariant to class reweighting; this is the correctness check.
  roc_spread <- all_metrics |>
    dplyr::filter(metric == 'roc_auc') |>
    dplyr::group_by(resolution, iteration, PAratio, distOrder, seed) |>
    dplyr::summarize(spread = max(estimate) - min(estimate), .groups = 'drop')
  message(sprintf('ROC-AUC spread across target prevalences: max %.2g (should be ~0).', max(roc_spread$spread)))

  summary_by_res_iter_pa <- summarise_metrics(all_metrics, target_prevalence, resolution, iteration, PAratio,
                                              drop_metrics = c('roc_auc', 'prevalence')) |>
    dplyr::arrange(target_prevalence, iteration, resolution, PAratio, metric)
  out_by_res_iter_pa <- file.path(results_dir, paste0(out_stem, '_summary_by_resolution_PAratio.csv'))
  write.csv(summary_by_res_iter_pa, out_by_res_iter_pa, row.names = FALSE)
  message('Written to: ', out_by_res_iter_pa)

  if(any(!is.na(all_metrics$split))){
    out_by_split <- file.path(results_dir, paste0(out_stem, '_summary_by_split_PAratio.csv'))
    summary_by_split <- summarise_metrics(all_metrics, target_prevalence, resolution, iteration, split, PAratio,
                                           drop_metrics = c('roc_auc', 'prevalence')) |>
      dplyr::arrange(target_prevalence, iteration, resolution, split, PAratio, metric)
    write.csv(summary_by_split, out_by_split, row.names = FALSE)
    message('Written to: ', out_by_split)
  }

  message('\n=== mean scaled Brier by PAratio (rows) x selected target prevalences (cols), iteration x resolution ===')
  summary_by_res_iter_pa |>
    dplyr::filter(metric == 'brier_scaled',
                  target_prevalence %in% c(0.01, 0.05, 0.10, 0.20, 0.34, 0.40)) |>
    dplyr::select(target_prevalence, iteration, resolution, PAratio, mean) |>
    tidyr::pivot_wider(names_from = target_prevalence, values_from = mean, names_prefix = 'pi=') |>
    as.data.frame() |>
    print(digits = 2)

} else {
  # pr_auc stored-vs-recomputed sanity check
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

  summary_by_res_iter    <- summarise_metrics(all_metrics, resolution, iteration,
                                              drop_metrics = c('pr_auc', 'roc_auc')) |>
    dplyr::arrange(iteration, resolution, metric)
  summary_by_res_iter_pa <- summarise_metrics(all_metrics, resolution, iteration, PAratio,
                                              drop_metrics = c('pr_auc', 'roc_auc')) |>
    dplyr::arrange(iteration, resolution, PAratio, metric)

  out_by_res_iter    <- file.path(results_dir, paste0(out_stem, '_summary_by_resolution.csv'))
  out_by_res_iter_pa <- file.path(results_dir, paste0(out_stem, '_summary_by_resolution_PAratio.csv'))
  write.csv(summary_by_res_iter, out_by_res_iter, row.names = FALSE)
  write.csv(summary_by_res_iter_pa, out_by_res_iter_pa, row.names = FALSE)
  message('Written to: ', out_by_res_iter)
  message('Written to: ', out_by_res_iter_pa)

  all_dca <- do.call(rbind, lapply(res, `[[`, 'dca'))
  dca_curve <- all_dca |>
    dplyr::group_by(resolution, iteration, PAratio, strategy, threshold) |>
    dplyr::summarize(mean_net_benefit = mean(net_benefit), sd_net_benefit = sd(net_benefit),
                     n = dplyr::n(), .groups = 'drop') |>
    dplyr::arrange(iteration, resolution, PAratio, strategy, threshold)
  out_dca_curve <- file.path(results_dir, 'iter0_iter1_dca_curve_by_resolution_PAratio.csv')
  write.csv(dca_curve, out_dca_curve, row.names = FALSE)
  message('Written to: ', out_dca_curve)

  message('\n=== Brier / Cox / DCA by resolution x iteration (across every PAratio/seed run) ===')
  print(as.data.frame(summary_by_res_iter))
}
