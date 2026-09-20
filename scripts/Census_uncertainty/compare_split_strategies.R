## Characterizes how much of Phase 1's density-model comparison (which brms
## family/ML candidate "wins") is genuine signal vs. split-to-split noise, and
## whether the twinning split splitData() uses by default actually earns its
## complexity over simpler alternatives - across 4 train/test split
## strategies (functions.R's split_indices(): 'twinning', 'classic',
## 'spatial_knn', 'population').
##
## Cost: one densityModeller() call per replicate is the full pipeline (RFE +
## 5 ML candidates + the 4-family x 2-CV-mode brms sweep), ~30min measured on
## the 1-3arc/PA1:1 product. 'twinning'/'classic' are stochastic (seeded) and
## use a racing-style stopping rule: run replicates 1-5, then only continue to
## 6-10 if the leading and runner-up Models' held-out MAE aren't yet clearly
## separated (see race_check()). 'spatial_knn' has exactly k=10 natural
## replicates (which CAST::knndm() fold is held out) and 'population' is the
## exhaustive leave-one-population-out set (as many replicates as
## populations) - neither has a stopping rule, both always run to completion.
##
## Run with the working directory set to scripts/ (this project's usual
## convention) - e.g. `Rscript Census_uncertainty/compare_split_strategies.R <method>`.

suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(terra); library(brms)
  library(caret); library(future); library(CAST); library(gstat); library(twinning)
  library(recipes); library(rsample); library(parsnip); library(tune)
  library(finetune); library(dials); library(yardstick); library(workflows)
  library(bonsai); library(lightgbm); library(xgboost)
})

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')
source('functions.R')
source('Census_uncertainty/density_bayes.R')
source('Census_uncertainty/density_conformal.R')

x <- '1-3arc-Iteration1-PA1:1DO:0-Seed1021-Pr.tif'
bn_base <- gsub('DO.*$', '', x)
fp <- file.path('..', 'results', 'count_models')

df <- wrapper(x, return_early = TRUE)
coords <- sf::st_coordinates(df)
df$Longitude <- coords[, 1]
df$Latitude  <- coords[, 2]

#' Compare the leading and runner-up Models' held-out MAE across completed replicates.
#'
#' @description A racing-style stop rule in the spirit of the ML candidates'
#' own `finetune::tune_race_anova()` (../functions.R's `poiss()`/`tweed()`/
#' `gbs()`), simplified to a transparent gap-vs-spread heuristic rather than a
#' formal ANOVA, appropriate at n=5-10 replicates: the leader/runner-up gap is
#' "resolved" once it exceeds the sum of their own standard deviations - i.e.
#' the two candidates' MAE distributions across replicates are not
#' obviously overlapping.
#' @param combined data.frame of stacked `densityModeller()$metrics`, tagged
#' with a `replicate` column.
#' @return list(table = data.frame(Model, mean_mae, sd_mae, n), resolved =,
#' leader =, runner_up =).
race_check <- function(combined){
  mae <- combined |>
    dplyr::filter(Metric == 'MAE') |>
    dplyr::group_by(Model) |>
    dplyr::summarize(mean_mae = mean(Value), sd_mae = stats::sd(Value), n = dplyr::n(), .groups = 'drop') |>
    dplyr::arrange(mean_mae)

  top2 <- mae[1:2, ]
  resolved <- (top2$mean_mae[2] - top2$mean_mae[1]) > (top2$sd_mae[1] + top2$sd_mae[2])

  list(table = mae, resolved = isTRUE(resolved), leader = top2$Model[1], runner_up = top2$Model[2])
}

#' Run one split-strategy replicate through the full densityModeller() pipeline.
#'
#' @param method,replicate,k as `split_indices()`.
#' @return `densityModeller()$metrics`, tagged with `replicate` and `n_test`
#' columns - `n_test` varies a lot across replicates (LOPO ranges from n=4 to
#' n=78 depending which population is held out; `spatial_knn`'s folds aren't
#' equal-sized either), so a replicate's MAE/MSE/RMSE isn't equally precise
#' across the board - carrying `n_test` through lets that be accounted for
#' (e.g. weighting, or just flagging) rather than silently averaged as if
#' every replicate were equally informative.
run_replicate <- function(method, replicate, k = 10, tune_metric = 'mae'){
  bn_i <- paste0(bn_base, '-', method, 'Rep', replicate)
  message(sprintf('\n===== %s replicate %s (bn = %s, tune_metric = %s) =====', method, replicate, bn_i, tune_metric))
  t0 <- Sys.time()
  res <- densityModeller(
    df, bn = bn_i, fp = fp,
    brms_refit_args = list(full_chains = 2, full_iter = 1000, full_warmup = 500),
    split_method = method, split_replicate = replicate, split_k = k, tune_metric = tune_metric
  )
  message(sprintf('[timing] replicate %s total: %.1fs', replicate, as.numeric(difftime(Sys.time(), t0, units = 'secs'))))
  res$metrics$replicate <- replicate
  res$metrics$n_test <- length(split_indices(df, method = method, replicate = replicate, k = k))
  res$metrics
}

#' Run a full split-strategy comparison (racing for stochastic methods,
#' exhaustive for 'spatial_knn'/'population'), incrementally saving progress.
#'
#' @param method as `split_indices()`.
#' @param out_csv path to write (and resume from) accumulated results.
#' @param tune_metric as `resolve_tune_metric()` - threaded through to every
#' replicate's `densityModeller()` call.
run_strategy <- function(method, out_csv, tune_metric = 'mae'){

  all_results <- if(file.exists(out_csv)){
    list(readr::read_csv(out_csv, show_col_types = FALSE))
  } else list()
  done <- if(length(all_results)) unique(all_results[[1]]$replicate) else integer(0)

  save_progress <- function(){
    combined <- dplyr::bind_rows(all_results)
    write.csv(combined, out_csv, row.names = FALSE)
    combined
  }

  replicate_ids <- switch(method,
    twinning = , classic = 1:5,
    spatial_knn = 1:10,
    # true-zero populations (e.g. CP/WMB - surveyed, never occupied) are
    # excluded as LOPO *test* targets: holding one out and predicting "0
    # plants" everywhere is trivially achievable and doesn't test a count
    # model's actual generalization the way holding out a real (nonzero)
    # population does. They still appear normally in every OTHER
    # replicate's training data (as real, informative zero rows) - only
    # excluded from being the held-out target itself.
    population = df |>
      sf::st_drop_geometry() |>
      dplyr::group_by(Lctn_bb) |>
      dplyr::summarize(has_presence = any(Prsnc_All > 0), .groups = 'drop') |>
      dplyr::filter(has_presence) |>
      dplyr::pull(Lctn_bb) |>
      unique(),
    stop("run_strategy(): unknown method '", method, "'"))

  for(r in replicate_ids){
    if(r %in% done) next
    all_results[[length(all_results) + 1]] <- run_replicate(method, r, tune_metric = tune_metric)
    save_progress()
  }

  combined <- save_progress()

  # racing extension: only for the stochastic methods, only if not yet resolved.
  if(method %in% c('twinning', 'classic')){
    rc <- race_check(combined)
    message('\n=== after replicates 1-5 (', method, ') ===')
    print(as.data.frame(rc$table))
    message('resolved: ', rc$resolved, ' (leader = ', rc$leader, ', runner_up = ', rc$runner_up, ')')

    if(!rc$resolved){
      message('\nGap not resolved after 5 replicates - continuing to 6-10...')
      for(r in 6:10){
        if(r %in% done) next
        all_results[[length(all_results) + 1]] <- run_replicate(method, r, tune_metric = tune_metric)
        save_progress()
      }
      combined <- save_progress()
    }
  }

  combined
}

# ---- entry point: method and tune_metric taken from the command line ----
# Rscript compare_split_strategies.R <method> [tune_metric]
# tune_metric: 'mae' (default), 'huber', or 'poisson' (resolve_tune_metric()).
args <- commandArgs(trailingOnly = TRUE)
method <- if(length(args) >= 1) args[[1]] else 'twinning'
tune_metric <- if(length(args) >= 2) args[[2]] else 'mae'
metric_suffix <- if(tune_metric == 'mae') '' else paste0('-', tune_metric)
out_csv <- file.path('..', 'results', 'tables', paste0(bn_base, metric_suffix, '-splitcompare-', method, '.csv'))

final <- run_strategy(method, out_csv, tune_metric = tune_metric)

message('\n=== FINAL (', method, ', tune_metric = ', tune_metric, ') ===')
print(as.data.frame(race_check(final)$table))
message('\nWritten to: ', out_csv)
