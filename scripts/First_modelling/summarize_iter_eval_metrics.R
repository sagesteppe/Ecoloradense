## One-off summary of PR-AUC / ROC-AUC across every finished iteration-0 and
## iteration-1 SDM production fit - i.e. every results/tables/*.csv modeller()
## has written from adaptive_PAratio_search()/distOrder_PAratio_simulator()'s
## normal (spatial_class_split or classic) train/test split, across whatever
## resolution/PAratio/seed combinations have been run so far. NOT the
## leave-populations-out comparison - see compare_split_population_group.R
## and its results/tables/*-splitcompare-population_group_sdm.csv output for
## that.
##
## Run with the working directory set to scripts/ (this project's usual
## convention) - e.g. `Rscript First_modelling/summarize_iter_eval_metrics.R`.

suppressPackageStartupMessages(library(dplyr))

results_dir <- file.path('..', 'results', 'tables')

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

# Some files predate a mid-project modeller() change that added brier_class/
# cox_intercept/cox_slope rows and a std_error column (see
# compare_split_population_group.R's 10m vs 30m/90m runs for a concrete case)
# - only pr_auc/roc_auc are pulled here since those are the only two metrics
# every eval file has.
read_one <- function(f){
  d <- tryCatch(read.csv(f), error = function(e) NULL)
  if(is.null(d) || !all(c('metric', 'estimate') %in% names(d))) return(NULL)
  d <- d[d$metric %in% c('pr_auc', 'roc_auc'), c('metric', 'estimate')]
  if(nrow(d) == 0) return(NULL)
  cbind(d, parse_fname(f))
}

all_metrics <- do.call(rbind, lapply(files, read_one))

summary_by_res_iter <- all_metrics |>
  dplyr::group_by(resolution, iteration, metric) |>
  dplyr::summarize(mean = mean(estimate), sd = sd(estimate),
                    min = min(estimate), max = max(estimate), n = dplyr::n(), .groups = 'drop') |>
  dplyr::arrange(iteration, resolution, metric)

summary_by_res_iter_pa <- all_metrics |>
  dplyr::group_by(resolution, iteration, PAratio, metric) |>
  dplyr::summarize(mean = mean(estimate), sd = sd(estimate),
                    min = min(estimate), max = max(estimate), n = dplyr::n(), .groups = 'drop') |>
  dplyr::arrange(iteration, resolution, PAratio, metric)

out_by_res_iter    <- file.path(results_dir, 'iter0_iter1_eval_summary_by_resolution.csv')
out_by_res_iter_pa <- file.path(results_dir, 'iter0_iter1_eval_summary_by_resolution_PAratio.csv')
write.csv(summary_by_res_iter, out_by_res_iter, row.names = FALSE)
write.csv(summary_by_res_iter_pa, out_by_res_iter_pa, row.names = FALSE)

message('\n=== PR-AUC / ROC-AUC by resolution x iteration (across every PAratio/seed run) ===')
print(as.data.frame(summary_by_res_iter))
message('\n=== PR-AUC / ROC-AUC by resolution x iteration x PAratio ===')
print(as.data.frame(summary_by_res_iter_pa))
message('\nWritten to: ', out_by_res_iter)
message('Written to: ', out_by_res_iter_pa)
