## Score highPA_grid.R's leave-one-population-out folds (field-absence design) by absence
## type: each fold's model is scored on its whole held-out zone, split into
##   presences vs field absences     - the near-population decision
##   presences vs background         - finding a new population in the landscape
##   the same, reweighted to 1%      - at a realistic landscape prevalence
## Test zones are identical at every LOPO ratio (ratio_train_only = TRUE), so ratios are
## compared on the same data.
##
## Usage (working directory scripts/): Rscript First_modelling/score_lopo_by_abs_type.R 90
## Writes results_highPA_populationsplit_fieldabs/tables/{res}-LOPO-by-absence-type.csv.

suppressPackageStartupMessages({
  library(sf)
  library(tidyverse)
  library(terra)
  library(ranger)
})
source('functions.R')

args <- commandArgs(trailingOnly = TRUE)
resolution <- as.numeric(args[1])
stopifnot('first argument must be a resolution: 90, 30 or 10' = resolution %in% c(90, 30, 10))
res_str <- res_string(resolution)
lopo_ratios <- c(4, 20)
root <- file.path(PROJ_ROOT, 'results_highPA_populationsplit_fieldabs')
p2dat <- file.path('..', 'data', 'Data4modelling')

# the modelling points, built as highPA_grid.R builds them
abs <- bind_rows(
  st_read(file.path(p2dat, 'iter1-pa.gpkg'), quiet = TRUE) |> select(Occurrence),
  st_read(file.path(p2dat, 'iter1-pa-highPA-reviewed.gpkg'), quiet = TRUE) |> select(Occurrence)
)
pres_file <- c('90' = '90m', '30' = '30m', '10' = '3m')[[as.character(resolution)]]
x <- st_read(file.path(p2dat, paste0(pres_file, '-presence-iter1.gpkg')), quiet = TRUE) |>
  rename(Occurrence = Presenc) |>
  mutate(field_abs = Occurrence == 0)
x <- bind_rows(x, mutate(abs, field_abs = FALSE))
x <- filter(x, st_is(x, 'POINT')) |> select(Occurrence, field_abs) |> mutate(.id = seq_len(n()))

# same zoning population_lopo_search() uses
areas <- st_read(file.path('..', 'data', 'GroundTruthPts', 'ground_truth_median-VISITED.gpkg'),
                 layer = 'areas', quiet = TRUE) |>
  st_transform(st_crs(x))
x <- assign_population(x, areas, max_dist = 500)
covs <- terra::extract(rastReader(paste0('dem_', res_str), file.path('..', 'data', 'spatial', 'processed')), vect(x))

score <- function(d, tp = NULL){
  if(length(unique(d$truth)) < 2) return(tibble(roc_auc = NA_real_, brier_scaled = NA_real_, cox_intercept = NA_real_))
  presenceScores(d$truth, d$p, target_prevalence = tp)$wide |> select(roc_auc, brier_scaled, cox_intercept)
}

out <- map_dfr(lopo_ratios, function(r){
  pp <- read.csv(file.path(root, 'tables', paste0(res_str, '-Iteration1-LOPO-PA', r, 'DO:0.csv')))
  map_dfr(seq_len(nrow(pp)), function(i){
    if(is.na(pp$pr_auc[i])) return(NULL) # fold population_lopo_search() had to skip
    # fold i is seeded 6000 + i (population_lopo_search(), base_seed = 6000)
    m <- readRDS(file.path(root, 'models', sprintf('%s-Iteration1-PA1:%sDO:0-Seed%d.rds', res_str, r, 6000 + i)))
    z <- which(x$population_zone == pp$population[i] & complete.cases(covs))
    d <- tibble(truth = x$Occurrence[z], field = x$field_abs[z],
                p = predict(m, data = covs[z, ])$predictions[, '1'])
    base <- function(subset, n_a) tibble(resolution = res_str, ratio = r, population = pp$population[i],
                                         subset = subset, n_p = sum(d$truth == 1), n_a = n_a)
    bind_rows(
      bind_cols(base('field abs vs presences', sum(d$field)), score(filter(d, truth == 1 | field))),
      bind_cols(base('background vs presences', sum(!d$field & d$truth == 0)), score(filter(d, truth == 1 | !field))),
      bind_cols(base('background vs presences @1%', sum(!d$field & d$truth == 0)), score(filter(d, truth == 1 | !field), 0.01)))
  })
})
write.csv(out, file.path(root, 'tables', paste0(res_str, '-LOPO-by-absence-type.csv')), row.names = FALSE)

options(width = 200)
message('Presence-weighted means by subset and ratio:')
out |>
  filter(!is.na(roc_auc)) |>
  group_by(subset, ratio) |>
  summarise(n_pops = n(), roc = weighted.mean(roc_auc, n_p), brier_scaled = weighted.mean(brier_scaled, n_p),
            cox_int_median = median(cox_intercept), .groups = 'drop') |>
  as.data.frame() |>
  print(digits = 2)
