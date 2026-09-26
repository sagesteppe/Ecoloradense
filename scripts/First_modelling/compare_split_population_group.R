## Leave-2-3-populations-out test/train split for the iteration-1 SDM
## (presence/absence RF classifier), emulating Census_uncertainty/
## compare_split_strategies.R's 'population_group' method for the count/
## density models - see functions.R's population_group_replicates_sdm()/
## population_group_split() for how it's adapted to SDM data (only the
## presence + paired field near-absence rows carry a real Lctn_bb population
## label; the shared background pseudo-absence pool doesn't, and is folded
## into a replicate's test set only via a spatial buffer around the held-out
## population's presence points - see those functions' docs for the full
## rationale).
##
## Cost: one distOrder_PAratio_simulator() call per replicate, with
## predict_surface = FALSE (fit + holdout eval only - Boruta + knndm CV +
## ranger tuning, no raster surface/AOA/SE prediction), matching how
## adaptive_PAratio_search() runs its own grid/replicate search. Exhaustive
## over population_group_replicates_sdm()'s replicates (no racing/stopping
## rule - unlike compare_split_strategies.R's stochastic 'twinning'/'classic'
## methods, there's nothing to race here).
##
## Results land under a separate `results_populationgroupsplit/` tree (never
## `results/`), so this never collides with or overwrites the production
## iteration-1 models. Run with the working directory set to scripts/ (this
## project's usual convention) - e.g.
## `Rscript First_modelling/compare_split_population_group.R [resolution]`.

suppressPackageStartupMessages({
  library(tidyverse); library(sf); library(terra)
  library(caret); library(ranger); library(Boruta); library(CAST); library(dismo)
})

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')
source('functions.R')

results_root <- file.path(PROJ_ROOT, 'results_populationgroupsplit')
for(d in c('models', 'modelsTune', 'evaluations', 'test_data', 'suitability_maps', 'tables')){
  dir.create(file.path(results_root, d), recursive = TRUE, showWarnings = FALSE)
}

# ---- resolution taken from the command line (numeric: 3, 10, 30, or 90) ----
# defaults to 10 (1/3 arc-second, 10m) - the same scale the count-model
# 'population_group' comparison used (1-3arc-Iteration1-PA1:1...).
args <- commandArgs(trailingOnly = TRUE)
resolution <- if(length(args) >= 1) as.numeric(args[[1]]) else 10
res_file <- switch(as.character(resolution),
  '90' = '90m-presence-iter1.gpkg',
  '30' = '30m-presence-iter1.gpkg',
  '10' = '3m-presence-iter1.gpkg',
  stop("resolution must be one of 3, 10, 30, 90 (got ", resolution, ")"))

## ---- build the combined presence+absence data, RETAINING Lctn_bb ----
## (First_modelling/Modelling.Rmd's own m90/m30/m10 immediately
## `select(Occurrence)`, dropping Lctn_bb - this comparison needs it kept.)
labeled <- sf::st_read(file.path('..', 'data', 'Data4modelling', res_file), quiet = TRUE) |>
  dplyr::rename(Occurrence = Presenc) |>
  sf::st_as_sf()
labeled <- labeled[sf::st_is(labeled, 'POINT'), ]
labeled <- dplyr::select(labeled, Occurrence, Lctn_bb)

bg_abs <- sf::st_read(file.path('..', 'data', 'Data4modelling', 'iter1-pa.gpkg'), quiet = TRUE) |>
  dplyr::mutate(Lctn_bb = NA_character_) |>
  dplyr::select(Occurrence, Lctn_bb)

x <- dplyr::bind_rows(labeled, bg_abs) |>
  dplyr::mutate(.id = seq_len(dplyr::n()))

message(sprintf('resolution = %s (%s): %d labeled rows (%d presence, %d near-absence), %d background absences',
                 resolution, res_string(resolution),
                 nrow(labeled), sum(labeled$Occurrence == 1), sum(labeled$Occurrence == 0),
                 nrow(bg_abs)))

replicates <- population_group_replicates_sdm(labeled)
message(sprintf('%d population_group replicates: %s', length(replicates), paste(names(replicates), collapse = ', ')))

bn_tag <- paste0(res_string(resolution), '-Iteration1-PA1:1')
out_csv <- file.path('..', 'results', 'tables', paste0(bn_tag, '-splitcompare-population_group_sdm.csv'))

## ---- run one replicate through distOrder_PAratio_simulator(), incrementally saving ----
run_replicate <- function(lbl, replicate, seed){
  message(sprintf('\n===== population_group_sdm replicate %s (seed = %d) =====', lbl, seed))
  t0 <- Sys.time()

  split <- population_group_split(x, replicate = replicate)

  res <- distOrder_PAratio_simulator(
    x, distOrder = 0, PAratio = 1, resolution = resolution, seed = seed,
    split = split, predict_surface = FALSE, results_root = results_root,
    p2proc = file.path('..', 'data', 'spatial', 'processed'),
    iteration = 1, train_split = 0.8, se_prediction = FALSE
  )

  message(sprintf('[timing] replicate %s total: %.1fs', lbl, as.numeric(difftime(Sys.time(), t0, units = 'secs'))))

  metrics <- res$metrics
  metrics$replicate <- lbl
  metrics$n_test_presence <- split$summary$n_test[split$summary$class == 'presence']
  metrics$n_test_absence  <- split$summary$n_test[split$summary$class == 'absence']
  metrics
}

all_results <- if(file.exists(out_csv)){
  list(readr::read_csv(out_csv, show_col_types = FALSE))
} else list()
done <- if(length(all_results)) unique(all_results[[1]]$replicate) else character(0)

save_progress <- function(){
  combined <- dplyr::bind_rows(all_results)
  write.csv(combined, out_csv, row.names = FALSE)
  combined
}

base_seed <- 2021
for(i in seq_along(replicates)){
  lbl <- names(replicates)[i]
  if(lbl %in% done) next
  all_results[[length(all_results) + 1]] <- run_replicate(lbl, replicates[[i]], seed = base_seed + i)
  save_progress()
}

final <- save_progress()

message('\n=== FINAL (population_group, SDM, resolution = ', resolution, ') ===')
print(as.data.frame(final))
message('\nWritten to: ', out_csv)
