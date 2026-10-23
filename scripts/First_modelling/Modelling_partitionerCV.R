## Iter-1 low-PA re-run under the Partitioner.R CV modes.
##
## Mirrors Modelling.Rmd (adaptive_PAratio_search at 90/30/10 m, iteration = 1), but
## swaps the outer train/test split and inner CV folds from spatial_class_split() +
## CAST::knndm() to Partitioner.R's build_partition_mapping() (spatial_mapping) or
## build_partition_lmpo() (LMPO). Results land in their own trees:
##
##   Results_spatialmapping_iter1a_low_pa/   or   Results_LMPO_iter1a_low_pa/
##
## so the original results/ and results_classicsplit/ / results_populationsplit/
## stay intact for comparison.
##
## Usage (working directory scripts/):
##   Rscript First_modelling/Modelling_partitionerCV.R spatial_mapping          # all 3 res
##   Rscript First_modelling/Modelling_partitionerCV.R LMPO
##   Rscript First_modelling/Modelling_partitionerCV.R spatial_mapping 90       # one res
##
## The adaptive search's bracketing/replicate schedule is unchanged - only the CV
## plumbing swaps. See adaptive_PAratio_search() in functions.R.

suppressPackageStartupMessages({
  library(sf)
  library(tidyverse)
  library(terra)
  library(caret)
  library(ranger)
})
source('functions.R')
source('final_sdm/Partitioner.R')
source('final_sdm/partition_entrypoints.R')

args <- commandArgs(trailingOnly = TRUE)
cv_mode <- if(length(args) >= 1) args[1] else 'spatial_mapping'
stopifnot('first argument must be spatial_mapping or LMPO' =
            cv_mode %in% c('spatial_mapping', 'LMPO'))
resolutions <- if(length(args) >= 2) as.numeric(args[-1]) else c(90, 30, 10)
stopifnot('resolutions must be one or more of 90/30/10' =
            all(resolutions %in% c(90, 30, 10)))

cv_tag <- c(spatial_mapping = 'spatialmapping', LMPO = 'LMPO')[[cv_mode]]
results_root <- file.path(PROJ_ROOT, paste0('Results_', cv_tag, '_iter1a_low_pa'))
for(d in c('models', 'modelsTune', 'test_data', 'tables', 'evaluations', 'suitability_maps'))
  dir.create(file.path(results_root, d), showWarnings = FALSE, recursive = TRUE)

p2dat  <- file.path('..', 'data', 'Data4modelling')
p2proc <- file.path('..', 'data', 'spatial', 'processed')

## Shared iter-1 absence pool (same file Modelling.Rmd uses) ------------------------------
abs <- st_read(file.path(p2dat, 'iter1-pa.gpkg'), quiet = TRUE)

presence_file <- c('90' = '90m-presence-iter1.gpkg',
                   '30' = '30m-presence-iter1.gpkg',
                   '10' = '3m-presence-iter1.gpkg')

for(resolution in resolutions){
  message('\n=== ', cv_mode, ' @ ', res_string(resolution), ' =============================')
  m <- st_read(file.path(p2dat, presence_file[[as.character(resolution)]]), quiet = TRUE) |>
    rename(Occurrence = Presenc) |>
    st_as_sf()
  m <- bind_rows(m, abs)
  m <- filter(m, st_is(m, 'POINT')) |> select(Occurrence)
  # adaptive_PAratio_search() tags .id/abs_k on `m` itself, so no pre-tagging needed here.

  t0 <- Sys.time()
  search <- adaptive_PAratio_search(
    x = m, distOrder = 0, resolution = resolution, iteration = 1,
    train_split = 0.8, p2proc = p2proc,
    cv_mode = cv_mode,
    results_root = results_root
  )
  mins <- as.numeric(difftime(Sys.time(), t0, units = 'mins'))
  message(sprintf('%s @ %s: best 1:%.2f in %.1f min (seed %d).',
                  cv_mode, res_string(resolution), search$best_ratio, mins, search$best_seed))

  # persist the partition that drove this resolution's search so plots / diagnostics can
  # reload it alongside the fits
  if(!is.null(search$part))
    saveRDS(search$part, file.path(results_root, 'tables',
                                   paste0(res_string(resolution), '-partition.rds')))
}
