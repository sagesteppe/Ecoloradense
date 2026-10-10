## Probability surfaces for the "far" median models (select_far_ratio.R), which highPA_grid.R
## fit with predict_surface = FALSE. Same single-tile path as modeller(): predict both class
## probabilities, keep layer 2 (presence) as {fname}-Pr.tif in the model's suitability_maps/.
## 3 and 1 arc-second only - 1/3 arc-second not run. No AOA (too slow, and not needed after
## the near/far analysis was dropped; see near_far_evaluation.R).
##
## Usage (working directory scripts/): Rscript First_modelling/predict_far_surfaces.R

suppressPackageStartupMessages({
  library(terra)
  library(ranger)
})
source('functions.R')

far <- read.csv(file.path(PROJ_ROOT, 'results_highPA_fieldabs', 'tables', 'far-ratio-median-models.csv'))
far <- far[far$resolution %in% c('3arc', '1arc'), ]
far$root <- ifelse(far$design == 'spatial', 'results_highPA_fieldabs', 'results_highPA_classicsplit_fieldabs')
far <- far[order(far$resolution != '3arc'), ] # fastest first

pr <- function(...) predict(..., type = 'response', num.threads = 1)$predictions

for(i in seq_len(nrow(far))){
  m <- far[i, ]
  fname <- sprintf('%s-Iteration1-PA1:%sDO:0-Seed%d', m$resolution, m$chosen_ratio, m$median_replicate)
  pout <- file.path(PROJ_ROOT, m$root, 'suitability_maps')
  out <- file.path(pout, paste0(fname, '-Pr.tif'))
  if(file.exists(out)){ message('Already done: ', m$root, '/', fname); next }

  t0 <- Sys.time()
  message(format(t0), ' - ', m$root, '/', fname)
  rf_model <- readRDS(file.path(PROJ_ROOT, m$root, 'models', paste0(fname, '.rds')))
  rast_dat <- rastReader(paste0('dem_', m$resolution), file.path('..', 'data', 'spatial', 'processed'))
  tmp <- file.path(pout, paste0(fname, '-both-classes.tif'))
  terra::predict(rast_dat, rf_model, f = pr, cores = 1, cpkgs = 'ranger',
                 filename = tmp, na.rm = TRUE, overwrite = TRUE)
  writeRaster(rast(tmp)[[2]], out, wopt = list(names = 'Probability'), overwrite = TRUE)
  file.remove(tmp)
  message(sprintf('  done in %.1f min', as.numeric(difftime(Sys.time(), t0, units = 'mins'))))
  rm(rf_model); gc(verbose = FALSE)
}
