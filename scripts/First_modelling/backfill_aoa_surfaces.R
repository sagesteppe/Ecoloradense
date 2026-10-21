## Rebuild the full area-of-applicability outputs (trainDI, and DI / AOA / LPD surfaces) for
## the iteration-1 production ("near") models, which were fit before aoaSurfaces() kept
## anything beyond the binary AOA layer. No refitting: reads each model's saved caret
## object. The rewritten single-layer {fname}-AOA.tif is identical in value to the old one.
## Fastest first (3 arc-second), so partial runs still leave complete resolutions.
## Started for the near/far analysis (dropped 2026-10-09); kept because the surfaces are
## useful in their own right. Resumable: models with surfaces already written are skipped.
##
## Usage (working directory scripts/): Rscript First_modelling/backfill_aoa_surfaces.R

suppressPackageStartupMessages({
  library(terra)
  library(CAST)
})
source('functions.R')

near_models <- data.frame(
  root  = c('results', 'results_classicsplit', 'results', 'results_classicsplit', 'results', 'results_classicsplit'),
  res   = c('3arc', '3arc', '1arc', '1arc', '1-3arc', '1-3arc'),
  seed  = c(1017, 1013, 1059, 1020, 1021, 1061)
)
p2proc <- file.path('..', 'data', 'spatial', 'processed')

for(i in seq_len(nrow(near_models))){
  m <- near_models[i, ]
  fname <- fname_build(m$res, 1, 1, 0, m$seed)
  root <- file.path(PROJ_ROOT, m$root)
  if(file.exists(file.path(root, 'suitability_maps', paste0(fname, '-AOA-surfaces.tif')))){
    message('Already done: ', m$root, '/', fname); next
  }
  t0 <- Sys.time()
  message(format(t0), ' - ', m$root, '/', fname)
  model <- readRDS(file.path(root, 'modelsTune', paste0(fname, '.rds')))
  aoaSurfaces(model, rastReader(paste0('dem_', m$res), p2proc), fname, root, LPD = TRUE)
  message(sprintf('  done in %.1f min', as.numeric(difftime(Sys.time(), t0, units = 'mins'))))
  rm(model); gc(verbose = FALSE)
}
