## Probability surfaces for the "far" median models (select_far_ratio.R), which highPA_grid.R
## fit with predict_surface = FALSE. Same prediction as modeller(): predict both class
## probabilities, keep layer 2 (presence) as {fname}-Pr.tif in the model's suitability_maps/.
## 1/3 arc-second is predicted in 2 x 2 tiles and mosaicked, as modeller() does; the tiled
## predictor stack is cut once and shared by all 1/3 arc models, then deleted.
## Targets: 1% (primary) and 5% landscape prevalence. Spatial 1/3 arc at 5% has no
## admissible ratio, so there is no model for it. No AOA (too slow, and not needed after
## the near/far analysis was dropped; see near_far_evaluation.R). Models whose -Pr.tif
## already exists are skipped, so the script can be restarted.
##
## Usage (working directory scripts/): Rscript First_modelling/predict_far_surfaces.R

suppressPackageStartupMessages({
  library(terra)
  library(ranger)
})
source('functions.R')

far <- read.csv(file.path(PROJ_ROOT, 'results_highPA_fieldabs', 'tables', 'far-ratio-median-models.csv'))
far <- far[far$target_prevalence %in% c(0.01, 0.05), ]
far$root <- ifelse(far$design == 'spatial', 'results_highPA_fieldabs', 'results_highPA_classicsplit_fieldabs')
far <- far[order(match(far$resolution, c('3arc', '1arc', '1-3arc')), far$target_prevalence), ] # fastest first
far$fname <- sprintf('%s-Iteration1-PA1:%sDO:0-Seed%d', far$resolution, far$chosen_ratio, far$median_replicate)
far <- far[!duplicated(paste(far$root, far$fname)), ] # one ratio can be chosen at both targets

p2proc <- file.path('..', 'data', 'spatial', 'processed')
pr <- function(...) predict(..., type = 'response', num.threads = 1)$predictions
tile_dir <- file.path(PROJ_ROOT, 'results_highPA_fieldabs', 'suitability_maps', 'tilesPR_1-3arc')

for(i in seq_len(nrow(far))){
  m <- far[i, ]
  pout <- file.path(PROJ_ROOT, m$root, 'suitability_maps')
  out <- file.path(pout, paste0(m$fname, '-Pr.tif'))
  if(file.exists(out)){ message('Already done: ', m$root, '/', m$fname); next }

  t0 <- Sys.time()
  message(format(t0), ' - ', m$root, '/', m$fname, ' (target ', m$target_prevalence, ')')
  rf_model <- readRDS(file.path(PROJ_ROOT, m$root, 'models', paste0(m$fname, '.rds')))
  rast_dat <- rastReader(paste0('dem_', m$resolution), p2proc)

  if(m$resolution == '1-3arc'){
    if(!dir.exists(tile_dir)){
      message('  making tiles')
      dir.create(tile_dir)
      template <- rast(nrows = 2, ncols = 2, extent = ext(rast_dat), crs = crs(rast_dat))
      terra::makeTiles(rast_dat, template, filename = file.path(tile_dir, 'tile_.tif'))
    }
    tiles <- list.files(tile_dir, full.names = TRUE)
    pr_tiles <- file.path(pout, paste0('pr_tiles_', m$fname))
    dir.create(pr_tiles, showWarnings = FALSE)
    for(j in seq_along(tiles)){
      f <- file.path(pr_tiles, paste0(j, '.tif'))
      if(!file.exists(f)) terra::predict(rast(tiles[j]), rf_model, f = pr, filename = f, overwrite = TRUE, na.rm = TRUE)
      message('  tile ', j, '/', length(tiles))
      gc(verbose = FALSE)
    }
    mos <- terra::mosaic(sprc(list.files(pr_tiles, full.names = TRUE)), fun = 'mean')
    writeRaster(mos[[2]], out, wopt = list(names = 'Probability'), overwrite = TRUE)
    unlink(pr_tiles, recursive = TRUE)
  } else {
    tmp <- file.path(pout, paste0(m$fname, '-both-classes.tif'))
    terra::predict(rast_dat, rf_model, f = pr, cores = 1, cpkgs = 'ranger',
                   filename = tmp, na.rm = TRUE, overwrite = TRUE)
    writeRaster(rast(tmp)[[2]], out, wopt = list(names = 'Probability'), overwrite = TRUE)
    file.remove(tmp)
  }
  message(sprintf('  done in %.1f min', as.numeric(difftime(Sys.time(), t0, units = 'mins'))))
  rm(rf_model); gc(verbose = FALSE)
}

# the shared 1/3 arc predictor tiles are scratch once every 1/3 arc model is done
if(all(file.exists(file.path(PROJ_ROOT, far$root, 'suitability_maps', paste0(far$fname, '-Pr.tif')))))
  unlink(tile_dir, recursive = TRUE)
