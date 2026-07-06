library(terra)

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')
# ensure files and their orders match
p2thresh <- file.path('..', 'results', 'evaluations')
rasts <- file.path('..', 'results', 'suitability_maps')

f_rasts <- list.files(rasts, pattern =  'Pr')
f_thresh <- list.files(p2thresh, pattern = 'thresh.rds')

r_name <- gsub('-Pr.*$', '', f_rasts); t_name <- gsub('-thresh.*', '', f_thresh)
r_name <- r_name[r_name%in%t_name]; t_name <- t_name[t_name%in%r_name]
r_name <- r_name[order(match(r_name,t_name))] 

f_rasts <- file.path(rasts, paste0(r_name, '-Pr.tif'))
f_thresh <- file.path(p2thresh, paste0(t_name, '-thresh.rds'))

rm(r_name, t_name, p2thresh, rasts)

## now create three raster layers per stack. 
# 1) of spec_sens, 2) equal_sense_spec, and 3) sensitivity


threshold_lyrs <- function(r, t, metrics, path){

  fname <- paste0(gsub('-thresh.*$', '', basename(t)), '-thresholds.tif')
  tv <- readRDS(t)
  r <- terra::rast(r)
  
  r_thresh <- terra::rast(
    terra::ext(r), 
    resolution=terra::res(r), 
    crs = terra::crs(r), 
    nlyrs = 3,
    names = metrics
    )
  
  for(i in seq_along(metrics)){
    r_thresh[[i]] <- terra::ifel(r > tv[[ metrics[[i]] ]], 1, NA)
  }
  
  # now connect this back to disk and start porting values over. 
  terra::writeRaster(
    r_thresh, 
    filename = file.path(path, fname),
    names = metrics,
    overwrite = TRUE, 
    datatype = 'INT1U'
  )
  
  gc()
  
}

p <- file.path('..', 'results', 'threshold_masks')
metrics <- c('spec_sens', 'equal_sens_spec', 'sensitivity')


purrr::map2(threshold_lyrs, .x =  f_rasts, .y = f_thresh, metrics = metrics, path = p)
