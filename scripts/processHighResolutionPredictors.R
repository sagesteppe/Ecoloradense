library(tidyverse)
library(sf)
library(terra)

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')

# first we will create a domain for all analysis. The 'closest' this bounding box is to a known 
# occurrence is 10 miles. The furthest distances vary. 
domain <- sf::st_read('../data/collections/occurrences_coloradense/occurrences.shp', quiet = T) %>% 
  st_union() %>% 
  st_transform(32613) %>% 
  st_buffer(16093) %>% 
  st_transform(4326) %>% 
  vect() 

ext(domain)
template <- rast(project(domain, 'EPSG:32613'), nrows = 5, ncols = 5)


# the 1m resolution data set is very large. We need to process it in batches
# to achieve the desired effect. 

project1m <- function(x){
  
  paths <- file.path(x, list.files(x, 'tif$'))
  
  for (i in seq_along(paths)){
    ob <- terra::project(terra::rast(paths[i]), 'EPSG:32613')
    terra::writeRaster(ob, gsub('USGS_1M_|USGS_one_meter_', '', paths[i]), overwrite = T)
    unlink(paths[i])
    gc()
  }
}
# project1m('../data/spatial/raw/dem_1m') # we ran this already.


# now we will create 'tiles' of the data set. 

tiles1m <- function(x, domain){
  
  paths <- file.path(x, list.files(x, 'tif$'))
  dir.create(tile_p <- file.path(x, 'tiles'))
  
  DEMs <- terra::sprc(paths); message('sprc assembled')
  DEMs <- terra::mosaic(DEMs); message('mosaic complete')
  
  names(DEMs) <- 'elevation'; message('Name set on layers')
  makeTiles(DEMs, y = template, filename = file.path(x, 'tiles', 'DEM_.tif'))
}

tiles1m('../data/spatial/raw/dem_1m', domain = domain)
# now crop to domain

tiles1m <- function(x, domain){
  
  paths <- file.path(x, list.files(x, 'tif$'))
  dir.create(tile_p <- file.path(x, 'tiles'))
  
  DEMs <- terra::sprc(paths); message('sprc assembled')
  DEMs <- terra::mosaic(DEMs); message('mosaic complete')
  
  # d <- terra::project(domain, crs(DEMs)) |>
  #   terra::ext()
  #  DEMs <- terra::crop(DEMs, d)
  names(DEMs) <- 'elevation'; message('Name set on layers')
  makeTiles(DEMs, y = template, filename = file.path(x, 'tiles', 'DEM_.tif'))
}

tiles1m('../data/spatial/raw/dem_1m', domain = domain)

# now we will resample these tiles to 3m resolution - usgs does not have a product
# at this resolution

tiles1m_crop <- function(x, domain){
  
  DEMs <- vrt(file.path(x, list.files(x, 'tif$')))
  
  d <- terra::project(domain, crs(DEMs)) |>
    terra::ext()
  DEMs <- terra::crop(DEMs, d)
  makeTiles(DEMs, y = template, filename = file.path(x, 'DEMc_.tif'))
}

tiles1m_crop('../data/spatial/raw/dem_1m/tiles', domain = domain)

tiles3m_make <- function(x, path_out){
  
  DEMs <- terra::vrt(file.path(x, list.files(x, 'tif$')))
  DEMs <- terra::aggregate(DEMs, fact = c(3, 3))
  terra::makeTiles(DEMs, y = template, filename = file.path(path_out, 'DEM_.tif'))
}

tiles3m_make(x = '../data/spatial/raw/dem_1m/tiles',
             path_out = '../data/spatial/raw/dem_3m/')

f <- paste0('../data/spatial/raw/dem_3m/', list.files('../data/spatial/raw/dem_3m/'))
dem3m <- vrt(f)
writeRaster(dem3m, '../data/spatial/processed/dem_3m/dem.tif')
