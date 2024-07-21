# SDM analyses will be performed at 4 spatial resolutions. 
# 90m or 3 arc seconds, 30m or 1 arc second, 10m or 1/3 arc second, and 3 meters. 
# These values will be compared via ground truthing, generating an interesting
# albeit small intermediary paper, how do these different resolution data perform 
# for predicting presence and absence? 

# In order for the most possible comparisons of these products, we will generate
# each predictor from a DEM native to the above generation using the same functions
# and processes. 

# Predictions from the more coarse grain products will be resampled to 3m resolution
# for equivalency in ground truthing - i.e. equal probabilities of detection of individuals!

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
# now we will assemble single DEM tifs which cover our domain

#' assmeble and crop DEM's' to an area. 
#' @param x path to directory contains tifs to combine for final product
#' @param domain vector data which to crop the extent to
DEMcrop <- function(x, domain){
  
  paths <- file.path(x, list.files(x, 'tif$'))
  
  if(length(paths) > 1){
    DEMs <- terra::sprc(paths)
    DEMs <- terra::mosaic(DEMs) } else {
      DEMs <- terra::rast(file.path(x, list.files(x, 'tif$')))}
  
  d <- terra::project(domain, crs(DEMs)) |>
    terra::ext()
  DEMs <- terra::crop(DEMs, d)
  names(DEMs) <- 'elevation'
  
  proj_dem_path <- file.path(
    '../data/spatial/processed', basename(x), 'dem.tif')
  
  terra::project(DEMs, "epsg:32613", filename = proj_dem_path, overwrite = TRUE)

}

DEMcrop('../data/spatial/raw/dem_3arc', domain = domain)
DEMcrop('../data/spatial/raw/dem_1arc', domain = domain)
DEMcrop('../data/spatial/raw/dem_1-3arc', domain = domain)
DEMcrop('../data/spatial/raw/dem_3m', domain = domain)

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
project1m('../data/spatial/raw/dem_1m')


# now we will create 'tiles' of the data set. 

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

#' morphoMaker
#' calculate many geomorphological features of a landscape
#' @param x a path to a digital elevation model
#' @param p an output path where all geomorphologic products should be saved. 
morphoMaker <- function(x, p){
  
  if(str_detect(x, 'tif$')){demIN <- x} else{
      demIN <- terra::mosaic(terra::sprc(file.path(x, list.files(x, 'tif$'))))
  }
  dir.create(p)
  whitebox::wbt_fill_depressions(demIN, output = file.path(p, 'DEM_filled.tif'))

  whitebox::wbt_aspect(file.path(p, 'DEM_filled.tif'), output = file.path(p, 'aspect.tif'))
  whitebox::wbt_slope(file.path(p, 'DEM_filled.tif'), output = file.path(p, 'slope.tif'))
  whitebox::wbt_ruggedness_index(file.path(p, 'DEM_filled.tif'), output = file.path(p, 'ruggedness.tif'))
  whitebox::wbt_geomorphons(file.path(p, 'DEM_filled.tif'), output = file.path(p, 'geomorphons.tif'))
  whitebox::wbt_pennock_landform_class(file.path(p, 'DEM_filled.tif'), output = file.path(p, 'Pennock.tif'))
  whitebox::wbt_relative_topographic_position(demIN, output = file.path(p, 'RTP.tif'))

  whitebox::wbt_downslope_index(file.path(p, 'DEM_filled.tif'), output = file.path(p, 'DSI.tif'))
  whitebox::wbt_profile_curvature(file.path(p, 'DEM_filled.tif'), output = file.path(p, 'profile_curv.tif'))
  whitebox::wbt_plan_curvature(file.path(p, 'DEM_filled.tif'), output = file.path(p, 'plan_curv.tif'))
  whitebox::wbt_maximal_curvature(file.path(p, 'DEM_filled.tif'), output = file.path(p, 'maximal_curv.tif'))
  whitebox::wbt_tangential_curvature(file.path(p, 'DEM_filled.tif'), output = file.path(p, 'tangential_curv.tif'))
  
  whitebox::wbt_d8_pointer(file.path(p, 'DEM_filled.tif'), esri_pntr = F,
                           output = file.path(p, 'D8pntr.tif'))
  whitebox::wbt_basins(file.path(p, 'D8pntr.tif'), esri_pntr = F, output = file.path(p, 'basins.tif')) 
  
}

morphoMaker(x = '../data/spatial/raw/dem_1m/tiles', 
            p = '../data/spatial/processed/dem_1m/geomorphology')




# cost surfaces between known occurrences can be generated via modelling
?wbt_hydrologic_connectivity() # how well connected are drainages?
?wbt_trace_downslope_flowpaths # where does alluvial drive dispersal go?

# and then the costs can be analyzed as such 
?wbt_cost_distance()
?wbt_cost_allocation()


# The NAIP data were processed as so, using bash, and leaving them at their native
# resolution 


# it unfortunately requires that we extract the .sid files from the .zip files!!! 
cd /media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/raw/NAIP/4band
  
for file in *zip; do

in_path='/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/raw/NAIP/4band/'
out_path='/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed/NAIP/'
infilename="$in_path${file%.zip}/${file%.zip}/${file%.zip}.sid"
outfilename="$out_path${file%.zip}.tif"

mkdir ${file%.zip}
unzip $file -d ${file%.zip}
~/MrSID_DSDK-9.5.1.4427-linux.x86-64.gcc48/Raster_DSDK/bin/mrsiddecode -i $infilename -o  $outfilename -s 1
rm -r ${file%.zip} 

now=$(date)
printf "${file%.zip} was processed at: $now\n"

done



################################################################################
## we will write these data to generate rasters using ClimateNA for the study area ##
# we will generate projections at 3arc. 

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')

## these will need to be in WGS 84, the same CRS as (epsg 4326) as PRISM data ##
project(rast('../data/spatial/processed/dem_3arc/dem.tif'), "EPSG:4326", threads = 16, 
        filename = '../data/spatial/processed/dem_3arc/dem_3arc-wgs84.asc', NAflag = -9999)
subst(
  rast('../data/spatial/processed/dem_3arc/dem_3arc.asc'), NA, -9999, 
  filename = '../data/spatial/processed/dem_3arc/dem_3arc.asc')

# and we will make the file work with the climateNA executable on windows.. 

temp <- readLines('../data/spatial/processed/dem_3arc/dem_3arc.asc')
write.table(temp,'../data/spatial/processed/dem_3arc/dem_3arc.asc',
            row.names=F,col.names=F,quote=F)

# once climateNA processes this file it will simply be 'cut' into smaller cells 
# for the remaining variables. 


############# Resample ClimateNA data for each Resolution ########################
# we ran climateNA on ~90m resolution data, which is roughly 11% of the native resolution
# of prism data. We will bring these data into the lower resolutions so they can be used
# with each modelling approach. 

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed/climateNA_Normal_1961_1990Y-3arc')
f <- rast(list.files())
crs(f) <- 'EPSG:4326'
f <- project(f, crs(arc3_template))

# write out the 3arc data here. 
resample(f, arc3_template, threads = 16, method = 'bilinear',
         filename = '../../processed/dem_3arc/climateNA.tif', overwrite = T)

resample(f, 
         rast('/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed/dem_1arc/dem.tif'), 
         threads = 16, method = 'bilinear',
         filename = '../../processed/dem_1arc/climateNA.tif', overwrite = T)

resample(f, 
         rast('/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed/dem_1-3arc/dem.tif'), 
         threads = 16, method = 'bilinear',
         filename = '../../processed/dem_1-3arc/climateNA.tif', overwrite = T)

resample(f, 
         rast('/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed/dem_3m/dem.tif'), 
         threads = 16, method = 'bilinear',
         filename = '../../processed/dem_3m/climateNA.tif', overwrite = T)

############ naip create products at each resolution ########
setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed/')

f <- vrt(file.path('NAIP', list.files('NAIP', pattern = 'tif')))
arc3_template <- rast(
  '/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed/dem_3arc/dem.tif')
arc1_template <- rast(
  '/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed/dem_1arc/dem.tif')
arc13_template <- rast(
  '/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed/dem_1-3arc/dem.tif')
m3_template <- rast(
  '/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed/dem_3m/dem.tif')

terra::resample(f, arc3_template, threads = 16, filename = './NAIP/3arc/NAIP.tif') # 6:40 start
terra::resample(f, arc1_template, threads = 16, filename = './NAIP/1arc/NAIP.tif') # 9:00 start
terra::resample(f, arc13_template, threads = 16, filename = './NAIP/1-3arc/NAIP.tif') # about 2 hours
terra::resample(f, m3_template, threads = 16, filename = './NAIP/3m/NAIP.tif') # about 2 hours

############ we will also create a GLCM texture for the NAIP imagery. ##########

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed/')
f <- vrt(file.path('NAIP', list.files('NAIP', pattern = 'tif')))

glcm <- glcm::glcm(raster::raster(rast('NAIP/3arc/NAIP.tif')), # 30 minutes or so
                   window = c(5,5), shift = list(c(0,1), c(1,1), c(1,0), c(1,-1)), 
                          min_x = 1, max_x = 255, na_opt = 'ignore', na_val = 0)
writeRaster(glcm, 'NAIP/3arc/GLCM.tif')

glcm <- glcm::glcm(raster::raster(rast('NAIP/1arc/NAIP.tif')), # 1.5 or so
                   window = c(5,5), shift = list(c(0,1), c(1,1), c(1,0), c(1,-1)), 
                   min_x = 1, max_x = 255, na_opt = 'ignore', na_val = 0)
writeRaster(glcm, 'NAIP/1arc/GLCM.tif')

glcm <- glcm::glcm(raster::raster(rast('NAIP/1-3arc/NAIP.tif')), # 3 hours or so
                   window = c(5,5), shift = list(c(0,1), c(1,1), c(1,0), c(1,-1)), 
                   min_x = 1, max_x = 255, na_opt = 'ignore', na_val = 0)
writeRaster(glcm, 'NAIP/1-3arc/GLCM.tif')

glcm <- glcm::glcm(raster::raster(rast('NAIP/3m/NAIP.tif')), # 10 hours or so
                   window = c(5,5), shift = list(c(0,1), c(1,1), c(1,0), c(1,-1)), 
                   min_x = 1, max_x = 255, na_opt = 'ignore', na_val = 0)
writeRaster(glcm, 'NAIP/3m/GLCM.tif')

rm(glcm)



############# make the vegetation cover match each grain ######################

# we want to combine the trees into a single dataset for 'forest'

p2d <- '/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/raw/VegCover'
forest <- rast(file.path(p2d, 'DecidiousBroadleafTrees.tif')) + 
  rast(file.path(p2d, 'MixedTrees.tif')) +
  rast(file.path(p2d, 'NeedleleafTrees.tif'))
shrubs <- rast(file.path(p2d, 'Shrubs.tif'))
herbs <- rast(file.path(p2d, 'HerbaceousVegetation.tif'))

veg <- c(forest, shrubs, herbs)
veg <- crop(veg, ext(project(arc3_template, crs(veg))))
veg <- project(veg, crs(arc3_template))
names(veg) <- c('Tree', 'Shrub', 'Herbaceous')

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed')

terra::resample(veg, arc3_template, threads = 16, filename = './dem_3arc/Vegetation.tif') 
terra::resample(veg, arc1_template, threads = 16, filename = 'dem_1arc/Vegetation.tif') 
terra::resample(veg, arc13_template, threads = 16, filename = 'dem_1-3arc/Vegetation.tif')
terra::resample(veg, m3_template, threads = 16, filename = 'dem_3m/Vegetation.tif')

rm(forest, shrubs, herbs, veg)



### Create Latitude and longitude layers for each resolution 

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed/')

coorder <- function(x){
  template <- terra::rast(x)
  
  Longitude <- terra::init(template, 'x') ; names(Longitude) <- 'Longitude'
  Latitude <- terra::init(template, 'y') ; names(Latitude) <- 'Latitude'
  
  coords <- c(Longitude, Latitude)
  terra::writeRaster(coords, file.path(dirname(x), 'Coordinates.tif'), overwrite = T)
}

lapply(
  c('/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed/dem_3arc/dem.tif',
  '/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed/dem_1arc/dem.tif',
  '/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed/dem_1-3arc/dem.tif',
  '/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed/dem_3m/dem.tif'),
  coorder)
