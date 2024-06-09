# SDM analyses will be performed at 4 spatial resolutions. 
# 90m or 3 arc seconds, 30m or 1 arc second, 10m or 1/3 arc second, and 3 meters. 
# These values will be compared via ground truthing, generating an interesting
# albeit small intermediary paper, how do these different resolution data perform 
# for predicting presence and absence. 

# In order for the most possible comparisons of these products, we will generate
# each predictor from a DEM native to the above generation using the same functions
# and processes. 

# Predictions from the more coarse grain products will be resampled to 3m resolution
# for equivalency in ground truthing - i.e. equal probabilities of detection of individuals!

library(tidyverse)
library(sf)
library(terra)


# first we will create a domain for all analysis. The 'closest' this bounding box is to a known 
# occurrence is 10 miles. The furthest distances vary. 
domain <- sf::st_read('../data/collections/occurrences_coloradense/occurrences.shp', quiet = T) %>% 
  st_union() %>% 
  st_transform(32613) %>% 
  st_buffer(16093) %>% 
  st_transform(4326) %>% 
  vect() 


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

# the domain for all analyses will be the Upper Gunnison, where 3m resolution data are 
# available. 
# From the RMBL Spatial Data Platform we will download an Digital Elevation Model
# https://rmbl-sdp.s3.us-east-2.amazonaws.com/data_products/released/release3/UG_dem_3m_v1.tif


#' morphoMaker
#' calculate many geomorphological features of a landscape
#' @param x a path to a digital elevation model
#' @param p an output path where all geomorphologic products should be saved. 
morphoMaker <- function(x, p){
  
  demIN <- x
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

morphoMaker(x = '../data/spatial/processed/dem_3m/dem.tif', 
            p = '../data/spatial/processed/dem_3m/geomorphology')




# cost surfaces between known occurrences can be generated via modelling
?wbt_hydrologic_connectivity() # how well connected are drainages?
?wbt_trace_downslope_flowpaths # where does alluvial drive dispersal go?

# and then the costs can be analyzed as such 
?wbt_cost_distance()
?wbt_cost_allocation()




# The NAIP data were processed as so, using bash, and leaving them at their native
# resolution 


# it unfortunately requires that we extract the .sid files from the .zip files!!! 
cd /media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/raw/NAIP
  
  for file in *zip; do

in_path='/media/sagesteppe/ExternalHD/NAIP/raw_imagery/'
out_path='/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed/NAIP'
infilename="$in_path${file%.zip}/${file%.zip}.sid"
outfilename="$out_path${file%.zip}.tif"

mkdir ${file%.zip}
unzip $file -d ${file%.zip}
~/MrSID_DSDK-9.5.1.4427-linux.x86-64.gcc48/Raster_DSDK/bin/mrsiddecode -i $infilename -o  $outfilename -s 4
rm -r ${file%.zip} # reduce by  a factor of 4 from 60cm -> 120 -> 240 -> 960 and resample to align with DEM

now=$(date)
printf "${file%.zip} was processed at: $now\n"

done














# From NAIP data we calculate the following metrics  
# Percent fractional rock/bedrock  
# Percent soil cover   
# Canopy height model 
# NDVI 

#' NAIPerville
#' Calculate a few derived metrics from NAIP
#' @param x input files

NAIPerville <- function(x){
  
  # NDVI
  # Soil Color Index
  # percent fractional rock/bedrock
  
  
}


LiDar <- function(x){
  
  # canopy height model
  # percent veg cover
  # percent fractional rock/bedrock
}

