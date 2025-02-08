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



### Create snow free days product from Sentinel 2 imagery ###

library(CDSE)
library(terra)
library(sf)

OAuthClient <- GetOAuthClient(
  id = Sys.getenv("CopernicusOAuthID"),
  secret = Sys.getenv("CopernicusOAuthSecret")
)

OAuthClient <- GetOAuthClient(
  id = Sys.getenv("CopernicusOAuthIDNORTHWESTERN"),
  secret = Sys.getenv("CopernicusOAuthSecretNORTHWESTERN")
)

script_file <- '/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts/NDSI_download.js'
aoi <- st_as_sfc(st_bbox(template))

catalog_results <- SearchCatalog(
  aoi = aoi, 
  from = as.Date('2017-05-01'), to = as.Date('2024-08-01'),
  # first relevant flights are in 2017, a was flying in 2016,
  # but we'll just skip that year to make writing methods easier
  collection = "sentinel-2-l2a",
  client = OAuthClient
)

cr <- catalog_results |>
  mutate(Month = as.numeric(str_remove_all(str_extract(acquisitionDate, '-[0-9]{2}-'), '-'))) |>
  filter(Month >= 5, Month <= 7) |> # we want to detect late laying snow packs
  group_by(acquisitionDate) |> 
  # only bother with tiles which have good coverage on that day! No need for random tiles. 
  mutate(
    totalArea = sum(areaCoverage), 
    n = n()) |>
  filter(n > 4 & totalArea >= 100) |> 
  # NDSI cannot see through clouds, we will drop dates with
  # very high cloud cover 
  group_by(acquisitionDate) |>
  mutate(
    meanCloud = mean(tileCloudCover)
  ) |>
  filter(meanCloud <= 40) |>

  # for visualizing whether we have an OK cloud drop off time
  mutate(
    year = str_extract(acquisitionDate, '[0-9]{4}'), 
    doy = yday(acquisitionDate)
  )

ggplot() + 
  geom_point(data = cr, aes(x = doy, y = year)) 

## what if we work on modelling the high snow years.... ##

cr_high <- cr |>
  filter(year %in% c(2017, 2019, 2023))

ggplot(cr_high) + 
  geom_point(aes(x = doy, y = year))


dates <- unique(cr_high$acquisitionDate)
p <- '/media/steppe/hdd/Geospatial_data/Sentinel2EriogonumColoradense/'         

# we can only download 2500x2500 of cells per time. Use tiles to accomplish this. 
fname <- paste0(p, 'rawTiles', '_.tif')
# template <- setValues(template, 1)
# template <- disagg(template, fact = 1000)
# template <- disagg(template, fact = 4)

setwd(paste0(p, 'template'))
terra::makeTiles(template, c(2500, 2500), filename = ".tif")

setwd(paste0(p, 'NDSI_raw'))
f <- paste0('../template/', list.files('../template/'))

aoImporter <- function(x){
  
  aoi <- sf::st_as_sfc(sf::st_bbox(terra::rast(x)))
  tileNO <- gsub('.tif', '', basename(x))
  for (i in seq_along(dates)){
    
    if(!file.exists(paste0(tileNO, '_', dates[i], '.tif'))) {
      message('Downloading: ', tileNO, '_', dates[i])
    
    CDSE::GetImageByAOI(
      aoi = aoi, 
      time_range = dates[i], 
      script = script_file,
      mask = TRUE, 
      collection = "sentinel-2-l2a", 
      format = "image/tiff",
      file = paste0(tileNO, '_', dates[i], '.tif'),
      mosaicking_order = "leastCC", # LEAST CLOUD COVER. 
      resolution = 10, 
      client = OAuthClient
    )
    }
  }
}

# lapply(f, aoImporter)

#### Now determine the last date snow was found in each cell each YEAR. #####
p <- '/media/steppe/hdd/Geospatial_data/Sentinel2EriogonumColoradense/'   
setwd(paste0(p, 'NDSI_raw'))

r <- terra::rast('22_2023-07-05.tif')
plot(r)
names(r) <- c('B02', 'B03', 'B04', 'B11')

r2 <- (r$B03 - r$B11) / (r$B03 + r$B11)
msk <- terra::ifel(r2 < 0.42, NA, 1)
r2 <- terra::mask(r2, msk)
plot(r2)

# group by TILE NUMBER, THEN YEAR, THEN DATE. 
# Calculate NDSI for each CELL X TIME INTERVAL 

FILES  <- list.files()
tiles <- data.frame(
  FNAME = FILES, 
  TILE = as.numeric(gsub('_.*', '', FILES)),
  YEAR = as.numeric(gsub('_|-', '', stringr::str_extract(FILES, '_[0-9]{4}-'))),
  DOY = lubridate::yday(
    as.Date(
      gsub('_|[.]', '', stringr::str_extract(FILES, '_.*[.]')), 
      format = "%Y-%m-%d"
      )
  )
) |>
  dplyr::arrange(TILE, YEAR, DOY)

tiles <- split(tiles, f = ~ TILE + YEAR)
tiles['22.2023']


# first we will create a cloud mask - it won't work perfectly but should function OK 
# for our purposes. 

# these are the terms for the braaten-cohen-yang cloud detection 
# index https://custom-scripts.sentinel-hub.com/custom-scripts/sentinel-2/cby_cloud_detection/ 

cloud_detection <- function(x){

  NDGR <- vector(mode = 'list', length = nrow(x))
  bRatio <- vector(mode = 'list', length = nrow(x))
  cloud <- vector(mode = 'list', length = nrow(x))
  
  for (i in 1:nrow(x)){
    r <- terra::rast(x$FNAME[i]) # read in each file 
    names(r) <- c('B02', 'B03', 'B04', 'B11') # name layers
    NDGR[[i]] <- (r$B03 - r$B04) / (r$B03 + r$B04)
    bRatio[[i]] <- ((r$B03 - 0.175) / (0.39 + 0.175))
    cloud[[i]] <- terra::ifel(r$B11 > 0.1 & bRatio[[i]] > 0 & NDGR[[i]] > 0, 1, 0)
  }
  stack <- terra::rast(cloud)
  names(stack) <- x$DOY
  return(stack)
}

ob <- lapply(tiles['22.2023'], cloud_detection)
plot(ob[['22.2023']])
writeRaster(ob[['22.2023']], '../test_tile_22.tif', overwrite = TRUE)





NDSI_summary <- function(x){
  
  stack <- vector(mode = 'list', length = nrow(x))

  for (i in 1:nrow(x)){
    r <- terra::rast(x$FNAME[i]) # read in each file 
    names(r) <- c('B02', 'B03', 'B04', 'B11') # name layers
    stack[[i]] <- (r$B03 - r$B11) / (r$B03 + r$B11) # calculate NDSI 
    msk <- terra::ifel(stack[[i]] < 0.42, NA, 1) # mask not snow pixels 
    stack[[i]] <- terra::mask(stack[[i]], msk) # 
    stack[[i]] <- terra::ifel(!is.na(stack[[i]]), 1, NA)
  }
  stack <- terra::rast(stack)
  names(stack) <- x$DOY
  
  return(stack)
}

ob <- lapply(tiles['22.2023'], NDSI_summary)
plot(ob[['22.2023']])

writeRaster(ob[['22.2023']], '../test_tile_22.tif', overwrite = TRUE)

################### NOW SELECT BANDS FOR NDVI CALCULATION #####################

# We want to calculate this during ~ July 10 - July 22
# DOY 192-204
# It will need to be in our high moisture years: 2017, 2019, 2023
# let's aim for the data set within those years with the absolute lowest amount of
# cloud cover. 

cr_high <- cr |>
  filter(year %in% c(2017, 2019, 2023, 2024))

# those years were no good... However, 2024 was an OK-GOOD year for precip, and
# a flight on the 4th of July had no clouds! (trust me.. i was there ;-))

filter(cr_high) %>%  
  ggplot() + 
  geom_point(aes(x = doy, y = year, size = meanCloud))

independence_day <- filter(cr_high, acquisitionDate == '2024-07-04')
dates <- unique(independence_day$acquisitionDate)

setwd(paste0(p, 'NDVI_raw'))
f <- paste0('../template/', list.files('../template/'))
script_file <- '/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts/NDVI_download.js'

aoImporter <- function(x){
  
  aoi <- sf::st_as_sfc(sf::st_bbox(terra::rast(x)))
  tileNO <- gsub('.tif', '', basename(x))
  for (i in seq_along(dates)){
    
    if(!file.exists(paste0(tileNO, '_', dates[i], '.tif'))) {
      message('Downloading: ', tileNO, '_', dates[i])
      
      CDSE::GetImageByAOI(
        aoi = aoi, 
        time_range = dates[i], 
        script = script_file,
        mask = TRUE, 
        collection = "sentinel-2-l2a", 
        format = "image/tiff",
        file = paste0(tileNO, '_', dates[i], '.tif'),
        mosaicking_order = "leastCC", # LEAST CLOUD COVER. 
        resolution = 10, 
        client = OAuthClient
      )
    }
  }
}

lapply(f, aoImporter)


########## Now calculate NDVI for each cell #########
setwd(paste0(p, 'NDVI_raw'))
f <- list.files()

ndvi_calc <- function(x){
  r <- terra::rast(x)
  names(r) <- c('Red', 'NIR')
  r2 <- (r$NIR - r$Red) / (r$NIR + r$Red)
  terra::writeRaster(r2, paste0('../NDVI_Calc/', x))
}

parallel::mclapply(f, ndvi_calc) # only like 1.1 gb, lets
# add them all together into a single template, and resample from there. 
setwd(paste0(p, 'NDVI_Calc'))
f <- list.files()

lapply(f, function(x){res(rast(x))})
res(rast(f[5]))
calced <- terra::vrt(f, '../NDVI_Products/NDVI_vrt')
r <- rast("../NDVI_Products/NDVI_vrt")
terra::writeRaster(r, "../NDVI_Products/NDVI.tif")
r <- terra::rast("../NDVI_Products/NDVI.tif")
names(r) <- 'NDVI'

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed')
terra::resample(r, arc3_template, threads = 16, filename = './dem_3arc/NDVI.tif') 
terra::resample(r, arc1_template, threads = 16, filename = './dem_1arc/NDVI.tif') 
terra::resample(r, arc13_template, threads = 16, filename = './dem_1-3arc/NDVI.tif')
terra::resample(r, m3_template, threads = 16, filename = './dem_3m/NDVI.tif')

# we can also calculate SAVI real quick
setwd(paste0(p, 'NDVI_raw'))
f <- list.files()

lapply(f, function(x){res(rast(x))})


savi_calc <- function(x){
  r <- terra::rast(x)
  names(r) <- c('Red', 'NIR')
  savi <- ((r['NIR'] - r['Red']) / (r['NIR'] + r['Red'] + 0.5)) * (1.5)
  terra::writeRaster(savi, paste0('../SAVI_Calc/', x))
}

parallel::mclapply(f, savi_calc)

setwd(paste0(p, 'SAVI_Calc'))
f <- list.files()

lapply(f, function(x){res(rast(x))})
res(rast(f[5]))
calced <- terra::vrt(f, '../SAVI_Products/SAVI_vrt')
r <- rast("../SAVI_Products/SAVI_vrt")
terra::writeRaster(r, "../SAVI_Products/SAVI.tif")
r <- terra::rast("../SAVI_Products/SAVI.tif")
names(r) <- 'SAVI'

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed')
terra::resample(r, arc3_template, threads = 16, filename = './dem_3arc/SAVI.tif') 
terra::resample(r, arc1_template, threads = 16, filename = './dem_1arc/SAVI.tif') 
terra::resample(r, arc13_template, threads = 16, filename = './dem_1-3arc/SAVI.tif')
terra::resample(r, m3_template, threads = 16, filename = './dem_3m/SAVI.tif')




## I want to make a few more geomorphology variables which were'nt made with the first 
# batch 
morphoMaker <- function(x, p){
  
  if(str_detect(x, 'tif$')){demIN <- x} else{
    demIN <- terra::mosaic(terra::sprc(file.path(x, list.files(x, 'tif$'))))
  }
  # these describe different dimensions of curavture at a cell relative to it's local
  # landscape 
  whitebox::wbt_minimal_curvature(demIN, output = file.path(p, "mininmal_curv.tif"))
  whitebox::wbt_mean_curvature(demIN,  output = file.path(p, "mean_curv.tif"))
  
}

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')

arc3 <- rast('../data/spatial/processed/dem_3arc/dem.tif')
morphoMaker(x = '../data/spatial/processed/dem_3arc/dem.tif', 
            p = '../data/spatial/processed/dem_3arc/geomorphology')

morphoMaker(x = '../data/spatial/processed/dem_1arc/dem.tif', 
            p = '../data/spatial/processed/dem_1arc/geomorphology')

morphoMaker(x = '../data/spatial/processed/dem_1-3arc/dem.tif', 
            p = '../data/spatial/processed/dem_1-3arc/geomorphology')

morphoMaker(x = '../data/spatial/raw/dem_3m', 
            p = '../data/spatial/processed/dem_3m/geomorphology')


# we will also decompose aspect into Northness and Eastness

x = c(
  './dem_3arc/geomorphology/aspect.tif', 
  './dem_1arc/geomorphology/aspect.tif', 
  './dem_1-3arc/geomorphology/aspect.tif'
)

decomposeAspect <- function(x){
  
  Aspect <- terra::rast(x)
  
  northness <- cos(Aspect * pi / 180)
  eastness <- sin(Aspect * pi / 180)
  names(northness) <- 'Northness'
  names(eastness) <- 'Eastness'
  
  file.path(dirname(x), 'northness.tif')
  writeRaster(northness, file.path(dirname(x), 'northness.tif'))
  writeRaster(eastness, file.path(dirname(x), 'eastness.tif'))
}

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed')
lapply(x, decomposeAspect)


#### we will also create a couple of products from spatialEco package which relate 
# to how much heat areas on a lanDscape accumulate
x = c(
  './dem_3arc/geomorphology/aspect.tif', 
  './dem_1arc/geomorphology/aspect.tif', 
  './dem_1-3arc/geomorphology/aspect.tif'
  )

heat <- function(x){
  
  r <- terra::rast(x)
  dahi_r <-  spatialEco::dahi(r)
  names(dahi_r) <- 'DAHI'
  terra::writeRaster(dahi_r, file.path(dirname(x), 'DAHI.tif'))
  
  hli <- spatialEco::hli(r)
  names(hli) <- 'HLI'
  terra::writeRaster(r, file.path(dirname(x), 'HLI.tif'))
  
}

setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed')
lapply(x, heat)






#############################################################################

# Convert geomorphon and pennock landform positions into a factor rasters ##

p <- '/media/steppe/hdd/EriogonumColoradenseTaxonomy/data/spatial/processed'
f <- list.files(p, recursive=TRUE, pattern = '.tif$')
f <- file.path(p, f[grep('Pennock|geomorphons',  f)])

geomorphon_df <- data.frame(
  Value	= 1:10, 
  Landform = c(
    'Flat', 'Summit', 'Ridge', 'Shoulder', 'Spur', 'Slope', 'Hollow', 
    'Footslope', 'Valley', 'Pit')
)

pennock_df <- data.frame(
  Value = c(1:7, 128),
  Class = c('CFS', 'DFS', 'CSH', 'DSH', 'CBS', 'DBS', 'L', NA)
)

makeFactors <- function(x){
  
  f_out <- file.path(
    dirname(x), 
    gsub('.tif', '-fact.tif',  basename(x))
  )
  
  r <- terra::rast(x)
  
  if(grepl('Pennock', basename(x)) == TRUE){
    levels(r) <- pennock_df
  } else {
    levels(r) <- geomorphon_df
  }

  terra::writeRaster(r, f_out)
  unlink(x)
}

lapply(f, makeFactors)
