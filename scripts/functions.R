rastReader <- function(x, p2proc){

  # create paths for the first level of files. 
  f <- file.path(p2proc, x, list.files(file.path(p2proc, x), pattern = '.tif$'))

  # create paths for the geomorphology files.
  geo_f <- file.path(p2proc, x, 'geomorphology',
                     list.files(file.path(p2proc, x, 'geomorphology')))
  
  geo_f <- geo_f[! grepl('D8pntr|basins|DEM', geo_f) ] # remove the D8pntr & Basins 
  
  # combine the file paths
  f <- c(f, geo_f)
  s <- terra::rast(f)
  
}
