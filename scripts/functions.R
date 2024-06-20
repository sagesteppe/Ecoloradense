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

#' identify the SW corner of a 3m raster cell for ground truthing. 
#' 
#' @param template
SWcorner <- function(x){
  
  
}



#' @param x input occurrence data
#' @param resolution list of paths to geodata at different resolutions. 
#' @param iteration numeric, which iteration of modelling is being performed? 
modeller <- function(x, resolution, iteration){
  
  rast_dat <- rastReader(paste0('dem_', resolution), p2proc) 
  
  df <- dplyr::bind_cols(
    dplyr::select(x, Occurrence), 
    dplyr::select(terra::extract(rast_dat, x), -ID), 
  ) |> 
    tidyr::drop_na() %>% 
    dplyr::mutate(
      Occurrence = as.factor(Occurrence), 
      Pennock = as.factor(Pennock), 
      geomorphons = as.factor(geomorphons), 
      Longitude = unlist(purrr::map(.$geometry,1)),
      Latitude  = unlist(purrr::map(.$geometry,2))) |> 
    sf::st_drop_geometry()
  
  # first we will perform boruta analysis, this will drop variables which have
  # no relationship to the marks at the resolution under analysis. 
  cores <- parallel::detectCores()
  BorutaRes <- Boruta::Boruta(Occurrence ~ ., data = df, num.threads = cores, doTrace = 0)
  importance <- Boruta::attStats(BorutaRes)
  rn <- rownames(importance[importance$decision %in% c('Confirmed'),])
  important_vars <- Boruta::getSelectedAttributes(BorutaRes, withTentative = F)
  
  df <- dplyr::select(df, dplyr::all_of(c('Occurrence', important_vars)))
  rm(BorutaRes, importance, rn, important_vars)
  
  # define the cross-fold validation 
  repeat_cv <- caret::trainControl(method = 'repeatedcv', number = 5, repeats = 10) 
  
  # split the input data into both an explicit train and test data set. 
  TrainIndex <- caret::createDataPartition(
    df$Occurrence, p = .8, list = FALSE, times = 1)
  Train <- df[ TrainIndex,]; Test <- df[-TrainIndex,]
  
  # perform the random forest modelling
  rf_model <- ranger::ranger(
    Occurrence ~ ., data = Train, probability = T, keep.inbag = TRUE)
  
  # save the model
  saveRDS(rf_model,
          file = paste0('../results/models/', resolution, '-Iteration', iteration, '.rds'))
  
  # save the confusion matrix.
  predictions <- predict(rf_model, Test, type = 'se', se.method = 'infjack', probability=TRUE)
  predictions$binary <- as.factor(if_else(predictions$predictions[,2] <= 0.49, 0, 1))
  cmRestrat <- caret::confusionMatrix(predictions$binary, Test$Occurrence)
  saveRDS(cmRestrat,
          file = paste0('../results/tables/', resolution, '-Iteration', iteration, '.rds'))
  
  # create these layers of coordinates for the models. 
  Longitude <- init(rast_dat, 'x') ; names(Longitude) <- 'Longitude'
  Latitude <- init(rast_dat, 'y') ; names(Latitude) <- 'Latitude'
  rast_dat <- c(rast_dat, Longitude, Latitude)
  
  pout <- '../results/suitability_maps'
  
  
  ###################      PREDICT ONTO SURFACE        #########################
  pr <- function(...) predict(..., type = 'response', num.threads = 1)$predictions 
  
  dir.create( file.path(pout, 'tiles'), showWarnings = F)
  dir.create( file.path(pout, 'pr_tiles'), showWarnings = F)
  if(resolution == '3arc'){ntile = 1}
  if(resolution == '1arc'){ntile = 1} else # has succeeded -narrowly- at 1 
    if(resolution == '1-3arc'){ntile = 2} else # 4 tiles 
      if(resolution == '3m'){ntile = 2} # 9 tiles. 
  
  if(ntile > 1){
    template <- rast(nrows = ntile, ncols = ntile, extent = ext(rast_dat), crs = crs(rast_dat))
    message('Making tiles for prediction')
    terra::makeTiles(rast_dat, template, filename = file.path(pout, 'tiles', "tile_.tif"))
    tiles <- file.path(pout, 'tiles', list.files(file.path(pout, 'tiles')))
    
    pb <- txtProgressBar(
      min = 0,  max = length(tiles), style = 3, width = 50, char = "+")
    for (i in seq_along(tiles)){
      terra::predict(
        rast(tiles[i]), rf_model, f = pr, 
        filename = file.path(pout, 'pr_tiles', paste0(i, '.tif')),
        overwrite = T, na.rm=TRUE)
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    mos <- terra::mosaic(sprc(
      file.path(pout, 'pr_tiles', list.files(file.path(pout, 'pr_tiles')))
    ), fun = 'mean')
    writeRaster(
      mos[[2]], 
      wopt = c(names = 'Probability'),
      filename = file.path(pout, paste0(resolution, '-Iteration', iteration, '-Pr.tif')))
    
    unlink(file.path(pout, 'tiles')); unlink(file.path(pout, 'pr_tiles'))
  } else { # don't want to spend time copying the contents of the predictor
    terra::predict( # rasters, skip straight to modelling. 
      type = 'prob', cores = 1, f = pr,
      rast_dat, rf_model, cpkgs = "ranger",
      filename = file.path(pout, paste0(resolution, '-IterationCoarse', iteration, '.tif')),
      wopt = c(names = 'predicted_suitability'), 
      overwrite = T)
    
    writeRaster(
      rast(file.path(pout, paste0(resolution, '-IterationCoarse', iteration, '.tif')))[[2]],
      file.path(pout, paste0(resolution, '-Iteration', iteration, '-Pr.tif'))
    )
    file.remove(file.path(pout, paste0(resolution, '-IterationCoarse', iteration, '.tif')))
  }
  
  unlink(file.path(pout, 'tiles')); unlink(file.path(pout, 'pr_tiles'))
  gc(verbose = FALSE)
  
  
  ###############      STANDARD ERROR PREDICTION        ########################
  # the SE predictions are absurdly memory hungry, We will create 'tiles' to 
  # predict onto, and then combine them once the predictions are complete
  
  dir.create( file.path(pout, 'tiles'), showWarnings = F)
  dir.create( file.path(pout, 'se_tiles'), showWarnings = F)
  if(resolution == '3arc'){ntile = 2} else # 4 tiles 
    if(resolution == '1arc'){ntile = 5} else # 25 tiles needed # 16 FAILED 
      if(resolution == '1-3arc'){ntile = 8} else # 81 tiles
        if(resolution == '3m'){ntile = 12} # 144 tiles. 
  
  template <- rast(nrows = ntile, ncols = ntile, extent = ext(rast_dat), crs = crs(rast_dat))
  message('Making tiles for prediction of standard errors')
  terra::makeTiles(rast_dat, template, filename = file.path(pout, 'tiles', "tile_.tif"))
  tiles <- file.path(pout, 'tiles', list.files(file.path(pout, 'tiles')))
  
  se <- function(...) predict(..., type = 'se', se.method = 'infjack',
                              predict.all = FALSE,
                              probability = TRUE, num.threads = cores)$se
  
  # # predict the standard error onto the tiles.
  # create and initialize progress bar
  pb <- txtProgressBar(
    min = 0,  max = length(tiles), style = 3, width = 50, char = "+")
  for (i in seq_along(tiles)){
    terra::predict(
      rast(tiles[i]), rf_model, f = se, 
      filename = file.path(pout, 'se_tiles', paste0('SE', i, '.tif')),
      overwrite = T, na.rm=TRUE)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  mos <- terra::mosaic(sprc(
    file.path(pout, 'se_tiles', list.files(file.path(pout, 'se_tiles')))
  ), fun = 'mean')
  writeRaster(
    mos[[2]], 
    wopt = c(names = 'standardError'),
    filename = file.path(pout, paste0(resolution, '-Iteration', iteration, '-SE.tif')))
  
  unlink(file.path(pout, 'tiles')); unlink(file.path(pout, 'se_tiles'))
  gc(verbose = FALSE)
  
}
