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
  # RandomForests do an excellent job of this - likely better than Boruta analysis... 
  # However what ends up happening is that many vars with very minor contributions to the
  # model are kept. WHen it comes time to predict these models onto gridded surfaces
  # these added terms (oftentimes requiring the re-reading of the rasters in a virtual context)
  # make the prediction take AGES. 
  # We want to avoid everything taking ages. 
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
      if(resolution == '3m'){ntile = 2} # 4 tiles. 
  
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


#' assign each random point to it's nearest trail
#' @param trails same as input to parent function `OrderRandomPtsIDS`
#' @param points same as input to parent function `OrderRandomPtsIDS`
point_assigner <- function(trails, points){
  segments <- sf::st_segmentize(trails, dfMaxLength = 150)
  vertices <- sf::st_cast(segments, 'POINT')
  
  segments <- lwgeom::st_split(segments, vertices) |>
    sf::st_collection_extract("LINESTRING") |>
    unique() |>
    dplyr::mutate(lineID = 1:n(), .before = geometry) |>
    dplyr::mutate(lineID = dplyr::case_when(
      lineID == min(lineID) ~ max(lineID), 
      lineID != min(lineID) ~ lineID - 1, 
      .default = as.numeric(lineID)
    )) |>
    arrange(lineID)
  
  points_cp <- sf::st_transform(points, sf::st_crs(segments))
  points_cp <- points_cp[ # determine which points are within the buffered areas  
    lengths(
      sf::st_intersects(points_cp, sf::st_buffer(segments, dist = 45)
      )
    ) > 0,
  ]
  points_cp <- points_cp |> 
    dplyr::mutate(
      Segment = st_nearest_feature(points_cp, segments),
      Trail = trails$Trail[1], 
      UID = 1:dplyr::n(),
      Pts = dplyr::n(), .before = geometry) |>
    dplyr::arrange(Segment) |>
    dplyr::mutate(TrailID = 1:n(), .before = geometry) |>
    dplyr::select(-Segment)
  
  return(points_cp)
}

#' snap occurrence data to trails which will be surveyed
#' @param trails same as input to parent function `OrderRandomPtsIDS`
#' @param presences same as input to parent function `OrderRandomPtsIDS`
#' @param final_pts the output of `point_assigner` within the fn `OrderRandomPtsIDS` 
#' after some minor cleaning
presenceSnapper <- function(trails, presences, final_pts){
  
  nearest <- sf::st_nearest_feature(presences, trails)
  Trail <- trails[nearest,] |>
    dplyr::pull(Trail)
  presences <- dplyr::mutate(presences, Trail = Trail, .before = geometry)
  
  Dist <- as.numeric(sf::st_distance(presences, trails[nearest,], by_element = T))  
  
  pres <- dplyr::mutate(presences, Dist = Dist,
                        Trail = dplyr::if_else(Dist <= 400, Trail, NA), .before = geometry)
  # now associate the presences points with the nearest random point along the 
  # trail and give it a letter associated with that plots TrailID. 
  
  revisit <- dplyr::filter(pres, ! is.na(pres$Trail))
  novisit <- dplyr::filter(pres, is.na(pres$Trail))
  revis_trails <- split(revisit, f = revisit$Trail)
  
  addTrailID <- function(revis_trails, random_pts){
    
    rnd_pts <- dplyr::filter(random_pts, Trail %in% 
                               unique(dplyr::pull(revis_trails, Trail)))
    TrailIDs <- sf::st_drop_geometry(rnd_pts)[sf::st_nearest_feature(revis_trails, rnd_pts), 'TrailID'] 
    
    revis_trails <- dplyr::mutate(revis_trails, TrailID = TrailIDs, .before = geometry)
    revis_trails <- dplyr::group_by(revis_trails, TrailID) |>
      mutate(TrailID = paste0(TrailID, LETTERS[1:n()])) |>
      dplyr::select(-Dist)
    
    return(revis_trails)
  }
  
  revis_trails <- dplyr::bind_rows(
    lapply(revis_trails, FUN = addTrailID, random_pts = final_pts)
  )
  return(revis_trails)
}


#' function to put a decent order to the ID's of random points based on trail positions
#' @param trails the trails network used to draw the random sample
#' @param points the random output points for which novel survey data should be collected
#' @param revisits the existing occurrence points which can be revisted as novel data are collected
OrderRandomPtsIDS <- function(trails, points, revisits){
  
  # first we split the data frame and lapply the process through each constituent
  
  split_trails <- split(trails, f = trails$Trail)
  raw_pts <- lapply(split_trails, FUN = point_assigner, points = points) |> 
    dplyr::bind_rows()
  
  # now ensure that any plots which are 'duplicated' across trails are reduced 
  # to a single representative. The plot will be retained in the trail which has
  # fewer points. 
  
  duplicate_pts <- raw_pts[lengths(sf::st_intersects(raw_pts)) > 1, ]
  singletons <- duplicate_pts |> 
    dplyr::group_by(geometry) |> 
    dplyr::slice_min(n = 1, order_by = Pts) 
  duplicate_pts <- duplicate_pts[! duplicate_pts$UID %in% singletons$UID, ]
  # now remove all but one copy of the duplicate points. 
  
  final_pts <- raw_pts[! raw_pts$UID %in% duplicate_pts$UID, 
                       -which(names(raw_pts) %in% c("Pts"))]
  final_pts$UID <- 1:nrow(final_pts)
  
  
  # now assign the historic occurrence points to the relevant trail #
  # since these were not generated off of our buffer trails, we will #
  # restrict the search to points within 400m of our trails, the other #
  # points we would need permits to access. #
  occurrence_revisits <- presenceSnapper(trails, occurrence_only, final_pts)
  occurrence_revisits <- dplyr::mutate(occurrence_revisits, 
                                       UID = paste0('P-', 1:n()), .before = geometry) 
  return(list(final_pts, occurrence_revisits))
}



#' Split text across multiple lines for use in a text grob
#' @params x a character vector which should be split out
#' @parms width the number of characters which each line should contain 
newliner <- function(x, width){
  
  # determine how many lines the text will be split across. 
  lines <- floor(nchar(x)/ width); wc_line <- ceiling(nchar(x)/lines)
  theo_cuts <- c(1, (wc_line * seq(lines))[1:lines-1])
  
  # determine sensible cuts for the words. Essentially just use the existing spaces.
  emp_cuts <- vector(mode = 'double', length = lines)
  for (i in seq(theo_cuts)){ # determine sensible cuts for the words. 
    emp_cuts[i] <- which(charToRaw(x) == '20')[
      which.min(abs(theo_cuts[i] - which(charToRaw(x) == '20')))]
  } 
  emp_cuts[1] <- 1
  end_cuts <- c(emp_cuts[2:length(emp_cuts)] - 1, nchar(x))
  
  # now we can just iterate through the string with a loop adding the newline
  # operator as required.
  new <- vector(mode = 'list', length = lines)
  for (i in seq(emp_cuts)){
    new[i] <- stringr::str_trim(rawToChar(charToRaw(x)[emp_cuts[i]:end_cuts[i]]))
  }
  newlined <- paste(new, collapse = '\n')
  
  # new <- list(emp_cuts, end_cuts)
  return(newlined)
  
}

