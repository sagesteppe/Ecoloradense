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


ensure_multipolygons <- function(X) { # @ stackoverflow
  tmp1 <- tempfile(fileext = ".gpkg")
  tmp2 <- tempfile(fileext = ".gpkg")
  st_write(X, tmp1)
  gdalUtilities::ogr2ogr(tmp1, tmp2, f = "GPKG", nlt = "MULTIPOLYGON")
  Y <- st_read(tmp2)
  st_sf(st_drop_geometry(X), geom = st_geometry(Y))
}


#' @param x input occurrence data
#' @param resolution list of paths to geodata at different resolutions. 
#' @param iteration numeric, which iteration of modelling is being performed? 
#' @param se_prediction boolean, whether to predict the SE surfaces or not, can add roughly
#' a week onto the prediction at 3m. 
#' @param p2proc path to the processed raster data. 
#' @param train_split Numeric. The proportion of data to use for train, at first iteration
#' a standard 0.8, good, but lot's of data are available for test later so up to 0.9 good. 
modeller <- function(x, resolution, iteration, se_prediction, p2proc, train_split){
  
  if(missing(se_prediction)){se_prediction <- FALSE}
  rast_dat <- rastReader(paste0('dem_', resolution), p2proc) 
  cores <- parallel::detectCores()
  
  df <- dplyr::bind_cols(
    dplyr::select(x, Occurrence), 
    dplyr::select(terra::extract(rast_dat, x), -ID), 
  ) |> 
    tidyr::drop_na() %>% 
    dplyr::mutate(
      Occurrence = as.factor(Occurrence)) |>#, 
 #     Pennock = as.factor(Pennock), 
  #    geomorphons = as.factor(geomorphons)) |> 
    sf::st_drop_geometry()

  #######################         MODELLING           ##########################
  # this is the fastest portion of this process. It will take at most a few minutes
  # at any resolution, the goal of this paper wasn't really focused on comparing
  # multiple models of families nor messing with parameters, so we use Random Forest
  # simply because they work very well with default settings, almost 'controlling' 
  # for stochasticity in this portion of the work. 
  
  # only rerun modelling if a saved .rds object doesn't exist. 
  if(file.exists(paste0('../results/models/', resolution, '-Iteration', iteration, '.rds'))){
    rf_model <- readRDS(paste0('../results/models/', resolution, '-Iteration', iteration, '.rds'))
    message(
      'An exising model for this resolution and iteration already exists; reloading it now for projection')
  } else {
  
  # split the input data into both an explicit train and test data set. 
  TrainIndex <- caret::createDataPartition(
    df$Occurrence, p = train_split, list = FALSE, times = 1)
  Train <- df[TrainIndex,]; Test <- df[-TrainIndex,]

  Train.sf <- sf::st_as_sf(Train, coords = c('Longitude', 'Latitude'), crs = 32613)

  # perform the random forest modelling using default settings. 
  indices_knndm <- CAST::knndm(Train.sf, rast_dat, k=10)
  
  tgrid <- expand.grid(
    mtry = 
      floor(ncol(Train)/ 2.4):floor(ncol(Train)/ 1.6),
    splitrule = 'gini',
    min.node.size = c(1:9, seq(10, 30, by = 5))
  )

  # calculate class weights for the factor levels. 
  wts = c(
    1 - sum(Train$Occurrence==0)/nrow(Train),
    1 - sum(Train$Occurrence==1)/nrow(Train)
  )
  
  message('Tuning hyperparameters.')
  model <- caret::train(
    x = Train[,-grep('Occurrence', colnames(Train))],
    y = Train$Occurrence,
    method = "ranger",
    num.trees = 750,
    tuneGrid = tgrid,
    metric = 'Accuracy', 
    keep.inbag = TRUE, 
    class.weights = wts,
    importance = 'permutation',
    trControl = caret::trainControl(
      method="cv",
      index = indices_knndm$indx_train,
      savePredictions = "final", 
      allowParallel = TRUE)
  )
  
  rf_model <- ranger::ranger(
    Occurrence ~ ., data = Train, probability = T, keep.inbag = TRUE, 
    mtry = model[['finalModel']][['mtry']], 
    min.node.size = model[['finalModel']][['min.node.size']],
    importance = 'permutation',
    class.weights = wts
    )
  
  # save the model
  saveRDS(rf_model,
          file = paste0('../results/models/', resolution, '-Iteration', iteration, '.rds'))
  
  # save the confusion matrix
  predictions <- predict(rf_model, Test, type = 'se', se.method = 'infjack', probability=TRUE)
  predictions$binary <- as.factor(if_else(predictions$predictions[,2] <= 0.49, 0, 1))
  
  cmRestrat <- caret::confusionMatrix(predictions$binary, Test$Occurrence, 
                                      positive = '1', mode = 'everything')
  saveRDS(cmRestrat,
          file = paste0('../results/tables/', resolution, '-Iteration', iteration, '.rds'))
  
  # we will also save the Pr-AUC as some folks find it a useful metric for class unbalanced
  # data sets. 
  df_prauc <- data.frame(
    truth = Test$Occurrence,
    Class1 = predictions$predictions[,1]
  ) 
  yardstick::pr_auc(df_prauc, truth, Class1) |>
    mutate(resolution = resolution, iteration = iteration) |>
    write.csv(
      paste0('../results/tables/', resolution, '-Iteration', iteration, '.csv'),
      row.names = F)
  
  rm(df, df_prauc, TrainIndex, Train, Test, predictions, cmRestrat)
  }
  
  pout <- '../results/suitability_maps'
  ###################      PREDICT ONTO SURFACE        #########################
  # this prediction is generally straightforward, it doesn't take an obscene
  # amount of RAM and can happen relatively quickly, just overnight for the 
  # higher resolution scenarios. 
  
    pr <- function(...) predict(..., type = 'response', num.threads = 1)$predictions 
  
    if(resolution %in% c('3arc', '1arc', '1-3arc', '3m')){ntile = 1} else
      if(resolution == '1m'){ntile = 4}  
  
    if(ntile > 1){
    
     tile_path_4pr <- file.path(pout, paste0('tilesPR', '_', resolution))
     if(file.exists(tile_path_4pr)){
       message('Tiles already created for this resolution, skipping to prediction!\n')
     } else {
  
        dir.create(tile_path_4pr, showWarnings = F); message('Making tiles for prediction')
        template <- rast(nrows = ntile, ncols = ntile, extent = ext(rast_dat), crs = crs(rast_dat))
        terra::makeTiles(rast_dat, template, filename = file.path(tile_path_4pr, "tile_.tif"))
        rm(template)
      }
    
     tiles <- file.path(pout, 'tilesPR', list.files(file.path(pout, 'tilesPR')))
     pr_tile_path <- file.path(pout, paste0('pr_tiles', '_', resolution, iteration))
     dir.create(pr_tile_path, showWarnings = F)
    
     message('Writing Predicted Probability to Raster using tiles - this while take some time')
     pb <- txtProgressBar(
       min = 0,  max = length(tiles), style = 3, width = 50, char = "+")
     for (i in seq_along(tiles)){
       terra::predict(
         rast(tiles[i]), rf_model, f = pr, 
         filename = file.path(pr_tile_path, paste0(i, '.tif')),
          overwrite = T, na.rm=TRUE)
        setTxtProgressBar(pb, i)
        gc(verbose = FALSE)
     }
     close(pb)
    
     mos <- terra::mosaic(sprc(
       file.path(pr_tile_path, list.files(pr_tile_path))
     ), fun = 'mean')
      writeRaster(
        mos[[2]], 
        wopt = c(names = 'Probability'),
       filename = file.path(pout, paste0(resolution, '-Iteration', iteration, '-Pr.tif')))
    
    rm(mos)
    unlink(file.path(pout, 'pr_tiles')) # can remove the tiles we used for the probability surface.
    
    } else { # in these cases, we only need to use a single tile for prediction, we can
     # just use the existing virtual raster to do this. 
    
      if(!file.exists(file.path(pout, paste0(resolution, '-Iteration', iteration, '-Pr.tif')))){
      message('Writing Predicted Probability to Raster')
        terra::predict( # rasters, skip straight to modelling. 
          cores = 1, f = pr, 
          rast_dat, rf_model, cpkgs = "ranger",
          filename = file.path(pout, paste0(resolution, '-IterationCoarse', iteration, '.tif')),
          wopt = c(names = 'predicted_suitability'), na.rm=TRUE,
          overwrite = T)
    
        writeRaster( # we wrote a raster with both predicted class probabilities onto it. 
        # we don't want both, we are going to 'rewrite' the raster so only one class remains. 
          rast(file.path(pout, paste0(resolution, '-IterationCoarse', iteration, '.tif')))[[2]],
          file.path(pout, paste0(resolution, '-Iteration', iteration, '-Pr.tif')), 
          overwrite = T
      )
      file.remove(file.path(pout, paste0(resolution, '-IterationCoarse', iteration, '.tif')))
      } else {'Predicted Probability Raster already exists - skipping.'}
    }
  
  unlink(file.path(pout, 'pr_tiles'))
  gc(verbose = FALSE)
  
  ################ AREA OF APPLICABILITY SURFACE    ############################
  
  message('Writing Area of Applicability to Raster')
  AOAfun <- function(rf_model, r) {
    CAST::aoa(r, model, LPD = FALSE, verbose=FALSE)$AOA
  }
  
  terra::predict( # rasters, skip straight to modelling. 
    cores = 1, f = AOAfun, 
    rast_dat, rf_model, cpkgs = "CAST",
    filename = file.path(pout, paste0(resolution, '-IterationAOA', iteration, '.tif')),
    wopt = c(names = 'AOA'), na.rm=TRUE,
    overwrite = T)
  
  ###############      STANDARD ERROR PREDICTION        ########################
  # the SE predictions are absurdly memory hungry, We will create 'tiles' to 
  # predict onto, and then combine them once the predictions are complete
  # Note I use the term SE here, because it's what the R package uses, however, 
  # upon further reading what they are actually calculating is a Confidence 
  # interval - so great naming... https://github.com/imbs-hl/ranger/issues/136 
  
  if(se_prediction == TRUE){
    
  tile_path_4se <- file.path(pout, paste0(resolution, 'TilesSE'))
  final_se_path <- file.path(pout, paste0(resolution, 'SETiles'))
  
  # it seems the tiles pretty much need to be under 200MB or so for safe CI prediction???
  if(resolution == '3arc'){ntile = 2} else # 4 tiles 
    if(resolution == '1arc'){ntile = 6} else # 36 tiles Needed # earlier iteration with 16 FAILED # iter w/ 25 failed on tile 18...
      if(resolution == '1-3arc'){ntile = 16}  # 144 /121/81 failed - follow under 200mb rule. 
  
  if(resolution == '3m'){
    'Confidence Intervals cannot be produced for data of these size;
    it would require over 2k tiles and 1 week of compute'}  else {
  
    # don't create the tiles if they already exist. 
    if(exists(tile_path_4se)){
      message('Tiles for SE exist at this resolution, skipping to predict.')} else {
     
        message('Making tiles for prediction of standard errors')
        dir.create(tile_path_4se, showWarnings = F)
        dir.create(final_se_path, showWarnings = F)
      
        template <- rast(nrows = ntile, ncols = ntile, extent = ext(rast_dat), crs = crs(rast_dat))
        terra::makeTiles(rast_dat, template, filename = file.path(tile_path_4se, "tile_.tif"))
        rm(template)
    }
  
    tiles <- file.path(tile_path_4se, list.files(tile_path_4se))
    se <- function(...) predict(..., type = 'se', se.method = 'infjack',
                               predict.all = FALSE,
                                probability = TRUE, num.threads = cores/4)$se
  
    # # predict the standard error onto the tiles.
    # create and initialize progress bar
    message('Predicting SE to tiles; this may take a really long time')
    pb <- txtProgressBar(
      min = 0,  max = length(tiles), style = 3, width = 50, char = "+")
  
    for (i in seq_along(tiles)){
      if(!file.exists(file.path(final_se_path, paste0('SE', i, '.tif')))){ # in case of crash, pick up where left off. 
      terra::predict(
        rast(tiles[i]), rf_model, f = se, 
        filename = file.path(final_se_path, paste0('SE', i, '.tif')),
        overwrite = T, na.rm=TRUE)
      gc(verbose = FALSE)}
      setTxtProgressBar(pb, i)
    }
    close(pb)
  
    mos <- terra::mosaic(sprc(
      file.path(final_se_path, list.files(final_se_path))
    ), fun = 'mean')
    writeRaster(
      mos[[2]], 
      wopt = c(names = 'standardError'),
      filename = file.path(pout, paste0(resolution, '-Iteration', iteration, '-SE.tif')),
      overwrite = T)
    gc(verbose = FALSE)
    }
    unlink(final_se_path)
  }
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


#' Classify historic observations records in groups based on event date and geographic proximity
#' 
#' @description Simply measure the distance between observations that were made 
#' on the same day and bin points within 50m of one another. 
#' @param x a list of sf/tibble/dataframes with more than one record per day. 
HistObsGrps <- function(x){
  
  # calculate distances between all points
  dists <- st_distance(x)
  
  # and converts from a distances object to a simple numeric matrix. which works
  # the fns below. 
  dists <- matrix(as.numeric(dists), nrow = nrow(dists), byrow = TRUE)
  
  # tag points which we want connected - those within 50m of each other. 
  adj_matrix <- as.matrix(dists <= 50)
  graph <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected") 
  
  # we can extract each 'grp' (including unconnected individuals)
  grps <- split(igraph::V(graph), igraph::components(graph)$membership)
  grps <- lapply(grps, \(x) setNames(
    data.frame(as.numeric(x)),
    nm = 'Indv')
  ) |>
    data.table::rbindlist(idcol = 'Obs.Grp') 
  
  grps <- grps[order(grps$Indv),]
  dplyr::mutate(x, 
                Obs.grp = grps$Obs.Grp, .before = 'geometry'
  )
}

#' Reduce the number of points per training/test set so that no cells are duplicated. 
#' 
#' @description Each of the different resolutions of raster data will require
#' different testing and training data to ensure that records (e.g. a raster grid 
#' cell) are not replicated. This function first calculates the hypotenuse of a the resolution
#' using the 'ideal' X& Y edges of the raster - which we know to be inaccurate by 
#' ca. an entire meter at the '10m' resolution scale - to give us a rough estimate
#' of how many points will be dropped.  Then it loads an actual layer of the 
#' raster stack at the relevant resolution, extracts the pixel ID's and removes
#' the duplicate cells. At this stage it will look at the absences, and remove
#' any which are in the same cell as an occurrence. 
#' 
#' @param x sf data set of presences and absences. 
#' @param res Numeric. Simplified resolution as meters, one of 3, 10, 30, 90. 
#' @param root Path to root of spatial data. 
#' @param mode Character. One of 'Presence', 'Count', defaults to 'Presence'. 
#' If running in 'Presence' mode, then select one record per raster cell, always
#' choosing a presence over an absence. If running in mode 'Count' sum all records 
#' per raster cell, and relative the sample quadrat 3x3 m, to the size of the cell. 
#' For example, if using a 90m cell, and their are three quadrats: 
#' 
subset_pts <- function(x, res, root, mode){
  
  if(missing(mode)){mode <- 'Presence'}
  res_string <- switch(as.character(res),
                       "3" = "3m",
                       "10" = "1-3arc",
                       "30" = "1arc",
                       "90" = "3arc",
                       stop("Input to `res` invalid. Need be one of: 3, 10, 30, 90")
  )
  
  within_hypotenuse <- sum( 
    as.numeric(
      sf::st_distance(
        x, 
        x[sf::st_nearest_feature(x),],
        by_element = TRUE)
    ) < sqrt(2 * (res^2))
  ) # these records are possibly in the same cell as another
  # record for the finest resolution modelling.  
  # doesn't fit very well. anyways on to the empirical method. 
  # message('As many as: ', within_hypotenuse/2, ' records could be dropped')
  
  rasta <- terra::rast( # we only need to read in one raster - by definition these
    # are all aligned to each other. 
    file.path(root, paste0('dem_' ,res_string), 'dem.tif')
  )
  
  x['RasterCell'] <- terra::extract(rasta, x, cells = TRUE)$cell
  RC <- split(x, f = x$RasterCell) 
  
  if(mode=='Presence'){
    # if there are multiple points per cell, 1st) discard the absences (if applicable)
    # 2) randomly sample out one of the presences.  
    
    select_rec <- function(x){
      # if only presenc==0 | presenc==1, sample randomly. 
      # if both presenc==0 and presenc==1 co-exist, sample from presenc==1, 
      if(nrow(x)>1){
        if(all(x$Presenc==0) | all(x$Presenc==1) == TRUE){
          # remove the historic record if present. 
          if(length(unique(x$Type)==2)){
            x <- x[x$Type=='Current',]
          }
          
          x <- x[sample(1:nrow(x), 1),]
        } else {
          # remove the historic record if present. 
          if(length(unique(x$Type)==2)){
            x <- x[x$Type=='Current',]
          }
          
          xsub <- x[x$Presenc==1,]
          x <- xsub[sample(1:nrow(xsub), 1),]}
      }
      return(x)
    }
    RC <- lapply(RC, select_rec)
    
  } else { # if mode == 'Count'
    # area cell, area quadrant * mean plants per life stage. 
    
    countR <- function(x){
      x$Prsnc_J <- ((res^2) / (3^2)) * mean(x$Prsnc_J)
      x$Prsnc_M <- ((res^2) / (3^2)) * mean(x$Prsnc_M)
      x$Prsnc_S <- ((res^2) / (3^2)) * mean(x$Prsnc_S, na.rm = TRUE)
      x <- x[sample(1:nrow(x), 1),] 
    }
    RC <- lapply(RC, countR)
  }
  dplyr::bind_rows(RC) |>
    dplyr::arrange(OBJECT) |>
    dplyr::select(-RasterCell) |>
    sf::st_as_sf()
}
