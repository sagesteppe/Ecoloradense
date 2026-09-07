PROJ_ROOT <- '/media/steppe/hdd/EriogonumColoradenseTaxonomy'

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
#' @param PAratio Numeric. The ratio of presence to absences - simply for appending to file name, 
#' the exact quantity are saved in model objects. "1:2" 
#' @param distOrder Character. The distance between the nearest absence to presence, as a multiple of the cell resolution
#' e.g. distOrder 1 at a 90m resolution indicates the nearest absence is >90m away from a presence, at 10m it indicates > 10m away.
#' @param seed Numeric. Passed to \code{set.seed()} before the train/test split and
#' model fit, so replicate calls at the same PAratio/resolution can use independent
#' seeds instead of silently reusing whichever model was fit first.
#' @param predict_surface Boolean. Defaults to TRUE. When FALSE, skip the probability
#' surface, area-of-applicability, and standard-error raster predictions entirely -
#' fitting and holdout evaluation still happen, and are fast, so this is what a grid/
#' replicate search over PAratio should use; only the final, chosen model needs the
#' (slow) full-surface predictions.
#' @param results_root Character, a directory with models/modelsTune/tables/test_data/
#' evaluations/suitability_maps subfolders, written and read exactly like
#' \code{PROJ_ROOT}'s own \code{results/} tree. Defaults to that tree; override to run a
#' parallel comparison (e.g. a different train/test split strategy) without reading or
#' overwriting anything under the default results.
modeller <- function(x, resolution, iteration, se_prediction, p2proc, train_split,
                     PAratio, distOrder, seed = 1, predict_surface = TRUE,
                     split = NULL, results_root = file.path(PROJ_ROOT, 'results')){

  if(missing(se_prediction)){se_prediction <- FALSE}
  set.seed(seed)
  rast_dat <- rastReader(paste0('dem_', resolution), p2proc)
  cores <- parallel::detectCores()

  df <- dplyr::bind_cols(
    dplyr::select(x, Occurrence, dplyr::any_of('.id')),
    dplyr::select(terra::extract(rast_dat, x), -ID),
  ) |>
    tidyr::drop_na() %>%
    dplyr::mutate(
      Occurrence = as.factor(Occurrence)) |>
    sf::st_drop_geometry()

  fname <- paste0(resolution, '-Iteration', iteration, '-PA', PAratio, distOrder, '-Seed', seed)
  #######################         MODELLING           ##########################
  # this is the fastest portion of this process. It will take at most a few minutes
  # at any resolution, the goal of this paper wasn't really focused on comparing
  # multiple models of families nor messing with parameters, so we use Random Forest
  # simply because they work very well with default settings, almost 'controlling' 
  # for stochasticity in this portion of the work. 
  
  # only rerun modelling if a saved .rds object doesn't exist. Reload both objects: the
  # cross-validated caret `model` for CAST::aoa()'s weighting/CV folds - not just `rf_model`.
  if(file.exists(paste0(results_root, '/models/', fname, '.rds'))){
    rf_model <- readRDS(paste0(results_root, '/models/', fname, '.rds'))
    model <- readRDS(paste0(results_root, '/modelsTune/', fname, '.rds'))
    message(
      'An existing model for this resolution and iteration already exists;
      reloading it now for prediction')
  } else {
  
  # split the input data into both an explicit train and test data set. If a
  # pre-computed spatial split is supplied (see spatial_class_split()), every ratio/seed
  # fit in a resolution's grid search shares the same class-balanced, prediction-domain-
  # matched partition; otherwise fall back to a random stratified split.
  if(is.null(split)){
    TrainIndex <- caret::createDataPartition(
      df$Occurrence, p = train_split, list = FALSE, times = 1)
    Train <- df[TrainIndex,]; Test <- df[-TrainIndex,]
  } else {
    Train <- df[df$.id %in% split$train_id, ]
    Test  <- df[df$.id %in% split$test_id, ]
  }
  Train <- dplyr::select(Train, -dplyr::any_of('.id'))
  Test  <- dplyr::select(Test, -dplyr::any_of('.id'))

  Train.sf <- sf::st_as_sf(Train, coords = c('Longitude', 'Latitude'), crs = 32613)
  
  ## remove totally uninformative features using Boruta analysis, these features
  # can only hinder the model, and make prediction onto raster surfaces take longer
  
  BorutaRes <- Boruta::Boruta(Occurrence ~ ., data = Train, num.threads = cores, doTrace = 1)
  importance <- Boruta::attStats(BorutaRes)
  rn <- rownames(importance)
  important_vars <- Boruta::getSelectedAttributes(BorutaRes, withTentative = TRUE)
  
  keep <- unique(c('Occurrence', important_vars))
  Train <- Train[ , names(Train) %in% keep]
  
  # develop a cross validation structure which is explicitly spatial  
  indices_knndm <- CAST::knndm(Train.sf, rast_dat, k=10)
  
  # Now we will tune the hyperparameters for this model. 
  
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
  
  # save the cv fold fit model for AOA preds
  saveRDS(model, 
          file = paste0(results_root, '/modelsTune/', fname, '.rds'))
  # save the model
  saveRDS(rf_model,
          file = paste0(results_root, '/models/', fname, '.rds'))
  
  # save the test data. 
  write.csv(Test, paste0(results_root, '/test_data/', fname, '.csv'))
  
  
  # threshold the model and use these values for binary classification estimates
  predictions <- predict(rf_model, Test, type = 'se', se.method = 'infjack', probability=TRUE)
  e <- dismo::evaluate(
    p = predictions$predictions[Test$Occurrence==1,2],
    a = predictions$predictions[Test$Occurrence==0,2]
  )
  
  th <- dismo::threshold(e)
  saveRDS(th, file = paste0(results_root, '/evaluations/', fname, '-thresh.rds'))
  saveRDS(e, file = paste0(results_root, '/evaluations/', fname, '-eval.rds'))
  
  # save the confusion matrix
  predictions <- predict(rf_model, Test, type = 'se', se.method = 'infjack', probability=TRUE)
  predictions$binary <- as.factor(if_else(predictions$predictions[,2] <= th$spec_sens, 0, 1))
  
  cmRestrat <- caret::confusionMatrix(predictions$binary, Test$Occurrence, 
                                      positive = '1', mode = 'everything')
  saveRDS(cmRestrat,
          file = paste0(results_root, '/tables/', fname, '.rds'))
  
  # we will also save the pr-auc and ROC-auc metrics. 
  df_auc <- data.frame(
    truth = as.factor(Test$Occurrence),
    Class1 = predictions$predictions[,2]
  ) 
  
  pr_auc_val <- yardstick::pr_auc(
    data = df_auc, truth = truth, Class1, event_level = 'second') 
  roc_auc_val <- yardstick::roc_auc(
    df_auc, truth, Class1, event_level = 'second')
  
  setNames(
    data.frame( 
      rbind(
        pr_auc_val, 
        roc_auc_val
      )
    ), nm = c('metric', 'estimator', 'estimate')
  ) |>
    mutate(resolution = resolution, iteration = iteration) |>
    write.csv(paste0(results_root, '/tables/', fname, '.csv'),
              row.names = F)
  
  rm(df, df_auc, TrainIndex, Train, Test, predictions, cmRestrat)
  }

  eval_path <- paste0(results_root, '/tables/', fname, '.csv')
  eval_tab <- if(file.exists(eval_path)) read.csv(eval_path) else NULL

  if(predict_surface){

  pout <- paste0(results_root, '/suitability_maps')
  ###################      PREDICT ONTO SURFACE        #########################
  # the prediction is generally straightforward, it doesn't take an obscene
  # amount of RAM and can happen relatively quickly; overnight for the 
  # higher resolution data products.  
  
  # there is some weird corner case where we could not skip to the post prediction steps:
  # AOA, and SE. 
  
  if(!file.exists(file.path(pout, paste0(fname, '-Pr.tif')))){
  
    pr <- function(...) predict(..., type = 'response', num.threads = 1)$predictions 
  
    if(resolution %in% c('3arc', '1arc')){ntile = 1} else 
      if(resolution == '1-3arc'){ntile = 2} else {ntile = 4}  
  
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
    
     tiles <- file.path(tile_path_4pr, list.files(tile_path_4pr))
     pr_tile_path <- file.path(pout, paste0('pr_tiles', '_', resolution, iteration))
     dir.create(pr_tile_path, showWarnings = F)
    
     message('Writing Predicted Probability to Raster using tiles - this while take some time')
     pb <- txtProgressBar(
       min = 0,  max = length(tiles), style = 3, width = 50, char = "+")
     for (i in seq_along(tiles)){
       
       if(!file.exists(file.path(pr_tile_path, paste0(i, '.tif')))){
       terra::predict(
         rast(tiles[i]), rf_model, f = pr, 
         filename = file.path(pr_tile_path, paste0(i, '.tif')),
          overwrite = T, na.rm=TRUE)
        setTxtProgressBar(pb, i)
        gc(verbose = FALSE)
     }
     close(pb)
    }
    
     mos <- terra::mosaic(sprc(
       file.path(pr_tile_path, list.files(pr_tile_path))
     ), fun = 'mean')
      writeRaster(
        mos[[2]], 
        wopt = list(names = 'Probability'),
       filename = file.path(pout, paste0(fname, '-Pr.tif')))
    
    rm(mos)

    # the Pr.tif mosaic above is now the durable output; these scratch tiles (the split
    # predictor stack, and the per-tile predictions) aren't needed once it exists.
    unlink(tile_path_4pr, recursive = TRUE)
    unlink(pr_tile_path, recursive = TRUE)

    } else { # in these cases, we only need to use a single tile for prediction, we can
     # just use the existing virtual raster to do this. 
    
      if(!file.exists(file.path(pout, paste0(fname, '-Pr.tif')))){
      message('Writing Predicted Probability to Raster')
        terra::predict( # rasters, skip straight to modelling. 
          cores = 1, f = pr, 
          rast_dat, rf_model, cpkgs = "ranger",
          filename = file.path(pout, paste0(resolution, '-IterationCoarse', iteration, '.tif')),
          wopt = list(names = 'predicted_suitability'), na.rm=TRUE,
          overwrite = T)
    
        writeRaster( # we wrote a raster with both predicted class probabilities onto it. 
        # we don't want both, we are going to 'rewrite' the raster so only one class remains. 
          rast(
            file.path(
              pout, paste0(resolution, '-IterationCoarse', iteration, '.tif')))[[2]],
          file.path(pout, paste0(fname, '-Pr.tif')), 
          overwrite = T
      )
      file.remove(
        file.path(
          pout, paste0(resolution, '-IterationCoarse', iteration, '.tif')
          )
        )
      } else {
        'Predicted Probability Raster already exists - skipping.'}
    }
  
  gc(verbose = FALSE)

  }
  
  
  
  ################ AREA OF APPLICABILITY SURFACE    ############################
  message('Writing Area of Applicability to Raster')
  if(!exists('rf_model')){
    rf_model <- readRDS(paste0(results_root, '/models/', fname, '.rds'))
  }
  AOAfun <- function(model, r) {
    CAST::aoa(r, model, LPD = FALSE, verbose=FALSE)$AOA
  }

  terra::predict( # rasters, skip straight to modelling.
    cores = 1, f = AOAfun,
    rast_dat, model = model, cpkgs = "CAST",
    filename = file.path(pout, paste0(fname, '-AOA', '.tif')),
    wopt = list(names = 'AOA'), na.rm=TRUE,
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
    if(file.exists(tile_path_4se)){
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
      wopt = list(names = 'standardError'),
      filename = file.path(pout, paste0(fname, '-SE.tif')),
      overwrite = T)
    gc(verbose = FALSE)
    }
    unlink(final_se_path, recursive = TRUE)
    unlink(tile_path_4se, recursive = TRUE)
  }
  }

  invisible(list(fname = fname, metrics = eval_tab))
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
#' using the 'ideal' X & Y edges of the raster - which we know to be inaccurate by 
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
    file.path(root, paste0('dem_', res_string), 'dem.tif')
  )
  
  x['RasterCell'] <- terra::extract(rasta, x, cells = TRUE)$cell
  RC <- split(x, f = x$RasterCell) 
  
    # if there are multiple points per cell, 1st) discard the absences (if applicable)
    # 2) randomly sample out one of the presences.  
    
    select_rec <- function(x){
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
    
    # area cell, area quadrant * mean plants per life stage. 
    
    countR <- function(x){
      x$Prsnc_J <- ((res^2) / (3^2)) * mean(x$Prsnc_J)
      x$Prsnc_M <- ((res^2) / (3^2)) * mean(x$Prsnc_M)
      x$Prsnc_S <- ((res^2) / (3^2)) * mean(x$Prsnc_S, na.rm = TRUE)
      x <- x[sample(1:nrow(x), 1),] 
    }

    # let's train the models on the raw count data,not the transformed values. 

  dplyr::bind_rows(RC) |>
    dplyr::arrange(OBJECT) |>
    dplyr::select(-RasterCell) |>
    sf::st_as_sf()
}

#' Map a nominal resolution to the file-naming string used across results/ and data/
res_string <- function(resolution){
  switch(
    as.character(resolution),
    "3" = "3m",
    "10" = "1-3arc",
    "30" = "1arc",
    "90" = "3arc",
    stop("Input to `resolution` invalid. Need be one of: 3, 10, 30, 90")
  )
}

#' Class-stratified, prediction-domain-matched spatial train/test split
#'
#' @description Splits presence and absence points into train/test sets independently
#' using \code{CAST::knndm()} (Linnenbrink et al. 2024), each matched against the same
#' prediction-domain raster, then unions the two class-wise splits. A single class-blind
#' \code{knndm()} call risks stranding most of the rare, spatially clustered presence
#' class in one partition, starving (or flooding) Test of positives; splitting per class
#' and unioning keeps both classes proportionally represented in Train/Test while still
#' spatially blocking each against the true prediction area. This is meant to be computed
#' once per resolution and reused across every ratio/seed fit in that resolution's grid
#' search, rather than re-split randomly per fit - a random per-fit split is a bigger
#' source of PR-AUC variance for spatially clustered, rare-presence data than anything
#' about the PA ratio itself. \code{knndm()}'s train/test mode is marked experimental
#' upstream, so the achieved split is summarized (n, achieved test proportion, and the W
#' match statistic per class) via \code{message()} for a quick sanity check on every run,
#' and the raw \code{knndm} objects are returned for a closer look
#' (e.g. \code{plot(split$presence_knndm)}) per the package's own recommendation to
#' visually check the match before trusting it.
#' @param x sf point object with an \code{Occurrence} (0/1) column and a \code{.id}
#' column - a stable identifier that survives any later subsetting of \code{x}.
#' @param rast_dat SpatRaster. The prediction domain, passed to \code{knndm()} as
#' \code{modeldomain}.
#' @param train_split Numeric in (0, 1). Target fraction of points in Train, per class.
#' @param test_tolerance Numeric. Passed to \code{knndm()}'s \code{test_tolerance} - how
#' far \code{knndm} may search from the target test proportion for a configuration that
#' matches the prediction-domain distance distribution at all. This needs to stay loose:
#' knndm's achievable test proportions are set by how the points happen to cluster, and
#' finding any valid configuration can require a surprisingly wide search (observed:
#' presences here needed ~0.1 before a single valid configuration existed, landing near
#' 27% rather than the ~20% target). The exact final proportion is enforced afterward
#' regardless of where knndm itself lands within that search.
#' @return A list: \code{train_id}, \code{test_id} (the \code{.id} values assigned to
#' each, pooled across both classes, after forcing the exact target proportion),
#' \code{presence_knndm}, \code{absence_knndm} (the raw per-class \code{knndm} objects,
#' pre-forcing), and \code{summary} (achieved n/proportion/W per class, post-forcing).
spatial_class_split <- function(x, rast_dat, train_split, test_tolerance = 0.1){

  test_prop <- 1 - train_split

  split_one_class <- function(cls){
    pts <- x[x$Occurrence == cls, ]
    k <- CAST::knndm(pts, modeldomain = rast_dat, test_prop = test_prop, test_tolerance = test_tolerance)
    train_id <- pts$.id[k$indx_train]
    test_id  <- pts$.id[k$indx_test]

    # knndm's achievable proportions are discrete and can land well off target even when
    # test_tolerance is loose enough to find a configuration at all - force the exact
    # requested proportion by randomly moving the minimum number of points across.
    # A remainder from rounding goes to test (ceiling), not train.
    n <- length(train_id) + length(test_id)
    target_n_test <- ceiling(n * test_prop)
    if(length(test_id) > target_n_test){
      moved <- sample(test_id, length(test_id) - target_n_test)
      test_id <- setdiff(test_id, moved)
      train_id <- c(train_id, moved)
    } else if(length(test_id) < target_n_test){
      moved <- sample(train_id, target_n_test - length(test_id))
      train_id <- setdiff(train_id, moved)
      test_id <- c(test_id, moved)
    }

    list(train_id = train_id, test_id = test_id, knndm = k)
  }

  pres <- split_one_class(1)
  abs  <- split_one_class(0)

  summary <- data.frame(
    class = c('presence', 'absence'),
    n_train = c(length(pres$train_id), length(abs$train_id)),
    n_test  = c(length(pres$test_id), length(abs$test_id)),
    test_prop_achieved = c(
      length(pres$test_id) / (length(pres$train_id) + length(pres$test_id)),
      length(abs$test_id) / (length(abs$train_id) + length(abs$test_id))
    ),
    W = c(pres$knndm$W, abs$knndm$W)
  )
  message('Spatial train/test split (target test_prop = ', round(test_prop, 3), '):')
  message(paste(utils::capture.output(print(summary, row.names = FALSE)), collapse = '\n'))

  list(
    train_id = c(pres$train_id, abs$train_id),
    test_id = c(pres$test_id, abs$test_id),
    presence_knndm = pres$knndm,
    absence_knndm = abs$knndm,
    summary = summary
  )
}

#' Join real-population identity onto occurrence/background points
#'
#' @description Phase 0 of \code{population_rarefaction_roadmap.md}: primary assignment
#' is point-in-polygon against the \code{areas} layer's named population/site polygons
#' (\code{data/GroundTruthPts/ground_truth_median-VISITED.gpkg}). That layer was built for
#' ground-truthing/transect work, not as an exhaustive population census, so coverage is
#' incomplete - checked empirically against the iter1 3m presence pool, ~75% of points
#' fall inside a polygon; of the rest, ~55% sit within 500m of one and ~71% within 1km, but
#' the tail runs out to ~27km, and those are real, distinct, unsurveyed populations rather
#' than boundary/GPS slop. \code{max_dist} defaults to 500m to recover the boundary-slop
#' cases via nearest-polygon fallback without silently annexing a genuinely separate,
#' unsurveyed population into whichever named site happens to be closest.
#' @param x sf points (presence and/or absence/background).
#' @param areas sf polygons with a \code{Site} name column - the \code{areas} layer.
#' @param max_dist Numeric, meters. Fallback radius for the strict \code{population} field.
#' @return \code{x} with two new columns. \code{population}: \code{Site} name if inside a
#' polygon or within \code{max_dist} of one, else \code{NA} - use this wherever true
#' population membership matters (Phase 0/2/3's per-population rarefaction; dropping
#' \code{NA} rows there is the right call, per the roadmap, rather than guessing).
#' \code{population_zone}: always the single nearest \code{Site}, regardless of distance -
#' a full-landscape, no-\code{NA} partition meant only as a spatial blocking key for
#' \code{population_lopo_search()}, not a membership claim.
assign_population <- function(x, areas, max_dist = 500){

  x <- sf::st_transform(x, sf::st_crs(areas))
  within <- sf::st_intersects(x, areas)
  nearest_idx <- sf::st_nearest_feature(x, areas)
  nearest_dist <- as.numeric(sf::st_distance(x, areas[nearest_idx, ], by_element = TRUE))

  in_idx <- vapply(within, function(i) if(length(i)) i[1] else NA_integer_, integer(1))
  strict_idx <- ifelse(!is.na(in_idx), in_idx, ifelse(nearest_dist <= max_dist, nearest_idx, NA_integer_))

  x$population <- areas$Site[strict_idx]
  x$population_zone <- areas$Site[nearest_idx]
  x
}

#' Leave-one-population-out CV (Phase 1's population-blocked robustness check)
#'
#' @description Every presence-bearing population gets a turn as the sole test
#' population - all other points, presence and absence/background alike, are train,
#' zoned by \code{population_zone} (see \code{assign_population()}: the nearest named
#' population regardless of distance, so every location on the landscape falls in exactly
#' one zone and no point is left unassigned by this scheme). This replaces an earlier
#' single-random-block design: reusing one random subset of populations as test across
#' every seed replicate (the same pattern \code{spatial_class_split()} uses for
#' \code{adaptive_PAratio_search()}) would only characterize model/seed stochasticity for
#' whichever populations happened to land in test on that one draw, not variability from a
#' different, equally arbitrary draw - which is the thing a "held-out population"
#' robustness check is actually supposed to speak to. Giving every population a turn
#' removes that source of variance instead of estimating the risk of a bad one.
#' This does NOT re-run \code{adaptive_PAratio_search()}'s PA-ratio bracketing search per
#' population - it fits one fixed \code{(PAratio, distOrder)} per fold, normally whichever
#' ratio a prior \code{adaptive_PAratio_search()} call already chose as best for this
#' resolution/iteration (its \code{best_ratio}). That's deliberate: the expensive part of
#' any single fit here is Boruta feature selection plus \code{modeller()}'s knndm-tuned
#' hyperparameter grid, paid once per fit regardless of what's being compared - the
#' adaptive search already pays that cost dozens of times across its own ratio/seed grid,
#' so ~15-20 more single fits (one per population) is a drop in the bucket by comparison,
#' not a multiplier on the expensive part of the pipeline.
#' @param x sf presence+absence pool with an \code{Occurrence} (0/1) column - not yet
#' \code{.id}-tagged or population-joined; both happen inside.
#' @param resolution,distOrder,iteration,p2proc as passed through to
#' \code{distOrder_PAratio_simulator()}/\code{modeller()}.
#' @param PAratio Numeric. The single fixed ratio fit at every fold - see Description.
#' @param areas sf polygons with a \code{Site} column - passed to \code{assign_population()}.
#' @param max_dist Numeric, meters - passed to \code{assign_population()}. Only affects
#' which points get \code{NA} in the strict \code{population} field used for the
#' \code{min_presence} count below; folds themselves key off \code{population_zone}, which
#' has no \code{NA}s (see \code{assign_population()}).
#' @param min_presence Numeric. Populations with fewer than this many presence records
#' don't get a fold - too little held-out signal for a meaningful PR-AUC/ROC-AUC from a
#' handful of points.
#' @param base_seed Numeric. Each fold gets \code{base_seed + } its index. This isn't just
#' replicate bookkeeping - \code{modeller()}'s cache file stem includes the seed but not
#' the held-out population, so folds sharing a seed would collide: the first fold's fit
#' would get silently reloaded and reused for every later fold instead of each being fit
#' and evaluated against its own held-out population.
#' @param results_root As \code{modeller()}'s; point this at its own tree (e.g.
#' \code{results_populationsplit/}) rather than overwriting the primary results.
#' @return A list: \code{per_population} (data.frame, one row per held-out population:
#' \code{population}, \code{n_train}, \code{n_presence_test}, \code{n_absence_test},
#' \code{pr_auc}, \code{roc_auc} - a population whose zone happens to catch only one class
#' in test gets \code{NA} metrics rather than a fit, see Description) and \code{summary}
#' (\code{summarize_lopo()}'s presence-n-weighted aggregate over the scored folds - see
#' its own docs for why an unweighted mean across populations this size-mismatched would
#' be misleading). Both are also written to
#' \code{results_root/tables/<res>-Iteration<iteration>-LOPO-PA<ratio>DO:<distOrder>.csv}
#' and the same path with a \code{-summary} suffix.
population_lopo_search <- function(x, resolution, distOrder, iteration, p2proc, PAratio,
                                    areas, max_dist = 500, min_presence = 3,
                                    base_seed = 1000,
                                    results_root = file.path(PROJ_ROOT, 'results_populationsplit')){

  x <- dplyr::mutate(x, .id = seq_len(nrow(x)))
  x <- assign_population(x, areas, max_dist = max_dist)
  x_df <- sf::st_drop_geometry(x)

  pres_n <- stats::aggregate(
    .id ~ population_zone, data = x_df[x_df$Occurrence == 1, ], FUN = length)
  names(pres_n)[names(pres_n) == '.id'] <- 'n'
  pops <- sort(pres_n$population_zone[pres_n$n >= min_presence])

  message(
    length(pops), ' of ', nrow(pres_n), ' populations have >= ', min_presence,
    ' presence records and get a LOPO fold.')

  results <- vector('list', length(pops))
  for(i in seq_along(pops)){
    pop <- pops[i]
    test_id  <- x_df$.id[x_df$population_zone == pop]
    train_id <- x_df$.id[x_df$population_zone != pop]
    n_presence_test <- sum(x_df$Occurrence[x_df$.id %in% test_id] == 1)
    n_absence_test  <- sum(x_df$Occurrence[x_df$.id %in% test_id] == 0)

    # a population whose zone happens to catch zero background/absence points (nearest-
    # population zoning isn't population size, so this does happen for small/edge zones)
    # can't be evaluated - pr_auc/roc_auc need both classes present in Test - so skip it
    # rather than let modeller()'s dismo::evaluate() fail on an empty class.
    if(n_presence_test == 0 || n_absence_test == 0){
      message(
        'Skipping ', pop, ': test fold has ', n_presence_test, ' presence / ',
        n_absence_test, ' absence records - need both classes to evaluate.')
      results[[i]] <- data.frame(
        population = pop, n_train = length(train_id),
        n_presence_test = n_presence_test, n_absence_test = n_absence_test,
        pr_auc = NA_real_, roc_auc = NA_real_)
      next
    }

    split <- list(train_id = train_id, test_id = test_id)
    out <- distOrder_PAratio_simulator(
      x = x, distOrder = distOrder, PAratio = PAratio, resolution = resolution,
      seed = base_seed + i, predict_surface = FALSE, split = split,
      results_root = results_root, iteration = iteration, se_prediction = FALSE,
      train_split = NA, p2proc = p2proc
    )

    results[[i]] <- data.frame(
      population = pop,
      n_train = length(train_id),
      n_presence_test = n_presence_test,
      n_absence_test  = n_absence_test,
      pr_auc  = out$metrics$estimate[out$metrics$metric == 'pr_auc'],
      roc_auc = out$metrics$estimate[out$metrics$metric == 'roc_auc']
    )
  }
  results <- do.call(rbind, results)
  summary <- summarize_lopo(results)

  fstem <- paste0(results_root, '/tables/', res_string(resolution), '-Iteration', iteration,
                   '-LOPO-PA', PAratio, 'DO:', distOrder)
  write.csv(results, paste0(fstem, '.csv'), row.names = FALSE)
  write.csv(summary, paste0(fstem, '-summary.csv'), row.names = FALSE)

  message(
    'LOPO summary (presence-n-weighted, ', summary$n_populations, ' populations, ',
    summary$total_n_presence_test, ' held-out presence records total): pr_auc = ',
    round(summary$weighted_pr_auc, 3), ', roc_auc = ', round(summary$weighted_roc_auc, 3))

  list(per_population = results, summary = summary)
}

#' Presence-n-weighted summary of a \code{population_lopo_search()} table
#'
#' @description Folds vary enormously in size - in the real iter1 90m pool, Cochetopa's
#' held-out presence n is 50 vs. Bellview-Rustler/Gothic's n=3 - so an unweighted mean
#' across populations would let a 3-point fold's PR-AUC (which, with that few positives,
#' can only take a handful of discrete values and is essentially a coin flip) count the
#' same as a 50-point fold's genuinely well-estimated one. Weights each population's
#' contribution by its own \code{n_presence_test} instead. Rows with \code{NA} metrics
#' (a zero-class fold \code{population_lopo_search()} had to skip) are dropped first.
#' @param x The \code{per_population} data.frame \code{population_lopo_search()} returns.
#' @return A one-row data.frame: \code{n_populations} (folds actually scored, after
#' dropping \code{NA} rows), \code{total_n_presence_test}, \code{weighted_pr_auc},
#' \code{weighted_roc_auc}.
summarize_lopo <- function(x){
  x <- x[!is.na(x$pr_auc) & !is.na(x$roc_auc), ]
  data.frame(
    n_populations = nrow(x),
    total_n_presence_test = sum(x$n_presence_test),
    weighted_pr_auc = stats::weighted.mean(x$pr_auc, w = x$n_presence_test),
    weighted_roc_auc = stats::weighted.mean(x$roc_auc, w = x$n_presence_test)
  )
}

#' Generate data sets for SD modelling at different configurations of distOrders and PAratios
#' 
#' @description This function helps subset the input data (x) to combinations of distOrders- where
#' points within x distance of the raster tile resolution are removed, and PAratios the ratio
#' of presence to absence points - which is used to control the excess of absence points. 
#' @param x an sf/data frame/tibble. All potential records which could be used for modelling
#' @param distOrder Numeric. The multiple of the resolution to remove near absences by. One of: 0, 1, 2, 4, 8 etc.
#' @param PAratio Numeric. The denominator of the ratio, Presence is always 1, so 1.5 would indicate
#' 100 presence records and 150 absence records. 
#' @param resolution Numeric. The approximate resolution of the input data set, one of: 3, 10, 30, 90.
#' @param seed Numeric. Seeds both the absence subsample drawn below and (via passthrough
#' to \code{modeller()}) the train/test split and model fit, so repeated calls at the same
#' PAratio use independent replicates rather than the identical absence subsample every time.
#' @return Whatever \code{modeller()} returns, with one extra element: \code{PAratio} - the
#' ratio actually used. This equals the requested \code{PAratio} unless the absence pool
#' (after the distOrder filter) is too small to supply it, in which case it's capped to
#' every available absence and this reports what was actually achieved, so callers label
#' and record the fit by what it really is rather than what was asked for.
distOrder_PAratio_simulator <- function(x, distOrder, PAratio, resolution, seed = 1, split = NULL, ...){

  set.seed(seed)

  res_str <- res_string(resolution)

  # subset to absence points at distance intervals from presence points.
  # This accomplishes the distOrder parameter.
  abs <- x[x$Occurrence==0,]
  prs <- x[x$Occurrence==1,]
  abs <- abs[as.numeric(sf::st_distance(
    abs,
    prs[sf::st_nearest_feature(abs, prs),], by_element = TRUE
  )) > (distOrder*resolution),]

  # now recombine the data, we will sample down to get the PA ratio we need.
  x <- dplyr::bind_rows(prs, abs)

  # results indicate that 3:1 is pretty good, but very conservative in terms of suitable habitat
  # PAratio may be fractional (e.g. from adaptive_PAratio_search's bracketing), hence round().
  # Cap at the available pool rather than erroring: a resolution can have too few absences
  # left (post distOrder filter) to supply a high requested ratio, and duplicating points via
  # sampling with replacement would just be fake data, not a real ratio of that size.
  n_requested <- round(nrow(prs) * PAratio)
  if(n_requested > nrow(abs)){
    message(
      'Requested PAratio ', PAratio, ' needs ', n_requested, ' absences but only ',
      nrow(abs), ' are available after the distOrder filter; capping to ',
      round(nrow(abs) / nrow(prs), 3), '.')
    n_requested <- nrow(abs)
    PAratio <- nrow(abs) / nrow(prs)
  }
  abs <- abs[sample(1:nrow(abs), size = n_requested, replace = F),]
  x <- dplyr::bind_rows(prs, abs)

  out <- modeller(x, PAratio = paste0("1:", PAratio),
           resolution = res_str, distOrder = paste0('DO:', distOrder), seed = seed, split = split, ...)
  out$PAratio <- PAratio
  out
}

#' Adaptive coarse-to-fine search over presence:absence ratio, with seed replication
#'
#' @description Fits a bracketing search over PAratio instead of a fixed value: an
#' initial integer grid (\code{initial_ratios}) is fit \code{initial_replicates} times
#' each, so every candidate - not just the eventual winner - has a variance estimate and
#' the coarse-stage comparison isn't decided by a single lucky/unlucky split; whichever
#' ratios land within \code{good_frac} of the best mean holdout metric define a "good"
#' band; the midpoints between adjacent tested ratios inside that band become the next
#' stage's candidates. Each successive stage in \code{replicate_schedule} fits more seed
#' replicates per candidate, so the search narrows toward the optimum while also building
#' up an estimate of seed-to-seed performance variance near it. Every fit in this function
#' uses \code{predict_surface = FALSE} (fit + holdout evaluation only, which is fast for a
#' random forest) - only the single winning ratio is promoted to a full surface/SE
#' prediction at the end, using whichever of its replicates scored closest to the median
#' (not the best-of-N replicate: selecting on the same metric being reported would bias it
#' optimistically). Results are not ensembled across replicates.
#'
#' @param x,distOrder,p2proc,train_split as passed through to \code{distOrder_PAratio_simulator()}.
#' @param resolution Numeric. One of 3, 10, 30, 90; the same resolution is used for every
#' stage, since the optimal ratio is expected to depend on resolution.
#' @param initial_ratios Numeric vector. First-stage PAratio grid.
#' @param initial_replicates Numeric. Seeds to fit per ratio in the first-stage grid, so
#' every candidate ratio - not only the winner - has an estimate of seed-to-seed variance.
#' @param replicate_schedule Numeric vector, one entry per refinement stage after the
#' first; each value is how many seeds to fit per candidate ratio at that stage.
#' @param final_replicates Numeric. Additional seeds to fit at the single overall
#' best-performing ratio, on top of \code{initial_replicates} (or more, if it was also
#' visited during refinement).
#' @param good_frac Numeric in (0, 1]. A ratio's mean metric must be >= good_frac * the
#' best mean metric seen so far to fall in the "good" band used to pick the next stage's
#' midpoints.
#' @param metric Character, 'pr_auc' or 'roc_auc' - which row of the eval CSV modeller()
#' writes to compare candidate ratios on. Defaults to 'pr_auc' since presences are rare.
#' @param base_seed Numeric. Seeds are drawn sequentially as base_seed + 1, 2, 3, ... across
#' every fit this function performs, so no two fits (any ratio, any stage) share a seed.
#' @param test_tolerance Passed to \code{spatial_class_split()}. The Train/Test partition
#' is computed once up front (see Details) and reused for every ratio/seed fit; the exact
#' \code{1 - train_split} proportion is enforced regardless of this value (see
#' \code{spatial_class_split()}), so this only controls how hard \code{knndm} searches for
#' a spatially valid configuration to start from. Ignored when \code{spatial_split = FALSE}.
#' @param spatial_split Boolean. Defaults to TRUE (see Details). Set FALSE to instead use
#' modeller()'s classic \code{caret::createDataPartition()} split - a fresh random draw
#' per fit rather than one shared partition - for a direct comparison against the spatial
#' split's effect on seed-to-seed variance.
#' @param results_root Passed through to \code{modeller()}; see its documentation. Use
#' this to point a \code{spatial_split = FALSE} comparison run at its own directory rather
#' than overwriting the default results tree.
#'
#' Population-blocked robustness (Phase 1 of \code{population_rarefaction_roadmap.md}) is
#' NOT a mode of this function - see \code{population_lopo_search()} instead. A single
#' population-block draw reused across every ratio/seed replicate here (the same pattern
#' \code{spatial_split} uses) would only characterize model/seed stochasticity for
#' whichever populations happened to land in test on that one draw, not variability from a
#' different, equally arbitrary draw - the thing a "held-out population" robustness check
#' is actually supposed to speak to. \code{population_lopo_search()} sidesteps that by
#' giving every population a turn instead of estimating the risk of a bad one.
#' @details The Train/Test split is not re-randomized per fit. It's computed once, before
#' any ratio/seed is tried, via \code{spatial_class_split()} - a class-stratified
#' \code{CAST::knndm()} split matched against this resolution's prediction-domain raster -
#' and every fit in the search (and the final promoted model) reuses that same partition.
#' A random per-fit split was found to be a much bigger source of PR-AUC variance for this
#' rare, spatially clustered presence class than the PA ratio itself; fixing the split
#' isolates the variance the seed loop is meant to characterize (model/resampling
#' stochasticity) from variance that was really just an artifact of a noisy evaluation set.
#' @return A list: \code{results} (data.frame of every ratio/seed/metric fit), \code{best_ratio},
#' \code{best_seed} (the promoted, median-performing replicate), \code{fname} (the
#' promoted model's file stem), and \code{split} (the \code{spatial_class_split()} object
#' used throughout - inspect \code{split$summary} or \code{plot(split$presence_knndm)}
#' to confirm the split is trustworthy).
adaptive_PAratio_search <- function(x, resolution, distOrder, p2proc, train_split,
                                     iteration = 1,
                                     initial_ratios = 1:4,
                                     initial_replicates = 3,
                                     replicate_schedule = c(3, 5),
                                     final_replicates = 10,
                                     good_frac = 0.9,
                                     metric = 'pr_auc',
                                     base_seed = 1000,
                                     test_tolerance = 0.1,
                                     spatial_split = TRUE,
                                     results_root = file.path(PROJ_ROOT, 'results')){

  # tag every row with a stable id before any filtering/resampling, so the one-time
  # spatial split below (if used) can be looked up by id no matter which ratio/seed
  # subsets x later.
  x <- dplyr::mutate(x, .id = seq_len(nrow(x)))
  if(spatial_split){
    rast_dat <- rastReader(paste0('dem_', res_string(resolution)), p2proc)
    split <- spatial_class_split(x, rast_dat, train_split, test_tolerance)
  } else {
    # fall back to modeller()'s own random caret::createDataPartition() split, re-drawn
    # per fit from each fit's own seed - the classic approach, for comparison against the
    # spatial split above.
    split <- NULL
  }

  seed_counter <- base_seed
  results <- data.frame(ratio = numeric(0), seed = numeric(0), metric = numeric(0))

  fit_one <- function(ratio, seed){
    out <- distOrder_PAratio_simulator(
      x = x, distOrder = distOrder, PAratio = ratio, resolution = resolution,
      seed = seed, predict_surface = FALSE, split = split, results_root = results_root,
      iteration = iteration, se_prediction = FALSE, train_split = train_split, p2proc = p2proc
    )
    # out$PAratio is the ratio actually achieved, not necessarily the one requested - it's
    # capped if the resolution's absence pool can't supply the request (see
    # distOrder_PAratio_simulator()). Record what was really fit, not the label asked for.
    list(ratio = out$PAratio, metric = out$metrics$estimate[out$metrics$metric == metric])
  }

  run_stage <- function(ratios, n_rep){
    for(r in ratios){
      for(i in seq_len(n_rep)){
        seed_counter <<- seed_counter + 1
        fit <- fit_one(r, seed_counter)
        results <<- rbind(results, data.frame(ratio = fit$ratio, seed = seed_counter, metric = fit$metric))
      }
    }
  }

  # midpoints between adjacent tested ratios that both fall within the current good band
  next_candidates <- function(){
    agg <- aggregate(metric ~ ratio, results, mean)
    best <- max(agg$metric)
    good <- sort(agg$ratio[agg$metric >= good_frac * best])
    if(length(good) < 2) return(numeric(0))
    tested <- sort(unique(results$ratio))
    in_band <- tested[tested >= min(good) & tested <= max(good)]
    if(length(in_band) < 2) return(numeric(0))
    mids <- (head(in_band, -1) + tail(in_band, -1)) / 2
    mids[!mids %in% tested]
  }

  # stage 0: coarse integer grid, replicated so every candidate (not just the
  # eventual winner) has a variance estimate to compare on.
  run_stage(initial_ratios, initial_replicates)

  # refinement stages: bisect toward the good band, replicating more each time
  for(n_rep in replicate_schedule){
    cands <- next_candidates()
    if(length(cands) == 0) break
    run_stage(cands, n_rep)
  }

  # final stage: replicate the single best-performing ratio to characterize
  # seed-to-seed uncertainty at the winner.
  agg <- aggregate(metric ~ ratio, results, mean)
  best_ratio <- agg$ratio[which.max(agg$metric)]
  run_stage(best_ratio, final_replicates)

  # promote the winning ratio's median-performing replicate - not the best-of-N - to a
  # full surface + SE prediction; picking the best score would optimistically bias the
  # metric we report, since it'd be chosen by the same statistic being reported. No
  # ensembling across replicates. distOrder_PAratio_simulator() re-derives the same
  # fname for this (ratio, seed) pair, so modeller() finds the cached fit and skips
  # straight to prediction rather than refitting.
  winners <- results[results$ratio == best_ratio, ]
  best_seed <- winners$seed[which.min(abs(winners$metric - median(winners$metric)))]

  # persist the full search history (every ratio/seed/metric tried) - this is the
  # evidence behind best_ratio, and the source data for a PR-AUC-vs-ratio plot - before
  # the raster prediction step below, which is the slowest and most failure-prone part
  # of this function (full-surface + AOA + SE prediction). If it crashes, the search
  # itself (the expensive, already-complete part) isn't lost along with it.
  write.csv(
    results, paste0(results_root, '/tables/', resolution, '-Iteration', iteration, '-PAratio-search.csv'),
    row.names = FALSE)

  final <- distOrder_PAratio_simulator(
    x = x, distOrder = distOrder, PAratio = best_ratio, resolution = resolution,
    seed = best_seed, predict_surface = TRUE, se_prediction = TRUE, split = split,
    results_root = results_root, iteration = iteration, train_split = train_split, p2proc = p2proc
  )

  list(results = results, best_ratio = best_ratio, best_seed = best_seed,
       fname = final$fname, split = split)
}

#' Convert a cross validation object generated by CAST into an rsample object
#' 
#' @description This function is for mix and matching Caret with tidymodels, essentially
#' allowing you to use carets RFE and then pushing to tidymodels for various racing
#' approaches for early exit while training hyperparameters. Also spatialsample
#' CV object is a massive disappointment, so keep yourself happy! Use CAST. 
#' @param x an object from CAST::knndm (others probably supported!)
#' @param train the data which CAST split from. 
CAST2rsample <- function(x, train){
  
  indices_L <- lapply(
    vector(mode = 'list', length = length(x$indx_train)), 
    \(x) setNames(vector(mode = 'numeric', 2), c('analysis', 'assessment'))
    )
  
  for (i in seq_along(indices_L)){
    indices_L[[i]]['analysis'] = sapply(x$indx_train, unlist)[i] 
    indices_L[[i]]['assessment'] = sapply(x$indx_test, unlist)[i]
  }
  
  splits <- lapply(indices_L, rsample::make_splits, data = train)
  rsample::manual_rset(splits, paste('Split', 1:length(x$indx_train)))
}

#' split data into test, and train, and generate CV folds too. 
#' @param dataframe which for modelling
#' @param bn file basename, this fn will strip out the resolution prefix and search for 
#' an existing file, and load the index from disk if found. Otherwise perform twinning, 
#' this is to be extra sure that each set of models us using at the least, the same test 
#' and train split (which should be accomplished alone by using the u1 argument to `twin`), 
#' but given how long hyperparam tuning takes we want to be extra sure. 
splitData <- function(df, fp, bn){
  
  # We will create three columns for our data, which can then be used
  # to separate the data sets into two sets, where the longitude, latitude, and response
  # variables are almost identical - like twins. Basically, we have a problem 
  # which makes splitting along the outcome variable, not an ideal choice. If you 
  # remember, the southern Cocheotopa Dome (CD) population has drastically higher 
  # counts of plants than the populations near CB, 10x individuals
  # not being uncommon! So when we 'split' our data, what we end up with is a gradient
  # which is actually a mix of CB:CD and then quickly just becomes CD. So our independent
  # test set isn't really what we see across the species, rather it's evaluating two
  # distinct components. 
  
  # fortunately for these predictions we are only including the plot level data and 
  # some 'near' absences. so we can try and 'ameliorate' this split gradients using
  # three steps: rescale longitude and latitude from 0:1, and do the same with our
  # counts. When we then 'combine' the three data sets, on paper, we should get
  # be able to incorporate a SMIDGE of each of these aspects to the split. Although
  # we will not entirely remove the CD area having more plants, I think we will slightly
  # mute the effect. 
  f <- file.path(fp, 'test_data', paste0('twin_indx-', gsub('-I.*$', '', bn), '.txt'))
  
  if(!file.exists(f)){
    
    twinning_dat <- df |>
      dplyr::mutate( # this algo PROBABLY rescales too, but we'll just feed em in to be sure. 
        x = scales::rescale(st_coordinates(df)[,1]), # the stats doc is rich the tech 
        y = scales::rescale(st_coordinates(df)[,2]),  # not so much 
        Prsnc = scales::rescale(Prsnc_All)
      ) |>
      sf::st_drop_geometry() |>
      dplyr::select(Prsnc, y, x)
    
    indx <- twinning::twin(twinning_dat, r=5, u1 = 2)
    cat(indx, file = f)
    
  } else {
    indx <- as.numeric(unlist(strsplit(readLines(f, warn = FALSE), ' ')))
  }
  
  train <- df[-indx,]
  test  <- df[indx,]
  
  # now that we have our training data, and our test data that we will compare our
  # final model predictions too, we need cross validation folds for model selection, 
  # hyperparameter tuning and to fit our model. 
  
  # We will use a spatial CV structure for this. Because we don't have that many
  # records to parse through we'll use nndm - but more on that later. 
  
  # the spatial CV will use the entire area of prediction when determining layouts
  # in space. Predicting XGBoost is pretty intense, and we really only have count
  # data from pretty limited areas. Let's restrict our predictions to adjacent areas
  # because the models will be saved, someone could re predict them in the future if
  # they wanted further testing. 
  
  km <- kmeans( df[,c('Longitude', 'Latitude')] |> sf::st_drop_geometry(),
                2, iter.max = 10, nstart = 1)
  df$Cluster <- km$cluster
  
  clusts <- split(df, f = df$Cluster)
  bbs <- lapply(clusts, \(x) sf::st_union(x) |>
                  sf::st_transform(5070) |>
                  sf::st_buffer(5000) |>
                  sf::st_bbox() |>
                  sf::st_as_sfc()
  ) |>
    dplyr::bind_rows() |>
    t() |>
    sf::st_as_sfc(crs=5070) |>
    sf::st_union()
  
  train.sf <- sf::st_as_sf(train, coords = c('Longitude', 'Latitude'), crs = 32613)
  # Fixed-k (not CAST::nndm()'s near-leave-one-out - fold count == nrow(train),
  # confirmed by testing, same finding density_bayes.R's brms_spatial_fold_ids()
  # made for the Bayesian candidates) - RFE (caret::rfeControl(index=...)) and the
  # ML candidates' spatial-CV resampling (CAST2rsample() below) both consume this,
  # and near-LOO scale made RFE alone take over a day on this dataset. knndm()'s
  # indx_train/indx_test are the same list-of-per-fold-indices shape nndm() returns,
  # so this is a drop-in swap - just k=10 folds instead of ~n.
  indices_nndm_CAST <- CAST::knndm(
    train.sf,
    sf::st_transform(bbs, sf::st_crs(train.sf)),
    k = 10,
    samplesize = 1000)

  return(
    list(
      nndm_indices = indices_nndm_CAST, 
      train = train, 
      test = test,
      train.sf = train.sf
    )
  )
}

#' fit poisson xgboost models to the data
#' 
#' @param rec a tidymodels recipe 
#' @param cv a cross validation structure, such as from `rsample` or `CAST`.  
#' @param train data split
#' @param test data split
#' @param tune_gr tuning grid. 
poiss <- function(rec, cv, train, test, tune_gr){
  
  xgb_poisson <- tune_gr |>
    parsnip::set_engine("xgboost", objective = "count:poisson") 
  
  xgb_poisson_gr <- xgb_poisson |>
    tune::extract_parameter_set_dials() |>
    dials::grid_regular(levels = 3)
  
  future::plan(multisession, workers = parallel::detectCores())
  params_xg_poisson <- xgb_poisson |>
    finetune::tune_race_anova(
      rec,
      metrics = yardstick::metric_set(yardstick::mae), 
      resamples = cv,
      grid = xgb_poisson_gr
    )
  
  best_xg_pois <- tune::select_best(params_xg_poisson, metric = "mae")
  xg_poisson <- xgb_poisson |>
    tune::finalize_model(best_xg_pois) |>
    fit(Prsnc_All ~ ., data = train)
  
  preds <- data.frame(
    'Observed' = test$Prsnc_All, 
    'Predicted' = stats::predict(xg_poisson, new_data = test),
    'Pr.suit' = test$Pr.SuitHab
  )
  preds <- setNames(preds, c('Observed', 'Predicted', 'Pr.suit'))
  
  return(
    list(
      Model = xg_poisson, 
      Predictions = preds
    )
  )
}

#' fit tweedie models to the data
#' 
#' @param rec a tidymodels recipe 
#' @param cv a cross validation structure, such as from `rsample` or `CAST`.  
#' @param train data split
#' @param test data split
#' @param tune_gr tuning grid. 
tweed <- function(rec, cv, train, test, tune_gr){
  
  xgb_tweedie_model <- tune_gr |>
    parsnip::set_engine("xgboost", objective = "reg:tweedie")
  
  xgb_tweedie_gr <- xgb_tweedie_model |> 
    tune::extract_parameter_set_dials() |> 
    dials::grid_regular(levels = 3)
  
  future::plan(multisession, workers = parallel::detectCores()) 
  xbg_tweedie_params <- xgb_tweedie_model |> 
    finetune::tune_race_anova(
      rec,
      resamples = cv,
      metrics = yardstick::metric_set(yardstick::mae), 
      grid = xgb_tweedie_gr
    )
  
  best_param_tweedie <- tune::select_best(xbg_tweedie_params, metric = "mae")
  
  xgb_tweedie_mod <- xgb_tweedie_model |>
    tune::finalize_model(best_param_tweedie)
  
  xgb_tweedie_fit <- xgb_tweedie_mod %>% 
    fit(Prsnc_All ~ ., data = train)
  
  preds <- data.frame(
    'Observed' = test$Prsnc_All, 
    'Predicted' = stats::predict(xgb_tweedie_fit, new_data = test),
    'Pr.suit' = test$Pr.SuitHab
  )
  preds <- setNames(preds, c('Observed', 'Predicted', 'Pr.suit'))
  
  return(
    list(
      Model = xgb_tweedie_fit, 
      Predictions = preds
    )
  )
}

#' calculate some summary values about the regression models 
mets <- function(x){
  data.frame(
    Metric = c('MAE', 'MSE', 'RMSE'),
    Value = c(
      MAE  = Metrics::mae(x$Observed, x$Predicted),
      MSE  = Metrics::mse(x$Observed, x$Predicted),
      RMSE = Metrics::rmse(x$Observed, x$Predicted)
    )
  )
}

#' Estimate the number of plants per raster cell. 
#'
#' @description
#' @param x Data frame of occurrences.
#' @param fp base file path to directory to save contents. 
#' @param bn Character. base file name which will unambiguously identify the objects saved from 
#' this function (the model, an evaluation table, variable selection object).
#' @param brms_families NULL (default: brms_cv_compare()'s own 4-family list)
#' or a named list of brms family objects, to scope Phase 1's Bayesian
#' comparison down for a cheap smoke test instead of the full sweep.
#' @param brms_refit_args list of extra arguments passed to
#' brms_promote_and_refit() (e.g. list(full_chains = 2, full_iter = 1000,
#' full_warmup = 500) for a lighter final refit than the 4/2000/1000 default).
densityModeller <- function(x, bn, fp, brms_families = NULL, brms_refit_args = list()){

  t_total <- Sys.time()

  # split the data into train/test and spatial CV.
  dsplit <- splitData(x, fp = fp, bn = bn)
  
  train <- dsplit$train
  test <- dsplit$test
  nndm_indices <- dsplit$nndm_indices
  train.sf <- dsplit$train.sf
  
  ##############       both null models are calculated up here    ##############
  # arithmetic means by population. A null model. 
  preds <- train |>
    dplyr::group_by(Lctn_bb) |>
    dplyr::summarize(Predicted = mean(Prsnc_All)) |>
    sf::st_drop_geometry()
  
  arith_mean <- test |> 
    dplyr::select(Observed = Prsnc_All, Lctn_bb) 
  mean_preds <- dplyr::left_join(arith_mean, preds, by = 'Lctn_bb')
  
  # kriging interpolation, a spatial null model. 
  krig_preds <- data.frame(
    'Observed' = test$Prsnc_All, 
    'Predicted' = gstat::krige(Prsnc_All ~ 1, train.sf, newdata = test)$var1.pred,
    'Pr.suit' = test$Pr.SuitHab
  )
  
  # this was the grouping variable required for the arithmetic mean, we drop it now.
  train <- sf::st_drop_geometry(train) |>
    select(-Lctn_bb)
  test <- sf::st_drop_geometry(test)

  # RFE below may drop Longitude/Latitude as non-predictive features; keep a copy
  # so the brms candidates' gp() term (needs coordinates, not a "selected" feature) can
  # reattach them later.
  train_coords <- train[, c('Longitude', 'Latitude')]

  # feature selection - timed: this is the step the near-LOO nndm() folds
  # made take >1 day before splitData() switched to fixed-k knndm().
  t0 <- Sys.time()
  if(!file.exists(file.path(fp, 'modelsTune', paste0(bn, '.rds')))){

    ctrl <- caret::rfeControl(
      functions = caret::treebagFuncs,
      method = "cv",
      index = nndm_indices$indx_train,
      allowParallel = TRUE)

    future::plan(future::multisession, workers = parallel::detectCores())
    rfProfile <- caret::rfe(
      x = train[,2:ncol(train)], y = train$Prsnc_All,
      rfeControl = ctrl, metric = 'MAE')

    saveRDS(rfProfile, file.path(fp, 'modelsTune', paste0(bn, '.rds')))
  } else {
    rfProfile <- readRDS(file.path(fp, 'modelsTune', paste0(bn, '.rds')))
  }
  message(sprintf('[timing] RFE: %.1fs', as.numeric(difftime(Sys.time(), t0, units = 'secs'))))

  train <- dplyr::select(train, all_of(c('Prsnc_All', predictors(rfProfile))))
  
  rec <- recipes::recipe(Prsnc_All ~ ., data = train) 
  rs <- rsample::bootstraps(train, 10)
  indx_nndm_rs <- CAST2rsample(nndm_indices, train)
  
  # now define a tuning grid for trees.  
  tune_gr <- parsnip::boost_tree( 
    mode = 'regression',
    trees = tune(),
    tree_depth = tune(),
    min_n = tune(), 
    loss_reduction = tune(),
    learn_rate = tune(),
    stop_iter = tune()
  )
  
  # tune hyper parameters and fit all models  - the hyper param tuning on occasion
  # is super slow, and may crash so we want to write to disk as they are completed.
  if(missing(fp)){fp <- file.path(PROJ_ROOT, 'results', 'count_models')}

  t0 <- Sys.time()
  f <- file.path(fp, 'models', paste0(bn, '-poisson_spat.rds'))
  if(!file.exists(f)){
    poiss_spat_cv <- poiss(rec, indx_nndm_rs, train, test, tune_gr)
    saveRDS(poiss_spat_cv, f)
  } else {poiss_spat_cv <- readRDS(f)}
  message(sprintf('[timing] XGB Poisson Spat.: %.1fs', as.numeric(difftime(Sys.time(), t0, units = 'secs'))))

  t0 <- Sys.time()
  f <- file.path(fp, 'models', paste0(bn, '-poisson.rds'))
  if(!file.exists(f)){
    poiss_cv <- poiss(rec, rs, train, test, tune_gr)
    saveRDS(poiss_cv, f)
  } else {poiss_cv <- readRDS(f)}
  message(sprintf('[timing] XGB Poisson: %.1fs', as.numeric(difftime(Sys.time(), t0, units = 'secs'))))

  t0 <- Sys.time()
  f <- file.path(fp, 'models', paste0(bn, '-tweedie_spat.rds'))
  if(!file.exists(f)){
    tweedie_spat_cv <- tweed(rec, indx_nndm_rs, train, test, tune_gr)
    saveRDS(tweedie_spat_cv, f)
  } else {tweedie_spat_cv <- readRDS(f)}
  message(sprintf('[timing] XGB Tweedie Spat.: %.1fs', as.numeric(difftime(Sys.time(), t0, units = 'secs'))))

  t0 <- Sys.time()
  f <- file.path(fp, 'models', paste0(bn, '-tweedie.rds'))
  if(!file.exists(f)){
    tweedie_cv <- tweed(rec, rs, train, test, tune_gr)
    saveRDS(tweedie_cv, f)
  } else {tweedie_cv <- readRDS(f)}
  message(sprintf('[timing] XGB Tweedie: %.1fs', as.numeric(difftime(Sys.time(), t0, units = 'secs'))))

  tune_gr <- parsnip::boost_tree(
    mode = 'regression',
    trees = tune(),
    min_n = tune(),
    tree_depth = tune()
  )

  t0 <- Sys.time()
  f <- file.path(fp, 'models', paste0(bn, '-lgbm-poisson_spat.rds'))
  if(!file.exists(f)){
    lgbm_cv <- gbs(rec, indx_nndm_rs, train, test, tune_gr, mode = 'regression', metric = 'mae',
                     engine = 'lightgbm', objective = 'poisson', resp = 'Prsnc_All')
    saveRDS(lgbm_cv, f)
  } else {lgbm_cv <- readRDS(f)}
  message(sprintf('[timing] LGBM Poisson Spat.: %.1fs', as.numeric(difftime(Sys.time(), t0, units = 'secs'))))

  # Bayesian candidates (census_uncertainty_roadmap.md Phase 1, scripts/Census_uncertainty/density_bayes.R).
  # RFE (above) may have dropped Longitude/Latitude from `train` since they're
  # coordinates, not predictive features - reattach them here, they're required by the
  # brms candidates' gp() term but aren't evaluated as XGBoost/LightGBM features.
  train_brms <- train
  train_brms$Longitude <- train_coords[,1]
  train_brms$Latitude  <- train_coords[,2]

  t0 <- Sys.time()
  brms_cv <- if (is.null(brms_families)) {
    brms_cv_compare(train_brms, test)
  } else {
    brms_cv_compare(train_brms, test, families = brms_families)
  }
  message(sprintf('[timing] brms_cv_compare (%d family x 2 CV modes): %.1fs',
                   if (is.null(brms_families)) 4 else length(brms_families),
                   as.numeric(difftime(Sys.time(), t0, units = 'secs'))))

  t0 <- Sys.time()
  brms_final <- do.call(brms_promote_and_refit, c(
    list(
      best_family = switch(brms_cv$best$family_name,
             Poisson = poisson(), NegBinomial = brms::negbinomial(),
             HurdlePoisson = brms::hurdle_poisson(), HurdleNegBinomial = brms::hurdle_negbinomial()),
      train = train_brms, seed = 1, fp = fp, bn = bn, family_label = tolower(brms_cv$best$family_name)
    ),
    brms_refit_args
  ))
  message(sprintf('[timing] brms_promote_and_refit (%s): %.1fs', brms_cv$best$family_name,
                   as.numeric(difftime(Sys.time(), t0, units = 'secs'))))

  # now calculate the evaluation statistics.
  namev <- c('Arithmetic Mean', 'Kriging',
             'XGB Poisson Spat.', 'XGB Poisson', 'XBG Tweedie Spat.',  'XGB Tweedie', 'LGBM Poisson Spat.')
  mods <- list(mean_preds, krig_preds, poiss_spat_cv$Predictions, poiss_cv$Predictions,
               tweedie_spat_cv$Predictions, tweedie_cv$Predictions, lgbm_cv$Predictions)

  metrrs <- lapply(mods, mets) |>
    dplyr::bind_rows() |>
    dplyr::mutate(Model = rep(namev, each = 3), .before = 1) |>
    dplyr::bind_rows(brms_cv$table)

  # now we will save the models, evaluation table, and information, note we
  # also save the variable selection object for AOA calculation

  write.csv(metrrs, file.path(fp, 'tables', paste0(bn, '.csv')), row.names = FALSE)

  message(sprintf('[timing] densityModeller() TOTAL: %.1fs', as.numeric(difftime(Sys.time(), t_total, units = 'secs'))))

  invisible(list(metrics = metrrs, brms_promoted = brms_final))
}

#' fit tweedie models to the data
#' 
#' @param rec a tidymodels recipe 
#' @param cv a cross validation structure, such as from `rsample` or `CAST`.  
#' @param train data split
#' @param test data split
#' @param response Character string. Name of the column holding the response variable. Unquoted. 
#' @param tune_gr tuning grid.
#' @param mode Character. Argument to `parsnip::boost_tree` mode, defaults to 'regression'. 
#' @param engine Character. Argument to `parsnip::set_engine`, engine. 
#' @param objective Character. Argument to `parsnip::set_engine`, objective.  
#' @param levels Numeric. Argument to `dials::grid_regular`  
#' @param metric Character. Evaluation metric for the model passed onto `yardstick::metric_set`
gbs <- function(rec, cv, train, test, resp, tune_gr, mode, engine, objective, levels, metric){
  
  if(missing(levels)){levels <-3}
  
  model <- parsnip::set_engine(tune_gr, engine = engine, objective = objective)
  
  gr <- model |>
    tune::extract_parameter_set_dials() |>
    dials::grid_regular(levels = levels)
  
  future::plan(multisession, workers = parallel::detectCores()/2)
  params <- model |>
    finetune::tune_race_anova(
      rec,
      resamples = cv,
      metrics = yardstick::metric_set(yardstick::mae),
      grid = gr
    )
  
  best_param <- tune::select_best(params, metric = metric)
  
  final_model <- model |>
    tune::finalize_model(best_param)
  
  fit <- final_model %>%
    workflows::fit(Prsnc_All ~ ., data = train)
  
  preds <- data.frame(
    'Observed' = test[[resp]], 
    'Predicted' = stats::predict(fit, new_data = test),
    'Pr.suit' = test[['Pr.SuitHab']]
  ) |>
    setNames(c('Observed', 'Predicted', 'Pr.suit'))
  
  return(
    list(
      Model = fit, 
      Predictions = preds
    )
  )
}

#' fit xgboosted models to predict plant count per space density.
#' @param x a vector of predicted suitable habitat rasters to use for input.
#' @param return_early logical, if TRUE (default) return the assembled plot-level
#' data.frame without fitting any density models - the behavior `ModelInterpretation.Rmd`
#' relies on for PDP-plot exploration. Set FALSE to actually run `densityModeller()`.
wrapper <- function(x, return_early = TRUE){

  res <- gsub('-I.*$', '', x)
  res_string <- switch(res,
                       "3m" = "3m",
                       "1-3arc" = "10m",
                       "1arc" = "30m",
                       "3arc" = "90m",
                       stop("Input to `res` invalid. Need be one of: 3, 10, 30, 90")
  )
  
  ct <- sf::st_read(
    file.path(PROJ_ROOT, 'data', 'Data4modelling', paste0(res_string, '-count-iter1.gpkg')
    ), quiet = TRUE
  ) |>
    dplyr::select(Prsnc_M, Prsnc_J, Lctn_bb)
  
  p2proc = paste0(PROJ_ROOT, '/data/spatial/processed')
  rast_dat <- rastReader(paste0('dem_', res), p2proc) 
  
  r <- terra::rast(file.path(PROJ_ROOT, 'results', 'suitability_maps', x))
  
  df <- dplyr::bind_cols(
    ct, 
    dplyr::select(terra::extract(rast_dat, ct), -ID), 
    Pr.SuitHab = terra::extract(r, ct)[,2]
  ) |> 
    filter(Lctn_bb != 'PABA') |>
    mutate(Prsnc_All = Prsnc_J + Prsnc_M, .before = 1) |>
    tidyr::drop_na() |>
    select(-Prsnc_J, -Prsnc_M)
  
  ## this just here for playing wiht PDP plots!!!
  if(return_early){ return(df) }

  coords <- sf::st_coordinates(df)
  df$Longitude <- coords[,1]
  df$Latitude <- coords[,2]

  densityModeller(df, fp = paste0(PROJ_ROOT, '/results/count_models'), bn = gsub('DO.*$', '', x))

}


patchAttributes <- function(x, region_k = 2, region_buffer_m = 5000, region_max_cells = 12e6){

  pr <- terra::rast(x[['Pr']])
  aoa <- terra::rast(x[['AOA']])
  threshs <- terra::rast(x[['thresholds']])

  # will iterate through each of the options.
  evals <- data.frame(
    eval = c('spec_sens', 'equal_sens_spec', 'sensitivity'),
    layer = 1:3
  )

  threshs <- terra::mask(threshs, aoa, maskvalues = 0)

  # --- crop to geographically disjoint regions around known occurrences,
  # same approach (and same clusterRegionExtents()) as CostDistances.R uses
  # for gdistance::transition(). The full 1-3arc domain (~330M cells) is far
  # too large for landscapemetrics::get_patches()/calculate_lsm() to run on
  # directly - occupied ground is a tiny fraction of the full grid, and an
  # unscoped run on this was observed to run 10+ hours with no guarantee of
  # ever finishing safely.
  occ_pts <- sf::st_read(
    file.path(PROJ_ROOT, 'data', 'Data4modelling', '3m-presence-iter1.gpkg'), quiet = TRUE) |>
    dplyr::filter(Presenc == 1) |>
    sf::st_transform(terra::crs(pr))

  region_exts <- clusterRegionExtents(
    occ_pts, k = region_k, buffer_m = region_buffer_m,
    max_cells = region_max_cells, res_m = mean(terra::res(pr)))

  landscape_tiles <- vector('list', length(region_exts))
  calcs_list      <- vector('list', length(region_exts))
  occ_list        <- vector('list', length(region_exts))
  central_list    <- vector('list', length(region_exts))
  id_offset       <- 0

  # read once - reused against every region below; a point outside a given
  # region's cropped extent just extracts NA and is dropped. Read via sf and
  # filtered to Presenc == 1, same convention as elsewhere in this pipeline
  # (e.g. sdm_groundtruth_design_1-3arc.R) - the file also carries 9 stray
  # GEOMETRYCOLLECTION rows (NA Presenc) alongside the 679 real POINT rows;
  # terra::vect() on the raw file guesses a single geometry type from that
  # mix and silently returns the *polygons* it finds inside those 9 rows,
  # dropping every real point. Filtering first sidesteps the type-guessing
  # entirely.
  pres <- sf::st_read(
    file.path(PROJ_ROOT, 'data', 'Data4modelling', '3m-presence-iter1.gpkg'), quiet = TRUE) |>
    dplyr::filter(Presenc == 1) |>
    terra::vect()

  for (r in seq_along(region_exts)) {

    reg_ext     <- region_exts[[r]]
    threshs_reg <- terra::crop(threshs, reg_ext)

    calcs_reg <- vector(mode = 'list', length = 3)
    names(calcs_reg) <- evals$eval

    for(i in seq_along(evals$eval)){

      calcs_reg[[i]] <- landscapemetrics::calculate_lsm(
        threshs_reg[[i]],
        what = c("lsm_p_area",  "lsm_p_enn", "lsm_p_cai", "lsm_p_para", "lsm_p_frac"),
        neighbourhood = 4, directions = 8) |>
        dplyr::select(id, metric, value)

    }

    patches_reg <- landscapemetrics::get_patches(threshs_reg, directions = 8)
    calcs_reg <- dplyr::bind_rows(calcs_reg, .id = 'evals')

    # kind of strange object, we'll just brute force it back to a happy terra object
    landscape_reg <- c(
      patches_reg$layer_1$class_1,
      patches_reg$layer_2$class_1,
      patches_reg$layer_3$class_1
    )
    names(landscape_reg) <- evals$eval

    # offset patch IDs so they stay globally unique once regions are remerged
    max_id_reg <- suppressWarnings(max(calcs_reg$id, na.rm = TRUE))
    if (is.finite(max_id_reg)) {
      landscape_reg <- landscape_reg + id_offset
      calcs_reg$id  <- calcs_reg$id + id_offset
      id_offset     <- id_offset + max_id_reg
    }

    landscape_tiles[[r]] <- landscape_reg
    calcs_list[[r]]      <- calcs_reg

    # --- occupied/unoccupied + central tendency, computed HERE while still
    # region-cropped. Running extract()/zonal() against the full 330M-cell
    # domain (post-extend()) is what OOM-killed (127GB RSS) the first pass
    # at this - zonal stats grouped by a huge-cardinality categorical raster
    # forces a full memory-resident scan. Confined to a cropped region, both
    # are cheap, and results across regions are just row-bound after.
    occ_list[[r]] <- terra::extract(landscape_reg, pres) |>
      dplyr::select(-ID) |>
      tidyr::pivot_longer(everything(), values_to = 'patch', names_to = 'threshold') |>
      tidyr::drop_na()

    pr_reg <- terra::crop(pr, reg_ext) |>
      terra::mask(terra::crop(aoa, reg_ext), maskvalues = 0)

    central_reg <- vector(mode = 'list', length = 3)
    for (i in seq_along(evals$eval)) {
      pr_sub <- terra::mask(pr_reg, threshs_reg[[i]], maskvalues = 0)
      mn  <- terra::zonal(pr_sub, landscape_reg[[i]], fun = 'mean',   na.rm = TRUE)
      mdn <- terra::zonal(pr_sub, landscape_reg[[i]], fun = 'median', na.rm = TRUE)
      central_reg[[i]] <- setNames(
        data.frame(cbind(mn, mdn[, 2], evals = evals$eval[[i]])),
        c('patchID', 'mean', 'median', 'evals'))
    }
    central_list[[r]] <- dplyr::bind_rows(central_reg)
  }

  landscape <- if (length(landscape_tiles) > 1) {
    terra::mosaic(terra::sprc(landscape_tiles), fun = 'first')
  } else {
    landscape_tiles[[1]]
  }
  landscape <- terra::extend(landscape, pr)
  names(landscape) <- evals$eval

  calcs <- dplyr::bind_rows(calcs_list)

  write.csv(
    calcs,
    file = file.path(PROJ_ROOT, 'results', 'patch_summaries', paste0(x[['version']][1], '-patches.csv')),
    row.names = F)

  writeRaster(
    landscape, overwrite = TRUE,
    filename = file.path(PROJ_ROOT, 'results', 'patches', paste0(x[['version']][1], '-patches.tif')))

  # now make a table noting which patches are occupied and which are unoccupied.
  occ_patches <- dplyr::bind_rows(occ_list) |>
    dplyr::distinct(patch, threshold, .keep_all = TRUE) |>
    dplyr::arrange(threshold, patch)

  write.csv(
    occ_patches,
    file = file.path(PROJ_ROOT, 'results', 'patch_summaries', paste0(x[['version']][1], '-occupied.csv')),
    row.names = FALSE)

  # aggregate predicted probability per patch ('mean'/'median' per patch).
  central <- dplyr::bind_rows(central_list)

  write.csv(
    central,
    file = file.path(PROJ_ROOT, 'results', 'patch_summaries', paste0(x[['version']][1], '-cntrlTend.csv')),
    row.names = FALSE)

}



patchDist <- function(x, patch_lkp, r_name, thresh_type){
  # simplify the patches - we will extract only the borders of the patches as 
  # any migration event would have to pass through these. 
  patch_borders <- terra::as.polygons(x) |> 
    sf::st_as_sf() 
  sf::st_agr(patch_borders) = "constant"
  
  patch_borders <- patch_borders |>
    sf::st_cast('MULTILINESTRING') |>
    sf::st_cast('LINESTRING') |>
    dplyr::rename('Patch' = 1)
  
  # now convert these lines into segments less than or (equal to) 90m in length #
  segs <- stplanr::line_segment(patch_borders, segment_length = 180, use_rsgeo = TRUE) |>
    dplyr::group_by(Patch) |>
    dplyr::mutate(ID = 1:dplyr::n(), .before = geometry)
  
  # from each segment sample one point along the lines length - near the center. 
  # this gives us control on where we calculate distances to/from.
  sf::st_agr(segs) = "constant"
  pts <- sf::st_point_on_surface(segs)
  
  # if there are many many points we will sample them to limit them to 100 points per
  # large polygon. 
  pts <- pts |>
    dplyr::group_by(Patch) |>
    dplyr::slice_sample(n = 100) |>
    dplyr::ungroup()
  
  # we focus on calculating a subset of distances between patches known to 
  # be occupied and patches not known to be occupied. We will separate the data set 
  # into the occupied and unoccupied patches. 
  
  focal_pts <- dplyr::filter(pts, Patch %in% patch_lkp$patch)
  nf_pts <- dplyr::filter(pts, !Patch %in% patch_lkp$patch)
  
  # calculate all distances from each of these points to to our focal occupied polygons #
  
  d_mat <- sf::st_distance(focal_pts, nf_pts)
  d_mat <- data.frame(apply(d_mat, MARGIN = 2, as.numeric))
  d_mat <- apply(d_mat, MARGIN = 2, round, 0)
  
  colnames(d_mat) <- paste0(nf_pts$Patch, '_', nf_pts$ID)
  rownames(d_mat) <- paste0(focal_pts$Patch, '_', focal_pts$ID)
  
  d_mat <- data.frame(
    cbind(
      t(d_mat), 
      setNames(
        data.frame(sf::st_coordinates(nf_pts)),
        c('Longitude', 'Latitude')
      )
    )
  )
  
  # prepare data for analysis and to return. 
  nms <- c('Occ_patch_ID', 'Occ_node_ID')
  
  d_long <- d_mat |>
    
    # create a long data set 
    tibble::rownames_to_column('UnOcc_patch') |>
    tidyr::pivot_longer(cols = starts_with('X'), values_to = 'distance', names_to = 'Occ_patch') |>
    dplyr::relocate(Longitude:Latitude, .after = last_col()) |>
    
    # split apart the patch and node IDs; we only need a fraction of distances between the 
    # occupied and unoccupied patches to get a sense of their spatial distances. 
    dplyr::mutate(Occ_patch = stringr::str_remove(Occ_patch, 'X')) |>
    tidyr::separate_wider_delim(Occ_patch, names = nms, delim = '_') |>
    tidyr::separate_wider_delim(UnOcc_patch, names = paste0('Un', nms), delim = '_') |>
    dplyr::mutate(across(.cols = everything(), as.numeric)) |>
    
    # reduce the number of patch X patch distances to at most 5 pairs of measurements. 
    # these theoretically can cover up to 400 meters of edge between both patches. 
    
    dplyr::group_by(UnOcc_patch_ID, Occ_patch_ID) |>
    dplyr::slice_min(order_by = distance, n  = 5) |> 
    
    # and make it so results will be in human readable order
    
    dplyr::arrange(UnOcc_patch_ID, Occ_patch_ID)
  
  # on paper, these measurements may more accurately measure the distance between patches - that is
  # how far a seed would need to travel on a crows foot to get from one patch to the other
  
  # it is useful because most patches are seriously elongated, and getting from the 
  # centroid of one to the other may be enormous! while the margin(al) habitats may
  # actually be very close
  
  # now using these data, we can identify the single closest linkages in euclidean 
  # space between patches.
  min_dist <- d_long |>
    dplyr::summarise(Min_distance = min(distance), .groups = 'drop_last')
  
  mean_dist_all_patches <- d_long |>
    dplyr::ungroup(Occ_patch_ID) |>
    dplyr::summarise(
      Mean_min_distance = mean(distance), 
      Quant25 = quantile(distance, 0.25), 
      Quant75 = quantile(distance, 0.25), 
      .groups = 'drop_last',
    )
  
  fp <- file.path(PROJ_ROOT, 'results', 'patch_distances')

  r_name <- paste0(r_name, '-', thresh_type)
  write.csv(min_dist, file.path(fp, paste0(r_name, 'MinDist.csv')), row.names = F)
  write.csv(mean_dist_all_patches, file.path(fp, paste0(r_name, 'MeanDist.csv')), row.names = F)
  write.csv(d_long, file.path(fp, paste0(r_name, 'Distances.csv')), row.names = F)

}


#' Extract border sample points from occupied patches
#'
#' Refactored from the interior of \code{patchDist()} so that occupied-patch
#' source points can be obtained independently — e.g. as inputs to
#' \code{leastCostDist()} and \code{leastCostDistAniso()}.
#'
#' @param patches_lyr      Single-layer SpatRaster of patch IDs for one
#'   threshold type (e.g. one layer from the output of \code{patchAttributes()}).
#' @param occ_ids          Integer vector of occupied patch IDs, as read from
#'   the \code{-occupied.csv} table produced by \code{patchAttributes()}.
#' @param max_pts_per_patch Integer. Maximum border sample points retained per
#'   patch after segment sampling (default 100).
#'
#' @return sf POINT object. One point near the centre of each border segment
#'   of every occupied patch, capped at \code{max_pts_per_patch} per patch.
getOccupiedBorderPts <- function(patches_lyr, occ_ids, max_pts_per_patch = 100) {

  occ_rast <- terra::ifel(patches_lyr %in% occ_ids, patches_lyr, NA)

  borders <- terra::as.polygons(occ_rast) |>
    sf::st_as_sf() |>
    sf::st_cast('MULTILINESTRING') |>
    sf::st_cast('LINESTRING')

  colnames(borders)[1] <- 'Patch'
  sf::st_agr(borders) <- 'constant'

  segs <- stplanr::line_segment(borders, segment_length = 180, use_rsgeo = TRUE) |>
    dplyr::group_by(Patch) |>
    dplyr::mutate(ID = dplyr::row_number())

  sf::st_agr(segs) <- 'constant'

  sf::st_point_on_surface(segs) |>
    dplyr::group_by(Patch) |>
    dplyr::slice_sample(n = max_pts_per_patch) |>
    dplyr::ungroup()
}


#' Subset sf points falling within a terra extent
#'
#' @param pts sf POINT object.
#' @param e   terra SpatExtent.
#'
#' @return sf POINT object, subset of \code{pts}.
sourcesInExt <- function(pts, e) {
  xy <- sf::st_coordinates(pts)
  e  <- as.vector(e)
  pts[xy[, 1] >= e[1] & xy[, 1] <= e[2] & xy[, 2] >= e[3] & xy[, 2] <= e[4], ]
}


#' Split source points into geographically disjoint buffered regions, capped
#' by cell count
#'
#' Used to keep \code{gdistance::transition()} matrices under the ~2.1B
#' nonzero-entry cap of the \code{Matrix} package's sparse matrix
#' representation (32-bit index limit) — and, in practice, far under that:
#' at 1-3arc (~9m) resolution, \code{gdistance::transition()}'s R-level
#' matrix assembly was observed to hold flat around 12GB RSS for ~49 minutes
#' on a ~94M-cell region and then blow up to 115GB+ and crash R within
#' seconds, on a machine with 125GB RAM. A synthetic 3M-cell benchmark and a
#' real-data benchmark up to 1.75M cells both stayed fast (<1min) and safe
#' (<10GB), so cost-distance layers are built per population cluster (e.g.
#' CB vs. Cochetopa Dome, ~50km apart) and, if a cluster is still too big,
#' recursively bisected (kmeans k=2) until every leaf's cropped extent is
#' under \code{max_cells}. Leaves are then mosaicked back together.
#'
#' kmeans() has no default seed, so which population lands in "cluster 1"
#' vs "cluster 2" (and thus which one a given for-loop iteration hits first)
#' is otherwise non-deterministic between runs — set \code{seed} for
#' reproducible tiling.
#'
#' @param pts       sf POINT object of source points (e.g. occupied-patch
#'   border points from \code{getOccupiedBorderPts()}).
#' @param k         Integer. Number of initial geographic clusters (kmeans).
#' @param buffer_m  Numeric. Buffer distance in metres added around each
#'   cluster's bounding box.
#' @param max_cells Integer or NULL. Cap on a leaf cluster's cropped-extent
#'   cell count (bbox+buffer area / \code{res_m}^2). Clusters over the cap
#'   are recursively bisected. NULL disables the cap (legacy behaviour:
#'   exactly \code{k} clusters, no recursion).
#' @param res_m     Numeric. Raster resolution in metres/cell. Required if
#'   \code{max_cells} is set.
#' @param seed      Integer or NULL. Passed to \code{set.seed()} before each
#'   kmeans() call for reproducible splits. NULL leaves the RNG untouched.
#'
#' @return list of terra SpatExtent, one per (possibly recursively split)
#'   cluster.
clusterRegionExtents <- function(pts, k = 2, buffer_m = 3000,
                                  max_cells = NULL, res_m = NULL, seed = 1) {
  stopifnot('res_m is required when max_cells is set' =
              is.null(max_cells) || !is.null(res_m))

  ext_of <- function(p) {
    sf::st_bbox(p) |>
      sf::st_as_sfc(crs = sf::st_crs(p)) |>
      sf::st_buffer(buffer_m) |>
      terra::vect() |>
      terra::ext()
  }
  ncell_of <- function(e) {
    v <- as.vector(e)
    ((v[2] - v[1]) / res_m) * ((v[4] - v[3]) / res_m)
  }
  kmeans_cluster <- function(p, centers) {
    if (!is.null(seed)) set.seed(seed)
    kmeans(sf::st_coordinates(p), centers = centers)$cluster
  }

  split_recursive <- function(p) {
    e <- ext_of(p)
    if (is.null(max_cells) || nrow(p) < 4 || ncell_of(e) <= max_cells) {
      return(list(e))
    }
    cl <- kmeans_cluster(p, centers = 2)
    if (length(unique(cl)) < 2) return(list(e))   # degenerate split, give up
    do.call(c, lapply(sort(unique(cl)), function(i) split_recursive(p[cl == i, ])))
  }

  cl0 <- kmeans_cluster(pts, centers = k)
  do.call(c, lapply(sort(unique(cl0)), function(i) split_recursive(pts[cl0 == i, ])))
}


#' Assemble a bi-directional landscape resistance surface
#'
#' All layers are min-max normalised to [0, 1] before weighting so weights
#' express relative influence directly. Slope and TRI use magnitude only
#' (no aspect dependency), keeping the surface symmetric. Suitability is
#' inverted so that unsuitable cells resist movement. Forest and shrub
#' proportions enter as-is (already [0, 1]). Water bodies receive a fixed
#' additive penalty large enough to act as a strong but non-infinite barrier,
#' preserving graph connectivity.
#'
#' Concavity/convexity layers can be added later by passing them through
#' the \code{...} mechanism once their role is better understood.
#'
#' @param slope       SpatRaster or file path. Slope magnitude in degrees.
#' @param tri         SpatRaster or file path. Topographic Ruggedness Index.
#' @param suitability SpatRaster or file path. Habitat suitability [0, 1].
#' @param forest      SpatRaster or file path. Proportion canopy cover [0, 1].
#' @param shrub       SpatRaster or file path. Proportion shrub cover [0, 1].
#' @param water       SpatVector, SpatRaster, or file path. Water body features.
#'   Rasterised to the template grid; intersecting cells receive a fixed
#'   additive penalty controlled by \code{water_cost}.
#' @param weights     Named numeric vector with elements: slope, tri,
#'   suitability, forest, shrub. Relative weights applied after normalisation.
#' @param water_cost  Numeric. Additive penalty for water cells on the same
#'   scale as the summed weighted layers — 10 is roughly 3-4x the maximum
#'   non-water resistance.
#'
#' @return SpatRaster. Resistance surface ready for \code{terra::costDist()}.
buildCostSurface <- function(
    slope, tri, suitability, forest, shrub,
    water      = NULL,
    weights    = c(slope = 1, tri = 1, suitability = 1, forest = 1, shrub = 0.5),
    water_cost = 10
) {

  load_lyr <- function(x) if (is.character(x)) terra::rast(x) else x
  slope       <- load_lyr(slope)
  tri         <- load_lyr(tri)
  suitability <- load_lyr(suitability)
  forest      <- load_lyr(forest)
  shrub       <- load_lyr(shrub)

  norm01 <- function(r) {
    mn <- terra::global(r, 'min', na.rm = TRUE)[[1]]
    mx <- terra::global(r, 'max', na.rm = TRUE)[[1]]
    (r - mn) / (mx - mn)
  }

  resistance <-
    norm01(slope)           * weights[['slope']]       +
    norm01(tri)             * weights[['tri']]         +
    norm01(1 - suitability) * weights[['suitability']] +
    forest                  * weights[['forest']]      +
    shrub                   * weights[['shrub']]

  if (!is.null(water)) {
    water <- load_lyr(water)
    water_r <- if (inherits(water, 'SpatVector')) {
      terra::rasterize(water, resistance, field = 1, background = 0)
    } else {
      terra::resample(water, resistance, method = 'near')
    }
    resistance <- resistance + (water_r * water_cost)
  }

  names(resistance) <- 'resistance'
  resistance
}


#' Build an isotropic (symmetric) transition matrix for least-cost distance
#'
#' Symmetric counterpart to \code{buildAnisotropicTransition()}. Conductance
#' between adjacent cells is the inverse of the mean resistance of the two
#' cells, with no directional component. Keeps both iso and aniso cost
#' distances on the same gdistance/accCost code path.
#'
#' @param resistance SpatRaster from \code{buildCostSurface()}.
#' @param directions Integer. Cell connectivity: 8 (default) or 4.
#'
#' @return gdistance TransitionLayer. Pass to \code{leastCostDist()}.
buildIsotropicTransition <- function(resistance, directions = 8) {
  load_lyr <- function(x) if (is.character(x)) terra::rast(x) else x
  resistance <- load_lyr(resistance)

  tr <- gdistance::transition(
    raster::raster(resistance),
    transitionFunction = function(x) 1 / mean(x),
    directions = directions,
    symm = TRUE
  )
  gdistance::geoCorrection(tr, type = 'c', scl = FALSE)
}


#' Compute least-cost distances from occupied patch borders to all raster cells
#'
#' Wraps \code{gdistance::accCost()}, matching the interface of
#' \code{leastCostDistAniso()} so both iso and aniso distances follow the
#' same code path.
#'
#' Returns raw accumulated cost, NOT rescaled — when tiling by region
#' (CostDistances.R), each region's cost range differs, so rescaling per-tile
#' before mosaicking would put tiles on different scales. Rescale once, after
#' mosaicking all regions together, on the merged raster.
#'
#' \code{gdistance::transition()} drops the CRS itself (confirmed: even
#' \code{raster::crs(transition)} is already NA before \code{accCost()} ever
#' runs), so it can't be recovered from \code{transition} downstream — it
#' must be passed in from the original resistance/DEM raster and reapplied
#' here as metadata (extent/res/dims are untouched by the round-trip, so this
#' is not a resample).
#'
#' @param transition gdistance TransitionLayer from
#'   \code{buildIsotropicTransition()}.
#' @param sources    SpatVector, sf, or two-column matrix (x, y).
#' @param crs        Character or crs object. CRS of the raster the
#'   transition was built from (e.g. \code{terra::crs(resistance)}) —
#'   reapplied to the output since gdistance/raster round-trips lose it.
#'
#' @return SpatRaster. Raw accumulated least-cost distance from the nearest
#'   source (same units as the resistance surface); Inf cells set to NA.
leastCostDist <- function(transition, sources, crs) {
  if (inherits(sources, 'SpatVector')) {
    sources <- sf::st_as_sf(sources) |> sf::as_Spatial()
  } else if (inherits(sources, c('sf', 'sfc'))) {
    sources <- sf::as_Spatial(sources)
  }
  r <- terra::rast(gdistance::accCost(transition, sources))
  terra::crs(r) <- crs
  terra::ifel(is.infinite(r), NA, r)
}


#' Build an anisotropic (directional) transition matrix for least-cost distance
#'
#' terra::costDist() is isotropic — it cannot distinguish uphill from downhill
#' movement. This function uses gdistance to build a directed transition matrix
#' where resistance between any two adjacent cells depends on the direction of
#' traversal: going uphill costs exponentially more, going downhill costs
#' exponentially less, reflecting gravity-driven seed dispersal.
#'
#' The slope asymmetry replaces the isotropic slope term from
#' \code{buildCostSurface()}. To avoid double-counting, build the base surface
#' with \code{weights[['slope']] = 0} before passing it here.
#'
#' The combined transition multiplies the directed slope conductance by the
#' isotropic base conductance at each directed edge — equivalent to a
#' multiplicative interaction between terrain direction and landscape
#' permeability.
#'
#' Requires a projected CRS (metres). The existing UTM 32613 projection is fine.
#'
#' @param dem        SpatRaster or file path. Digital elevation model. Used
#'   only for the directional slope component.
#' @param resistance SpatRaster from \code{buildCostSurface()} with
#'   \code{weights[['slope']] = 0}. All non-slope resistance layers.
#' @param uphill_k   Numeric. Exponential penalty rate for uphill movement.
#'   Conductance scales as \code{exp(-uphill_k * rise_per_metre)}. Adapted
#'   from Tobler's hiking function (default 3.5); higher values steepen the
#'   uphill penalty.
#' @param downhill_k Numeric. Exponential facilitation rate for downhill
#'   movement. Conductance scales as \code{exp(downhill_k * drop_per_metre)}.
#'   Lower than uphill_k (default 1.5) reflects that gravity assistance is
#'   more modest than gravity resistance — steep descents also risk seeds
#'   overshooting suitable habitat.
#' @param directions Integer. Cell connectivity: 8 (default, includes
#'   diagonals) or 4.
#'
#' @return gdistance TransitionLayer. Pass to \code{leastCostDistAniso()}.
buildAnisotropicTransition <- function(
    dem,
    resistance,
    uphill_k   = 3.5,
    downhill_k = 1.5,
    directions = 8
) {

  load_lyr <- function(x) if (is.character(x)) terra::rast(x) else x
  dem        <- load_lyr(dem)
  resistance <- load_lyr(resistance)

  cell_m <- mean(terra::res(dem))

  # --- directed slope conductance ---
  # transitionFunction receives [elev_from, elev_to] for each directed edge.
  # symm = FALSE preserves directionality: A→B != B→A when there is relief.
  slope_trans <- gdistance::transition(
    raster::raster(dem),
    transitionFunction = function(x) {
      rise <- (x[2] - x[1]) / cell_m   # positive = uphill, negative = downhill
      ifelse(
        rise >= 0,
        exp(-uphill_k   *  rise),       # uphill:   conductance falls steeply
        exp( downhill_k * -rise)        # downhill: conductance rises moderately
      )
    },
    directions = directions,
    symm = FALSE
  )
  slope_trans <- gdistance::geoCorrection(slope_trans, type = 'c', scl = FALSE)

  # --- isotropic base conductance (all non-slope layers) ---
  # mean resistance of the two adjacent cells → conductance for that edge.
  base_trans <- gdistance::transition(
    raster::raster(resistance),
    transitionFunction = function(x) 1 / mean(x),
    directions = directions,
    symm = TRUE
  )
  base_trans <- gdistance::geoCorrection(base_trans, type = 'c', scl = FALSE)

  slope_trans * base_trans
}


#' Compute anisotropic least-cost distances from source points to all cells
#'
#' Wraps \code{gdistance::accCost()}, which returns the minimum accumulated
#' cost from any source point to each raster cell — the directed analogue of
#' \code{leastCostDist()}. Result is returned as a terra SpatRaster to stay
#' consistent with the rest of the pipeline.
#'
#' Returns raw accumulated cost, NOT rescaled — see \code{leastCostDist()}
#' for why (per-tile rescaling before mosaicking puts regions on different
#' scales; rescale once after mosaicking instead).
#'
#' Also see \code{leastCostDist()} for why \code{crs} must be passed in
#' explicitly: \code{gdistance::transition()} drops the CRS itself, so it
#' cannot be recovered from \code{transition} downstream.
#'
#' @param transition gdistance TransitionLayer from
#'   \code{buildAnisotropicTransition()}.
#' @param sources    SpatVector, sf, or two-column matrix (x, y) in the same
#'   CRS as the transition layer. Source points — typically border sample
#'   points of known-occupied patches.
#' @param crs        Character or crs object. CRS of the raster the
#'   transition was built from (e.g. \code{terra::crs(dem)}) — reapplied to
#'   the output since gdistance/raster round-trips lose it.
#'
#' @return SpatRaster. Raw accumulated anisotropic cost distance from the
#'   nearest source, in the same units as the resistance surface; Inf cells
#'   set to NA.
leastCostDistAniso <- function(transition, sources, crs) {

  if (inherits(sources, 'SpatVector')) {
    sources <- sf::st_as_sf(sources) |> sf::as_Spatial()
  } else if (inherits(sources, c('sf', 'sfc'))) {
    sources <- sf::as_Spatial(sources)
  }

  r <- terra::rast(gdistance::accCost(transition, sources))
  terra::crs(r) <- crs
  terra::ifel(is.infinite(r), NA, r)
}


#' Mosaic per-region raw cost-distance tiles into one globally-rescaled surface
#'
#' Companion to \code{leastCostDist()}/\code{leastCostDistAniso()}, which now
#' return raw (unscaled) accumulated cost precisely so this step can rescale
#' once, globally, after mosaicking — not per-tile before it. Adjacent
#' regions' buffered extents (see \code{clusterRegionExtents()}) can overlap;
#' those cells are resolved by minimum cost (\code{fun = 'min'}), i.e. a cell
#' reachable more cheaply from a neighbouring region's sources wins, rather
#' than an arbitrary first-tile-wins pick. Finally extended onto
#' \code{template}'s grid — tiles and template share origin/res (mosaicking
#' only shrinks the extent to the clustered regions' neighbourhoods), so this
#' is a pad with NA, not a resample.
#'
#' @param tiles    list of SpatRaster, raw accumulated cost, one per region,
#'   from \code{leastCostDist()}/\code{leastCostDistAniso()}.
#' @param template SpatRaster. Target grid to extend onto (e.g. the version's
#'   \code{-Pr.tif} suitability raster).
#'
#' @return SpatRaster, rescaled to [0, 65535] (ready for INT2U), on
#'   \code{template}'s grid.
mosaicCostTiles <- function(tiles, template) {
  merged <- if (length(tiles) > 1) {
    terra::mosaic(terra::sprc(tiles), fun = 'min')
  } else {
    tiles[[1]]
  }
  mn <- terra::global(merged, 'min', na.rm = TRUE)[[1]]
  mx <- terra::global(merged, 'max', na.rm = TRUE)[[1]]
  terra::round((merged - mn) / (mx - mn) * 65535L) |>
    terra::extend(template)
}
