library(xgboost)
library(iml)
library(tidyverse)

setwd('/home/sagesteppe/Documents/Ecoloradense/scripts')
#setwd('/media/steppe/hdd/EriogonumColoradenseTaxonomy/scripts')
mod <- readRDS('../results/CountModels/models/1-3arc-Iteration1-PA1_2.7-poisson_spat.rds')

test_indices <- as.numeric(
  read.delim(
    '../results/CountModels/test_data/twin_indx-1-3arc.txt', header = F, sep = ' ')
  )

test_dat <- sf::st_read('../data/Data4modelling/10m-count-iter1.gpkg')[test_indices,] |>
  rowwise() |>
  mutate(Prsnc_All = sum(Prsnc_M, Prsnc_J)) |>
  select(Prsnc_All)

#######
# this is where we should extract values from a raster to points
#######


## for practice code we have this in here instead


rm(test_indices, rast_vals)
# Calculate SHAP values for a single observation 

shapley <- Shapley$new(predictor, x.interest = tdat[1, ]) 

mod <- readRDS('../results/CountModels/models/1-3arc-Iteration1-PA1:2.7-tweedie.rds')$Model
vars <- as.character(attr(mod[["preproc"]][["terms"]], 'variables')) # fns will run on all 
vars <- vars[2:length(vars)] # vars even if not in model - but produce null wts as expected. 

test_indices <- as.numeric(
  read.delim(
    '../results/CountModels/test_data/twin_indx-1-3arc.txt', header = F, sep = ' ')
)

test_dat <- d_w_vars[test_indices,] |>
  select(any_of(vars)) |>
  st_drop_geometry()

predict.fun <- function(object, newdata){
  newData_x = xgb.DMatrix(data.matrix(newdata), missing = NA)
  results<-predict(mod, newData_x)
}

subby <- test[which(names(test) != "Prsnc_All")]
predictor <- Predictor$new(
  model = mod,
  data = subby, 
  y = test_dat$Prsnc_All, 
  predict.function = predict.fun
)

# this does not work... 
shapley <- Shapley$new(predictor, x.interest = test) 

private$sampler$feature.names

# these both run   
int_out <- Interaction$new(predictor, grid.size = 10)
plot(int_out)

prsuit_int <- Interaction$new(predictor, feature = "Pr.SuitHab", grid.size = 15)
plot(prsuit_int)

######## try to get shap with a different package. 
library(kernelshap)
library(shapviz)

s <- kernelshap::kernelshap(mod, X = test_dat) 
sv <- shapviz::shapviz(s)
shapviz::sv_importance(sv, kind = "bee") +
  theme_bw()

shapviz::sv_importance(sv, kind = "bar") +
  theme_bw()

shapviz::sv_dependence(sv, v = colnames(test_dat[, -1]))

shapviz::sv_waterfall(sv) +
  theme(axis.text = element_text(size = 11))









#####################################3 original tutorial

data("Boston", package = "MASS")

library("iml")
library("randomForest")
data("Boston", package = "MASS")
rf <- randomForest(medv ~ ., data = Boston, ntree = 10)


X <- Boston[which(names(Boston) != "medv")]
predictor <- Predictor$new(rf, data = X, y = Boston$medv)

shapley <- Shapley$new(predictor, x.interest = X[1, ], sample.size = 50)
shapley$plot()

interact <- Interaction$new(predictor, grid.size = 15)
plot(interact)

interacts_with <- Interaction$new(predictor, feature = "lstat", grid.size = 15)
plot(interacts_with)
