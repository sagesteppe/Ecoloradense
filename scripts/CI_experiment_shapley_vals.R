library(xgboost)
library(iml)
library(tidyverse)

setwd('/home/sagesteppe/Documents/Ecoloradense/scripts')
mod <- readRDS('../results/CountModels/models/1-3arc-Iteration1-PA1_2.7-poisson_spat.rds')

test_indices <- as.numeric(
  read.delim(
    '../results/CountModels/test_data/twin_indx-1-3arc.txt', header = F, sep = ' ')
  )

train_dat <- sf::st_read('../data/Data4modelling/10m-count-iter1.gpkg')[test_indices,] |>
  rowwise() |>
  mutate(Prsnc_All = sum(Prsnc_M, Prsnc_J)) |>
  select(Prsnc_All)

#######
# this is where we should extract values from a raster to points
#######


## for practice code we have this in here instead

# these are old and from something else. so the results won't be very sensible but should allow for testing
tdat <- read.csv('../results/test_data/1-3arc-Iteration1-PA1_2.7DO_0.csv')


predictor <- Predictor$new(mod, data = tdat, y = tdat$Prsnc_All)


rm(test_indices, rast_vals)
# Calculate SHAP values for a single observation 
shapley <- Shapley$new(predictor, x.interest = tdat[1, ]) 










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
