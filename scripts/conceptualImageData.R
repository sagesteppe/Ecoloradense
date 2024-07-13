library(terra)
library(sf)
library(tidyverse)
library(fields)
set.seed(20)
setwd('~/Documents/Ecoloradense/scripts')

# A conceptual figure should showcase the interrelated goals of modelling 
# 1) Predicting suitable habitat 
# 2) Showcasing the distance between Occupied and unoccupied patches 
# 3) Emphasizing the size and shapes of patches 
# 4) Displaying the population borders 
# 5) Predicting plant density per cell  

# we will create around 5 population clumps across the example domain. 
data <- data.frame(
  x = c( 
    rnorm(n = 10, mean = 20, sd = 5), 
    rnorm(n = 10, mean = 40, sd = 10), 
    rnorm(n = 10, mean = 100, sd = 7), 
    rnorm(n =  5, mean = 45, sd = 15), # very disjunct
    rnorm(n = 15, mean = 165, sd = 20), # large population
    rnorm(n = 25, mean = 225, sd = 10) # big population. 
  ), 
  y = c(
    rnorm(n = 10, mean = 85, sd = 5), 
    rnorm(n = 10, mean = 40, sd = 7), 
    rnorm(n =  5, mean = 100, sd = 8), 
    rnorm(n = 10, mean = 225, sd = 10), # very disjunct
    rnorm(n = 15, mean = 65, sd = 10), # large population
    rnorm(n = 25, mean = 175, sd = 15) # big population . 
  ), 
  pr_suit = c(
    rnorm(n = 10, mean = 95, sd = 10), 
    rnorm(n = 10, mean = 85, sd = 5), 
    rnorm(n =  5, mean = 90, sd = 10), 
    rnorm(n = 10, mean = 80, sd = 5), 
    rnorm(n = 15, mean = 95, sd = 5), 
    rnorm(n = 25, mean = 95, sd = 5)
  ), 
  pr_dens = c(
    rnorm(n = 10, mean = 20, sd = 3), 
    rnorm(n = 10, mean = 15, sd = 2), 
    rnorm(n =  5, mean = 10, sd = 3), 
    rnorm(n = 10, mean = 5, sd = 5), 
    rnorm(n = 15, mean = 25, sd = 5), 
    rnorm(n = 25, mean = 25, sd = 5)
  ), 
  patch = c(
    rep(LETTERS[1], times = 10),
    rep(LETTERS[2], times = 10),
    rep(LETTERS[3], times = 5),
    rep(LETTERS[4], times = 10),
    rep(LETTERS[5], times = 15),
    rep(LETTERS[6], times = 25)
  )
) |> 
 mutate(
    pr_suit = if_else(pr_suit > 100, 100, pr_suit)) |>
 st_as_sf(coords = c('x', 'y'), crs = 32613)

# we will create an empty template raster
r <- rast(ext(data), nrow = 250, ncols = 250, crs = 'EPSG:32613')

# we'll create background points- those with low probability of suitable habitat
# and low plant density estimates. 

bg <- st_bbox(data) |>
  st_sample(type = 'regular', size = 100) |>
  st_as_sf() |> 
  mutate(
    pr_suit = sample(1:45, size = n(), replace = TRUE),
    pr_dens = sample(0:2, size = n(), replace = TRUE), 
    Patch = NA, .before = x) |> 
  rename(geometry = x) |> 
  st_as_sf() 

ch <- data |> 
  st_buffer(3) |> 
  group_by(patch) |> 
  summarise(geometry = st_union(geometry)) |> 
  st_convex_hull()

bg <- bg[ lengths( st_disjoint(bg, st_union(ch)) ) > 0, ]
rm(ch)
  
# use kriging to simply interpolate these values across the domain
data_sp <- bind_rows(data, bg)

r_pr_suit <- interpolate(r, Tps(st_coordinates(data_sp), data_sp$pr_suit)) 
r_pr_suit <- ifel(r_pr_suit > 100, 100, r_pr_suit) 
r_pr_suit <- ifel(r_pr_suit < 0, 0, r_pr_suit) 

# we will smooth out the initial points which went into the kriging... 
r_pr_suit_l <- focal(
  ifel(r_pr_suit < 60, r_pr_suit, NA), w = 7, fun = 'min', na.rm=TRUE)

r_pr_suit_h <- focal(
  ifel(r_pr_suit >= 60, r_pr_suit, NA), w = 7, fun = 'mean', na.rm=TRUE) 

r_pr_suit <- mean(r_pr_suit_l, r_pr_suit_h, na.rm = T)

rm(r_pr_suit_l, r_pr_suit_h)
# now repeat the process for the plant density
r_pr_dens <- interpolate(r, Tps(st_coordinates(data_sp), data_sp$pr_dens)) 

# combine into the same raster 
r <- c(r_pr_suit, r_pr_dens) 
names(r) <- c('pr_suit', 'pr_dens') 

rm(bg, data_sp, r_pr_dens, r_pr_suit) 

# Create the example patches for plotting
r <- c(r, patches(ifel(r$pr_suit < 60, NA, 1))) 

# we will define that some of these areas are occupied, and some of them lack plants #
occ <- data.frame( 
  patch = LETTERS[1:6], 
  Occupied = c(0, 1, 0, 0, 1, 1) 
) 

data <- left_join(data, occ) 
patch_occupancy <- bind_cols(data, 
  terra::extract(r$patches, vect(data), ID = FALSE) 
) |> 
  distinct(patch, Occupied, patches) |> 
  drop_na() 

# detect the 'boundaries' of the populations 
boundaries <- boundaries(r$patches) 
boundaries <- ifel(boundaries == 0, NA, 1) |> 
  as.polygons() |> 
  st_as_sf() 

rm(data)
# 6 rows in one column, with text to the side. 
writeRaster(r, '../data/ConceptualFigure/DensityPatchSuit.tif', overwrite = T)
st_write(boundaries, '../data/ConceptualFigure/boundaries.shp', append = F)
write.csv(patch_occupancy, row.names = F, '../data/ConceptualFigure/occupancy.csv')
