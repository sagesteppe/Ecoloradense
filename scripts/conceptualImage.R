library(terra)

set.seed(20)

# A conceptual figure should showcase the interrelated goals of modelling
# 1) Predicting suitable habitat
# 2) Noting the Occupied habitat patches
# 3) Showcasing the distance between Occupied and unoccupied patches
# 4) Emphasizing the size and shapes of patches
# 5) Displaying the population borders
# 6) Predicting plant density per cell. 

# it is likely to require 100 cells or so.  
r <- rast(nrow = 250, ncols = 250)

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
 st_as_sf(coords = c('x', 'y'), crs = 32613)

# we'll create background points- those with low probability of suitable habitat
# and low plant density estimates. 

bg <- st_bbox(data) |>
  st_sample(type = 'regular', size = 100) |>
  st_as_sf() |> 
  mutate(
    pr_suit = sample(1:65, size = n(), replace = TRUE),
    pr_dens = sample(0:2, size = n(), replace = TRUE), .before = x) |>
  rename(geometry = x) |>
  st_as_sf()

ch <- data |>
  st_buffer(3) |>
  group_by(patch) |>
  summarise(geometry = st_union(geometry)) |>
  st_convex_hull()

bg <- bg[ lengths( st_disjoint(bg, st_union(ch)) ) > 0, ]
rm(ch)
  
library(gstat)

# generate the background data set. 
x <- seq(1, 250, length.out = 250)
y <- seq(1, 250, length.out = 250)
grd <- expand.grid(x = x, y = y) |>
  st_as_sf(coords = c('x', 'y'), crs = 32613) |>
  as_Spatial()

# ensure that it doesn't overlap with the existing points

dt.vgm <- variogram(pr_suit ~ 1, data)
dt.fit <- fit.variogram(dt.vgm, model = vgm(1, "Lin", 10, 1)) # fit model
lzn.kriged <- krige((pr_suit) ~ 1, data, grd, model = dt.fit) |>
  as.data.frame()

r <- setValues(r, lzn.kriged$var1.pred)
plot(r)


