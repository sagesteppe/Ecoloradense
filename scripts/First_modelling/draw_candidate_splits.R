## Draw candidate spatial splits for highPA_grid.R without fitting anything, and report
## which populations each puts in test - for resolutions where the default draw holds a
## whole population out (1/3 arc-second: seed 5010 put all 56 Cochetopa presences in
## test, none in train).
##
## Each candidate is saved as results_highPA_fieldabs/tables/{res}-split-seed{seed}.rds
## with its by-site table beside it. To use one, copy it over {res}-split.rds in that
## directory (highPA_grid.R reuses a cached split as long as it was built on the same
## points) and record the seed - the cached rds is the reproducible record of the draw.
##
## Usage (working directory scripts/):
##   Rscript First_modelling/draw_candidate_splits.R 10 6010 7010 8010

suppressPackageStartupMessages({
  library(sf)
  library(tidyverse)
  library(terra)
})
source('functions.R')

args <- commandArgs(trailingOnly = TRUE)
resolution <- as.numeric(args[1])
split_seeds <- as.integer(args[-1])
res_str <- res_string(resolution)
p2dat <- file.path('..', 'data', 'Data4modelling')
out_dir <- file.path(PROJ_ROOT, 'results_highPA_fieldabs', 'tables')

# the modelling points, built exactly as highPA_grid.R builds them (same .id order)
abs <- bind_rows(
  st_read(file.path(p2dat, 'iter1-pa.gpkg'), quiet = TRUE) |> select(Occurrence, ID),
  st_read(file.path(p2dat, 'iter1-pa-highPA-reviewed.gpkg'), quiet = TRUE) |> select(Occurrence, ID)
)
pres_file <- c('90' = '90m', '30' = '30m', '10' = '3m')[[as.character(resolution)]]
x <- st_read(file.path(p2dat, paste0(pres_file, '-presence-iter1.gpkg')), quiet = TRUE) |>
  rename(Occurrence = Presenc) |>
  st_as_sf() |>
  mutate(field_abs = Occurrence == 0)
x <- bind_rows(x, mutate(abs, field_abs = FALSE))
x <- filter(x, st_is(x, 'POINT')) |>
  select(Occurrence, field_abs) |>
  mutate(.id = seq_len(n()))

areas <- st_read(file.path('..', 'data', 'GroundTruthPts', 'ground_truth_median-VISITED.gpkg'),
                 layer = 'areas', quiet = TRUE) |>
  st_transform(st_crs(x))
x$Site <- areas$Site[st_nearest_feature(x, areas)]
rast_dat <- rastReader(paste0('dem_', res_str), p2proc = file.path('..', 'data', 'spatial', 'processed'))

for(s in split_seeds){
  set.seed(s)
  split <- spatial_class_split(select(x, Occurrence, .id), rast_dat, train_split = 0.8, test_tolerance = 0.1)
  saveRDS(split, file.path(out_dir, sprintf('%s-split-seed%d.rds', res_str, s)))

  set_of <- ifelse(x$.id %in% split$test_id, 'test', 'train')
  by_site <- x |>
    st_drop_geometry() |>
    mutate(set = set_of,
           class = case_when(Occurrence == 1 ~ 'presence', field_abs ~ 'field_abs', TRUE ~ 'background')) |>
    filter(class != 'background') |>
    count(Site, class, set) |>
    pivot_wider(names_from = c(class, set), values_from = n, values_fill = 0)
  write.csv(by_site, file.path(out_dir, sprintf('%s-split-seed%d-by-site.csv', res_str, s)), row.names = FALSE)

  held_out <- by_site$Site[by_site$presence_test >= 5 & by_site$presence_train == 0]
  message(sprintf('\n=== split seed %d: %d test presences; whole population(s) held out: %s',
                  s, sum(by_site$presence_test), if(length(held_out)) paste(held_out, collapse = ', ') else 'none'))
  print(as.data.frame(by_site))
}
