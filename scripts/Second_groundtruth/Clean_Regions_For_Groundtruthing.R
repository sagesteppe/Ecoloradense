library(dplyr)
source('functions.R')

sf::sf_use_s2(FALSE)
region_v <- sf::st_read('../data/hikingTrails/Feasible_to_Groundtruth.gpkg', quiet = TRUE) |>
  sf::st_make_valid() |>
  sf::st_transform(32613) |>
  sf::st_make_valid() |>
  dplyr::group_by(Site) |>
  dplyr::summarise(geometry = sf::st_union(geometry), .groups = 'drop')
sf::sf_use_s2(TRUE)

trails = sf::st_read('../data/hikingTrails/hiking_trails_2026_all.gpkg')|>
  sf::st_make_valid() |>
  sf::st_transform(32613) |>
  sf::st_make_valid() |>
  dplyr::filter(eligible == 'TRUE') |>
  dplyr::group_by(Trail) |>
  dplyr::summarise(geometry = sf::st_union(geom), .groups = 'drop') |>
  sf::st_buffer(25) |>
  dplyr::rename(Site = Trail) |>
  mutate(Site = stringr::str_replace(Site, 'Yule$', 'Yule Pass'))

obbies = bind_rows(region_v, trails)
unioned = obbies |>
  group_by(Site) |>
  dplyr::summarise(geometry = sf::st_union(geometry), .groups = 'drop')  |>
  sf::st_make_valid()

unioned$Site


## remove non public access areas

domain <- sf::st_read('../data/collections/occurrences_coloradense/occurrences.shp', quiet = T) %>% 
  st_union() %>% 
  st_transform(32613) %>% 
  st_buffer(16093) %>% 
  st_transform(4326) %>% 
  st_bbox()

pub <- st_read('../data/spatial/PADUS4_0_State_CO_GDB/PADUS4_0_StateCO.gdb',
               layer = 'PADUS4_0Fee_State_CO', quiet = T) |>
  ensure_multipolygons() |> 
  st_make_valid() |>
  st_transform(4326) |>
  st_crop(domain) |>
  st_transform(32613) |>
  dplyr::summarise(geometry = sf::st_union(geom))


unioned <- sf::st_intersection(unioned, pub)

sf::st_write(
  unioned, 
  file.path('..', 'data', 'hikingTrails', 'GroundTruch_Areas.gpkg'),
  append = FALSE
)
