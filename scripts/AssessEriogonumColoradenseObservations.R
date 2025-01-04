library(tidyverse)
library(sf)

setwd('~/Documents/Ecoloradense/scripts')

inat <- read.csv('../data/collections/iNaturalistObservations_coloradense/observations-417874.csv') %>% 
  filter(positional_accuracy <= 10) %>% 
  st_as_sf(coords = c('longitude', 'latitude'), crs = 4326) %>% 
  select(id, observed_on, user_login, quality_grade) 

soro <- readr::read_csv('../data/collections/SoroSymbiota_coloradense/occurrences.csv', show_col_types = FALSE) %>% 
  drop_na(decimalLatitude, decimalLongitude) %>% 
  st_as_sf(coords = c('decimalLongitude', 'decimalLatitude'), crs = 4326) %>% 
  select(georeferenceRemarks, georeferencedBy,locality, eventDate, recordedBy, 
         identifiedBy, identificationRemarks, catalogNumber)

st_write(inat, '../data/collections/manual_cleaning_geodata/inat.shp', append = F)
st_write(soro, '../data/collections/manual_cleaning_geodata/soro.shp', append = F)

Soro_pts_remove <- c('WSC0003417', '109996', 'KHD00057673', '48427', 'RMBL0004279',
                     '102217', '3076163', '3076165', '3076162', '719290', 'KHD00074171',
                     'RMBL004279', 'ASU0121771', '3076166', '48360', 'KHD00039107')

WSC0003417 # this point point on a west aspect slope near veg - discard
109996 # this point put in irwin graveyard
KHD00057673 # not yule pass. 
48427 # at kettle ponds, not emerald lake
102217 # judd falls, not avery - paradise
KHD00039107 # likely along trail just before copper lake
3076163 & 3076162 # this is meant to be cocheotopa area
3076165 # no info
RMBL004279 # should be further up slope
ASU0121771 # not enough email
3076166 # not enough info, area unlikely
48360 # wrong position and covered by rick
719290 # should be up near peak, down on 401 instead
RMBL0004279 # close but likely up the ridge a smidge
KHD00074171 # Mis-identified?

soro <- filter(soro, ! catalogNumber %in% Soro_pts_remove)

# combine the data sets to provide a unified data set. 

bind_rows(
  select(soro, eventDate) ,
  select(inat, eventDate = observed_on) %>% 
    mutate(eventDate = as.Date(eventDate))
) %>% 
  st_write('../data/collections/occurrences/occurrences.shp', append = F)

rm(inat, soro, Soro_pts_remove)

