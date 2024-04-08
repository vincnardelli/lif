library(dplyr)
library(sf)
library(lubridate)
library(spdep)
library(eurostat)
library(tmap)

map <- get_eurostat_geospatial(
  resolution = "60",
  nuts_level = "2",
  year = 2016
) %>% 
  filter(CNTR_CODE !="IS") %>% 
  st_crop(xmin=-12, ymin=30, xmax=45, ymax=80 ) 

nb <- poly2nb(map, queen=TRUE)
map <- map %>% 
  mutate(neigh = card(nb)) %>% 
  filter(neigh > 0)
nb <- poly2nb(map, queen=TRUE)
listw<-nb2listw(nb, zero.policy = F)

tm_shape(map) +
  tm_polygons() +
  tm_grid()


data_raw <- eurostat::get_eurostat_json(
  id = "hlth_cd_ysdr2"
) %>% 
  filter(icd10=="C50")


data <- map %>%
  left_join(data_raw, by = "geo") %>% 
  na.omit() %>% 
  select(time, values, geometry)

save(data, file="data/eurostat_hlth_cd_ysdr2_c50.Rdata")

