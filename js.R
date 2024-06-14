library(dplyr)
library(sf)
library(ggplot2)
library(patchwork)
library(lubridate)
library(spdep)
library(tmap)
library(robspdep)
library(readxl)
library(spData)

map <- read_sf("data/john_snow/polys.shp")
map$indicatore <- map$Death_dens
pumps <- read_sf("data/john_snow/Pumps.shp")
streets <- read_sf("data/john_snow/streets_js.shp")


centroids <- st_centroid(map)
knn <- knn2nb(knearneigh(st_coordinates(centroids), k = 5)) # Adjust 'k' as needed
listw <- nb2listw(knn)

# Your existing code
t <- tm_shape(map) +
  tm_fill("Death_dens",
          style = "quantile",
          n = 7, 
          palette = "Reds", 
          title = "Death rate") +
  tm_borders(col = "black", lwd = 0.2) +
  tm_shape(streets) + 
  tm_lines(col = "black", lwd = 1) +  # Add streets as black lines
  tm_shape(pumps) +
  tm_symbols(col = "black", size = 1) +  # Add pumps as points, modify color and size as needed
  tm_layout(frame = FALSE,    # Add your map title here
            title.position = c("center", "top"),
            legend.position = c("right", "bottom"),  # Position of the legend, can be "left", "right", "top", "bottom"
            legend.outside = TRUE)         # Place the legend outside of the map

# Print or save the map
t
tmap_save(t, "results/js_cholera_map.png")



a <- ggplot(map) +
  geom_density(aes(Death_dens)) +
  theme_minimal() 
b <- ggqqplot(map$Death_dens)
a | b
ggsave("results/js_cholera_density.png", width=1280, height=720, units = "px", scale = 2, bg="white")

x <- map$indicatore
m <- moran.mc(x, listw, 999)
m$statistic
m$p.value


mh <- gk(x, listw)
mh$statistic
mh$p.value