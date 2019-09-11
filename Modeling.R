#Mapping the creosote bush
#https://www.gbif.org/species/7453543

library (spocc)
library (mapr)
library (dismo)
library(ENMeval)
library(raster)
library(rgdal)
library(rgbif)

creosote_df = occ(query = "Larrea tridentata", 
               from = "gbif",
               limit = 1000,
               has_coords = TRUE)

map_leaflet(creosote_df)

#Filter out outlier locations
clim = raster::stack(list.files('/usr/share/data/wc2.0/bio2.5/', pattern='.tif', full.names = TRUE))
#crop climate data down to north america
#plot data points on a climate map