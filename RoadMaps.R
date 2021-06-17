
#### SETUP ####

## Packages
# install required packages if needed
packagelist <- c("rgdal", "osmdata","rgeos","dplyr","rgee","raster")
new.packages <- packagelist[!(packagelist %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load required packages
lapply(packagelist, require, character.only = TRUE)


#### DATA IMPORT ####
getwd()

# manually digitized roads from IZW
roads_2014 <- readOGR("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/2014.shp")
aoi <- readOGR("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/prediction_area.shp")


## OSM - roads data

# convert bb to an overpass query object (API)
roads_osm <- opq(aoi@bbox) %>%
  add_osm_feature(key = "highway") %>%  
  # output as simple features object (or R spatial (sp) - osmdata_sp()) - advantage of sf: use geom_sf() for ggplot2
  osmdata_sp()

# extract roads as lines
roads_osm_lines <- roads_osm$osm_lines


## GEE - Hansen Global Forest Change

# needs prior python installation
# # It is necessary just once
# ee_install()
# 
# # Initialize Earth Engine!
# ee_Initialize()

# aoi_sf <- st_as_sf(aoi)
# aoi_sf <- st_set_crs(aoi_sf, 4326)
# 
# ee_aoi <- aoi_sf %>% 
#  # st_read(quiet = TRUE) %>% 
#   sf_as_ee()
# region <- mask$geometry()$bounds()
# 
# ee_gfc <- ee$
#   Image('UMD/hansen/global_forest_change_2020_v1_8')$
#   clip(ee_aoi)
# 
# gfc <- ee_as_raster(ee_gfc, maxPixels=200000000)
#plot(gfc)

gfc_00 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_00.tif")
gfc_05 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_05.tif")
gfc_10 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_10.tif")
gfc_15 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_15.tif")
gfc_20 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_20.tif")
# dir <- list.files(path="C:/Users/carob/Dropbox/EAGLE/SS21/Conservation", pattern="*.tif")
# for (i in 1:length(dir)) assign(dir[i], raster(dir[i]))


#### PRE-PROCESSING ####

# clip roads to aoi 
# byid = T to keep single features (roads)
roads_osm_clip <- gIntersection(roads_osm_lines, aoi, byid = T)
roads_2014_clip <- gIntersection(roads_2014, aoi, byid = T)

plot(aoi)
plot(roads_osm_clip, add = T, col ="pink")
plot(roads_2014_clip, add = T, col ="blue")

# 2 more roads in OSM, one more for IZW data


# mask loss areas in gfc
hist(gfc_20)
gfc_20_loss <- gfc_20
gfc_20_loss[gfc_20_loss < 1] <- NA
plot(gfc_20_loss)


# union roads from OSM & IZW
roads_union <- union(roads_osm_clip, roads_2014_clip)

plot(aoi)
plot(roads_union, add = T, col ="pink")


# mask gfc according to roads
gfc_20_roads <- mask(gfc_20_loss, roads_union)

plot(gfc_20_roads, col = "purple")


