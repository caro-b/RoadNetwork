
#### SETUP ####

## Packages
# install required packages if needed
packagelist <- c("dplyr","osmdata","raster","rgee","rgdal","rgeos","sf","tidyverse")
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

gfc_01 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_01.tif")
gfc_05 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_05.tif")
gfc_10 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_10.tif")
gfc_15 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_15.tif")
gfc_20 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_20.tif")
# dir <- list.files(path="C:/Users/carob/Dropbox/EAGLE/SS21/Conservation", pattern="*.tif")
# for (i in 1:length(dir)) assign(dir[i], raster(dir[i]))


#### ANALYSIS ####

# clip roads to aoi 
# byid = T to keep single features (roads)
roads_osm_clip <- gIntersection(roads_osm_lines, aoi, byid = T)
roads_2014_clip <- gIntersection(roads_2014, aoi, byid = T)

plot(aoi)
plot(roads_osm_clip, add = T, col ="pink")
plot(roads_2014_clip, add = T, col ="blue")

# 2 more roads in OSM, one more for IZW data

# union roads from OSM & IZW
roads_union <- union(roads_osm_clip, roads_2014_clip)
crs(roads_union) <- crs(gfc_01)

## mask loss areas in gfc data

# function for masking loss areas in gfc loss data according to road data
maskForestLoss <- function(forest_loss, roads_union) {
  forest_loss[forest_loss < 1] <- NA
  forest_loss_roads <- mask(forest_loss, roads_union)
  return(forest_loss_roads)
}

# year 2001
gfc_01_roads <- maskForestLoss(gfc_01, roads_union)
plot(roads_union, col ="pink")
plot(gfc_01_roads, col = "purple", main = "Forest Loss 2001", add=T)

# year 2005
gfc_05_roads <- maskForestLoss(gfc_05, roads_union)
plot(gfc_05_roads, col = "yellow", main = "Forest Loss 2005")

# year 2010
gfc_10_roads <- maskForestLoss(gfc_10, roads_union)
plot(gfc_10_roads, col = "blue", main = "Forest Loss 2010")

# year 2015
gfc_15_roads <- maskForestLoss(gfc_15, roads_union)
plot(gfc_15_roads, col = "purple", main = "Forest Loss 2015")

# year 2020
gfc_20_roads <- maskForestLoss(gfc_20, roads_union)
plot(gfc_20_roads, col = "pink", main = "Forest Loss 2020")

plot(aoi)
plot(roads_union, add = T, col ="pink")
plot(gfc_01_roads, add=T, col = "purple")


## extract roads on which forest loss pixels lie
# raster to polygons
# year 2001
gfc_01_roads_pol <- rasterToPolygons(gfc_01_roads, dissolve=T)
# year 2005
gfc_05_roads_pol <- rasterToPolygons(gfc_05_roads, dissolve=T)
# year 2010
gfc_10_roads_pol <- rasterToPolygons(gfc_10_roads, dissolve=T)
# year 2015
gfc_15_roads_pol <- rasterToPolygons(gfc_15_roads, dissolve=T)
# year 2020
gfc_20_roads_pol <- rasterToPolygons(gfc_20_roads, dissolve=T)



## year 2001
# roads which intersect with forest loss
intersect_01 <- roads_union[intersect(roads_union, gfc_01_roads_pol)]

# remove intersecting roads from road network for analysis of next year
roads_union_01 <- roads_union[is.na(over(roads_union, gfc_01_roads_pol))]


## year 2005
intersect_05 <- roads_union_01[intersect(roads_union_01, gfc_05_roads_pol)]
roads_union_05 <- roads_union_01[is.na(over(roads_union_01, gfc_05_roads_pol))]


## year 2010
intersect_10 <- roads_union_05[intersect(roads_union_05, gfc_10_roads_pol)]
roads_union_10 <- roads_union_05[is.na(over(roads_union_05, gfc_10_roads_pol))]


## year 2015
intersect_15 <- roads_union_10[intersect(roads_union_10, gfc_15_roads_pol)]
roads_union_15 <- roads_union_10[is.na(over(roads_union_10, gfc_15_roads_pol))]


## year 2020
intersect_20 <- roads_union_15[intersect(roads_union_15, gfc_20_roads_pol)]
roads_union_20 <- roads_union_15[is.na(over(roads_union_15, gfc_20_roads_pol))]



#### VISUALIZATION ####

# plot road development per year
plot(aoi, main= "Road Development per year")
plot(intersect_01, add=T, col = "purple", main = "year 2001")
plot(intersect_05, add=T, col = "pink", main = "year 2005")
plot(intersect_10, add=T, col = "blue", main = "year 2010")
plot(intersect_15, add=T, col = "orange", main = "year 2015")
plot(intersect_20, add=T, col = "green", main = "year 2020")
legend("topleft",legend=c("Year 2001", "Year 2005","Year 2010","Year 2015","Year 2020"),
       col=c("purple","pink","blue","orange","green"),lty = 1)
