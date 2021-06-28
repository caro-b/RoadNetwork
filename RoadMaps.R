
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
gfc_02 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_02.tif")
gfc_03 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_03.tif")
gfc_04 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_04.tif")
gfc_05 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_05.tif")
gfc_06 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_06.tif")
gfc_07 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_07.tif")
gfc_08 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_08.tif")
gfc_09 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_09.tif")
gfc_10 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_10.tif")
gfc_11 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_11.tif")
gfc_12 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_12.tif")
gfc_13 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_13.tif")
gfc_14 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_14.tif")
gfc_15 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_15.tif")
gfc_16 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_16.tif")
gfc_17 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_17.tif")
gfc_18 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_18.tif")
gfc_19 <- raster("C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/forest_loss_19.tif")
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

# function for function for masking loss areas in gfc loss data according to road data
maskForestLoss <- function(forest_loss, roads_union) {
  # masking forest loss areas in gfc loss data according to road data
  forest_loss[forest_loss < 1] <- NA
  forest_loss_roads <- mask(forest_loss, roads_union)
  # raster to polygons
  forest_loss_roads_pol <- rasterToPolygons(forest_loss_roads, dissolve=T)
  return(forest_loss_roads_pol)
}


## year 2001
# mask loss areas in gfc loss data according to road data
gfc_01_roads_pol <- maskForestLoss(gfc_01, roads_union)

# roads which intersect with forest loss
intersect_01 <- roads_union[intersect(roads_union, gfc_01_roads_pol)]

# remove intersecting roads from road network for analysis of next year
roads_union_01 <- roads_union[is.na(over(roads_union, gfc_01_roads_pol))]


## year 2002
gfc_02_roads_pol <- maskForestLoss(gfc_02, roads_union)

intersect_02 <- roads_union_01[intersect(roads_union_01, gfc_02_roads_pol)]
roads_union_02 <- roads_union_01[is.na(over(roads_union_01, gfc_02_roads_pol))]


## year 2003
gfc_03_roads_pol <- maskForestLoss(gfc_03, roads_union)

intersect_03 <- roads_union_02[intersect(roads_union_02, gfc_03_roads_pol)]
roads_union_03 <- roads_union_02[is.na(over(roads_union_02, gfc_03_roads_pol))]


## year 2004
gfc_04_roads_pol <- maskForestLoss(gfc_04, roads_union)

intersect_04 <- roads_union_03[intersect(roads_union_03, gfc_04_roads_pol)]
roads_union_04 <- roads_union_03[is.na(over(roads_union_03, gfc_04_roads_pol))]


## year 2005
gfc_05_roads_pol <- maskForestLoss(gfc_05, roads_union)

intersect_05 <- roads_union_04[intersect(roads_union_04, gfc_05_roads_pol)]
roads_union_05 <- roads_union_04[is.na(over(roads_union_04, gfc_05_roads_pol))]


## year 2006
gfc_06_roads_pol <- maskForestLoss(gfc_06, roads_union)

intersect_06 <- roads_union_05[intersect(roads_union_05, gfc_06_roads_pol)] # no intersecting roads
roads_union_06 <- roads_union_05[is.na(over(roads_union_05, gfc_06_roads_pol))]


## year 2007
gfc_07_roads_pol <- maskForestLoss(gfc_07, roads_union)

intersect_07 <- roads_union_06[intersect(roads_union_06, gfc_07_roads_pol)]
roads_union_07 <- roads_union_06[is.na(over(roads_union_06, gfc_07_roads_pol))]


## year 2008
gfc_08_roads_pol <- maskForestLoss(gfc_08, roads_union)

intersect_08 <- roads_union_07[intersect(roads_union_07, gfc_08_roads_pol)]
roads_union_08 <- roads_union_07[is.na(over(roads_union_07, gfc_08_roads_pol))]


## year 2009
gfc_09_roads_pol <- maskForestLoss(gfc_09, roads_union)

intersect_09 <- roads_union_08[intersect(roads_union_08, gfc_09_roads_pol)]
roads_union_09 <- roads_union_08[is.na(over(roads_union_08, gfc_09_roads_pol))]


## year 2010
gfc_10_roads_pol <- maskForestLoss(gfc_10, roads_union)

intersect_10 <- roads_union_09[intersect(roads_union_09, gfc_10_roads_pol)] # no intersecting roads
roads_union_10 <- roads_union_09[is.na(over(roads_union_09, gfc_10_roads_pol))]


## year 2011
gfc_11_roads_pol <- maskForestLoss(gfc_11, roads_union)

intersect_11 <- roads_union_10[intersect(roads_union_10, gfc_11_roads_pol)]
roads_union_11 <- roads_union_10[is.na(over(roads_union_10, gfc_11_roads_pol))]


## year 2012
gfc_12_roads_pol <- maskForestLoss(gfc_12, roads_union)

intersect_12 <- roads_union_11[intersect(roads_union_11, gfc_12_roads_pol)]
roads_union_12 <- roads_union_11[is.na(over(roads_union_11, gfc_12_roads_pol))]


## year 2013
gfc_13_roads_pol <- maskForestLoss(gfc_13, roads_union)

intersect_13 <- roads_union_12[intersect(roads_union_12, gfc_13_roads_pol)]
roads_union_13 <- roads_union_12[is.na(over(roads_union_12, gfc_13_roads_pol))]


## year 2014
gfc_14_roads_pol <- maskForestLoss(gfc_14, roads_union)

intersect_14 <- roads_union_13[intersect(roads_union_13, gfc_14_roads_pol)]
roads_union_14 <- roads_union_13[is.na(over(roads_union_13, gfc_14_roads_pol))]


## year 2015
gfc_15_roads_pol <- maskForestLoss(gfc_15, roads_union)

intersect_15 <- roads_union_14[intersect(roads_union_14, gfc_15_roads_pol)]
roads_union_15 <- roads_union_14[is.na(over(roads_union_14, gfc_15_roads_pol))]


## year 2016
gfc_16_roads_pol <- maskForestLoss(gfc_16, roads_union)

intersect_16 <- roads_union_15[intersect(roads_union_15, gfc_16_roads_pol)]
roads_union_16 <- roads_union_15[is.na(over(roads_union_15, gfc_16_roads_pol))]


## year 2017
gfc_17_roads_pol <- maskForestLoss(gfc_17, roads_union)

intersect_17 <- roads_union_16[intersect(roads_union_16, gfc_17_roads_pol)] # no intersecting roads
roads_union_17 <- roads_union_16[is.na(over(roads_union_16, gfc_17_roads_pol))]


## year 2018
gfc_18_roads_pol <- maskForestLoss(gfc_18, roads_union)

intersect_18 <- roads_union_17[intersect(roads_union_17, gfc_18_roads_pol)]
roads_union_18 <- roads_union_17[is.na(over(roads_union_17, gfc_18_roads_pol))]


## year 2019
gfc_19_roads_pol <- maskForestLoss(gfc_19, roads_union)

intersect_19 <- roads_union_18[intersect(roads_union_18, gfc_19_roads_pol)]
roads_union_19 <- roads_union_18[is.na(over(roads_union_18, gfc_19_roads_pol))]


## year 2020
gfc_20_roads_pol <- maskForestLoss(gfc_20, roads_union)

intersect_20 <- roads_union_19[intersect(roads_union_19, gfc_20_roads_pol)] # no intersecting roads
roads_union_20 <- roads_union_19[is.na(over(roads_union_19, gfc_20_roads_pol))]



plot(aoi)
plot(roads_union, add = T, col ="pink")
plot(gfc_01_roads, add=T, col = "purple")



#### VISUALIZATION ####

# plot road development per year
plot(aoi, main= "Road Development per year")
plot(intersect_01, add=T, col = "purple", main = "year 2001")
plot(intersect_02, add=T, col = "brown", main = "year 2002")
plot(intersect_03, add=T, col = "red", main = "year 2003")
plot(intersect_04, add=T, col = "yellow", main = "year 2004")
plot(intersect_05, add=T, col = "pink", main = "year 2005")

plot(intersect_10, add=T, col = "blue", main = "year 2010")
plot(intersect_15, add=T, col = "orange", main = "year 2015")
plot(intersect_20, add=T, col = "green", main = "year 2020")
legend("topleft",legend=c("Year 2001","Year 2002","Year 2003","Year 2004","Year 2005","Year 2010","Year 2015","Year 2020"),
       col=c("purple","brown","red","yellow","pink","blue","orange","green"),lty = 1)
