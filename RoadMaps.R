
#### SETUP ####

## Packages
# install required packages if needed
packagelist <- c("dplyr","gganimate","ggplot2","ggthemes","maps","osmdata","raster","rgee","rgdal","rgeos","sf","tidyverse")
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
#dir <- list.files(path="C:/Users/carob/Dropbox/EAGLE/SS21/Conservation", pattern="*.tif")
#for (i in 1:length(dir)) assign(dir[i], raster(dir[i]))


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



#### VISUALIZATION ####

# plot road development per time period of 5 years
plot(aoi, main= "Road Development per year")
plot(intersect_01, add=T, col = "magenta", main = "year 2001")
plot(intersect_02, add=T, col = "magenta", main = "year 2002")
plot(intersect_03, add=T, col = "magenta", main = "year 2003")
plot(intersect_04, add=T, col = "magenta", main = "year 2004")

plot(intersect_05, add=T, col = "yellow", main = "year 2005")

plot(intersect_07, add=T, col = "yellow", main = "year 2007")
plot(intersect_08, add=T, col = "yellow", main = "year 2008")
plot(intersect_09, add=T, col = "yellow", main = "year 2009")

plot(intersect_11, add=T, col = "blue", main = "year 2011")
plot(intersect_12, add=T, col = "blue", main = "year 2012")
plot(intersect_13, add=T, col = "blue", main = "year 2013")
plot(intersect_14, add=T, col = "blue", main = "year 2014")

plot(intersect_15, add=T, col = "orange", main = "year 2015")
plot(intersect_16, add=T, col = "orange", main = "year 2016")

plot(intersect_18, add=T, col = "orange", main = "year 2018")
plot(intersect_19, add=T, col = "orange", main = "year 2019")

legend("topleft",legend=c("Year 2001-2004","Year 2005-2009","Year 2011-2014","Year 2015-2019"),
       col=c("magenta","yellow","blue","orange"),lty = 1)


## Animation

# first convert data to sf & add year column
#intersect_01_sf <- st_as_sf(intersect_01)
intersect_01_sf$year <- 2001
 
# # convert to SpatialLinesDataFrame & add year
# intersect_01_sl <- SpatialLinesDataFrame(intersect_01, as.data.frame(intersect_01_sf$year), match.ID = F)


intersect_02_sf <- st_as_sf(intersect_02)
intersect_02_sf$year <- 2002

intersect_03_sf <- st_as_sf(intersect_03)
intersect_03_sf$year <- 2003

intersect_04_sf <- st_as_sf(intersect_04)
intersect_04_sf$year <- 2004

intersect_05_sf <- st_as_sf(intersect_05)
intersect_05_sf$year <- 2005

intersect_07_sf <- st_as_sf(intersect_07)
intersect_07_sf$year <- 2007

intersect_08_sf <- st_as_sf(intersect_08)
intersect_08_sf$year <- 2008

intersect_09_sf <- st_as_sf(intersect_09)
intersect_09_sf$year <- 2009

intersect_11_sf <- st_as_sf(intersect_11)
intersect_11_sf$year <- 2011

intersect_12_sf <- st_as_sf(intersect_12)
intersect_12_sf$year <- 2012

intersect_13_sf <- st_as_sf(intersect_13)
intersect_13_sf$year <- 2013

intersect_14_sf <- st_as_sf(intersect_14)
intersect_14_sf$year <- 2014

intersect_15_sf <- st_as_sf(intersect_15)
intersect_15_sf$year <- 2015

intersect_16_sf <- st_as_sf(intersect_16)
intersect_16_sf$year <- 2016

intersect_18_sf <- st_as_sf(intersect_18)
intersect_18_sf$year <- 2018

intersect_19_sf <- st_as_sf(intersect_19)
intersect_19_sf$year <- 2019



# combine years into one dataframe
intersect <- rbind(intersect_01_sf, intersect_02_sf, intersect_03_sf, intersect_04_sf, intersect_05_sf, intersect_07_sf, intersect_08_sf, intersect_09_sf,
                    intersect_11_sf, intersect_12_sf, intersect_13_sf, intersect_14_sf, intersect_15_sf, intersect_16_sf, intersect_18_sf, intersect_19_sf)



# basic plot
area <- 
  ggplot(data=intersect) +
  geom_sf(aes(fill = as.factor(year))) +
  #scale_fill_hue() +
  borders(aoi, colour = "gray85") +
  theme_map() 

plot(area)



# animation - road development per year
map <- area +
  theme(legend.position = "none") +
  # Here comes the gganimate part
  labs(title = "Road Development 2000-2020", subtitle = "Year: {closest_state}") +
  transition_states(year) +
  ease_aes('linear') 

# rather consecutive road development of previous years?

