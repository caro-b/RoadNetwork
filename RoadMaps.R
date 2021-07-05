
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


## input OSM - roads data

# convert bb to an overpass query object (API)
roads_osm <- opq(aoi@bbox) %>%
  add_osm_feature(key = "highway") %>%  
  # output as simple features object (or R spatial (sp) - osmdata_sp()) - advantage of sf: use geom_sf() for ggplot2
  osmdata_sp()

# extract roads as lines
roads_osm_lines <- roads_osm$osm_lines


## input Hansen GFC raster files 

# set directory where downloaded files from GEE are in
dir <- "C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/"

# read all files from directory into list
files <- list.files(path= dir, pattern="*.tif", full.names=T)

# stack raster images
gfc_stack <- stack(files)

# convert to raster brick to shorten processing time
gfc <- brick(gfc_stack)




#### PRE-PROCESSING ####

# set crs of aoi
crs(aoi) <- crs(gfc$forest_loss_01)
crs(roads_osm_lines) <- crs(gfc$forest_loss_01)
crs(roads_2014) <- crs(gfc$forest_loss_01)

# clip roads to aoi 
# byid = T to keep single features (roads)
roads_osm_clip <- gIntersection(roads_osm_lines, aoi, byid = T)
roads_2014_clip <- gIntersection(roads_2014, aoi, byid = T)


# discrepancies between OSM & IZW road data
plot(aoi)
plot(roads_2014_clip, add = T, col ="blue")
plot(roads_osm_clip, add = T, col ="orange")
legend("topleft",legend=c("Roads IZW (2014)","Roads OSM (2020)"),
       col=c("blue","orange"),lty = 1)
# take OSM data (as better for scaling)


# combine roads from OSM & IZW

# remove roads where OSM & IZW roads intersect (only include additional roads from IZW which aren't in OSM data)
roads_2014_new <- roads_2014_clip[is.na(over(roads_2014_clip,roads_osm_clip))]

# only use osm roads, add additional IZW roads
roads_union <- rbind(roads_osm_clip,roads_2014_new)



## buffer around roads (as OSM only provides lines, not the original street width)
# highway width = 22m
# crs needs to be in long/lat for buffer in meters
roads_union@proj4string #WGS84 (EPSG: 4326)

# set crs to EPSG:3405 VN-2000 / UTM zone 48N
# https://spatialreference.org/ref/epsg/vn-2000-utm-zone-48n/
crs(roads_union) <- CRS("+proj=utm +zone=48 +ellps=WGS84 +units=m +no_defs")

roads_union_buffer <- buffer(roads_union, width = 0.0005, dissolve=F) #~50m buffer

plot(roads_union)
plot(roads_union_buffer, add=T)



#### ANALYSIS ####

## function for masking forest loss areas in gfc loss data according to road data
maskForestLoss <- function(forest_loss, roads) {
  # masking forest loss areas in gfc loss data according to road data
  forest_loss[forest_loss < 1] <- NA
  forest_loss_roads <- mask(forest_loss, roads)
  # raster to polygons
  forest_loss_roads_pol <- rasterToPolygons(forest_loss_roads, dissolve=F)
  return(forest_loss_roads_pol)
}


#### question for Martin: why does area(gfc_01_roads_pol) correspond to pixel size (900m²), but width oly works with degrees?? ####

# function to make buffer around forest loss pixels to check if forest loss pixels actually correspond to road development
bufferPixels <- function(gfc_roads_pol) {
  # buffering for SpatialPolygons can only deal with planar coordinate reference systems (eg UTM)
  gfc_buffer <- buffer(gfc_roads_pol, width = 0.001, dissolve=F) # ~100m buffer around forest loss pixels
  # check if more than x forest loss polygons in buffer - to only include forest loss polygons which correspond to road development
  # as each polygon corresponds to one pixel (and they have the similar size), possible to sum up pixels
  # count number of polygons in per buffer
  polygonsInBuffer <- over(gfc_roads_pol, gfc_buffer, fn = sum)
  #sort(test$forest_loss_01)
  # remove polygons which don't have 3 or more forest loss pixels per buffer 
  # (1 pixel = ~30m x 30m = 900m² forest loss), (3 pixels =~ 1800m² = 1.8km)
  # buffer (100m) / 3 * 30m
  roadDevPixels <- gfc_roads_pol[polygonsInBuffer$forest_loss_01 > 2,]
  return(roadDevPixels)
}
#### TODO: altern. with area instead of sum ? ####


## year 2001
# mask loss areas in gfc loss data according to road data (using function)
gfc_01_roads_pol <- maskForestLoss(gfc$forest_loss_01, roads_union_buffer)

plot(gfc_01_roads_pol)
plot(gfc_01_buffer, add = T)

# function to buffer around forest loss pixels to get forest loss pixels corresponding to road development
gfc_01_roadDev <- bufferPixels(gfc_01_roads_pol)

# roads which intersect with forest loss
intersect_01 <- roads_union_buffer[intersect(roads_union_buffer, gfc_01_roadDev)]

# remove intersecting roads from road network for analysis of next year
roads_union_01 <- roads_union_buffer[is.na(over(roads_union_buffer, gfc_01_roadDev))]


#### TODO: subdivide roads into ~km (?) segments??? ####
# length of spatial lines
sort(gLength(roads_union, byid = T))



## year 2002
gfc_02_roads_pol <- maskForestLoss(gfc$forest_loss_02, roads_union)

intersect_02 <- roads_union_01[intersect(roads_union_01, gfc_02_roads_pol)]
roads_union_02 <- roads_union_01[is.na(over(roads_union_01, gfc_02_roads_pol))]


## year 2003
gfc_03_roads_pol <- maskForestLoss(gfc$forest_loss_03, roads_union)

intersect_03 <- roads_union_02[intersect(roads_union_02, gfc_03_roads_pol)]
roads_union_03 <- roads_union_02[is.na(over(roads_union_02, gfc_03_roads_pol))]


## year 2004
gfc_04_roads_pol <- maskForestLoss(gfc$forest_loss_04, roads_union)

intersect_04 <- roads_union_03[intersect(roads_union_03, gfc_04_roads_pol)]
roads_union_04 <- roads_union_03[is.na(over(roads_union_03, gfc_04_roads_pol))]


## year 2005
gfc_05_roads_pol <- maskForestLoss(gfc$forest_loss_05, roads_union)

intersect_05 <- roads_union_04[intersect(roads_union_04, gfc_05_roads_pol)]
roads_union_05 <- roads_union_04[is.na(over(roads_union_04, gfc_05_roads_pol))]


## year 2006
gfc_06_roads_pol <- maskForestLoss(gfc$forest_loss_06, roads_union)

intersect_06 <- roads_union_05[intersect(roads_union_05, gfc_06_roads_pol)] # no intersecting roads
roads_union_06 <- roads_union_05[is.na(over(roads_union_05, gfc_06_roads_pol))]


## year 2007
gfc_07_roads_pol <- maskForestLoss(gfc$forest_loss_07, roads_union)

intersect_07 <- roads_union_06[intersect(roads_union_06, gfc_07_roads_pol)]
roads_union_07 <- roads_union_06[is.na(over(roads_union_06, gfc_07_roads_pol))]


## year 2008
gfc_08_roads_pol <- maskForestLoss(gfc$forest_loss_08, roads_union)

intersect_08 <- roads_union_07[intersect(roads_union_07, gfc_08_roads_pol)]
roads_union_08 <- roads_union_07[is.na(over(roads_union_07, gfc_08_roads_pol))]


## year 2009
gfc_09_roads_pol <- maskForestLoss(gfc$forest_loss_09, roads_union)

intersect_09 <- roads_union_08[intersect(roads_union_08, gfc_09_roads_pol)]
roads_union_09 <- roads_union_08[is.na(over(roads_union_08, gfc_09_roads_pol))]


## year 2010
gfc_10_roads_pol <- maskForestLoss(gfc$forest_loss_10, roads_union)

intersect_10 <- roads_union_09[intersect(roads_union_09, gfc_10_roads_pol)] # no intersecting roads
roads_union_10 <- roads_union_09[is.na(over(roads_union_09, gfc_10_roads_pol))]


## year 2011
gfc_11_roads_pol <- maskForestLoss(gfc$forest_loss_11, roads_union)

intersect_11 <- roads_union_10[intersect(roads_union_10, gfc_11_roads_pol)]
roads_union_11 <- roads_union_10[is.na(over(roads_union_10, gfc_11_roads_pol))]


## year 2012
gfc_12_roads_pol <- maskForestLoss(gfc$forest_loss_12, roads_union)

intersect_12 <- roads_union_11[intersect(roads_union_11, gfc_12_roads_pol)]
roads_union_12 <- roads_union_11[is.na(over(roads_union_11, gfc_12_roads_pol))]


## year 2013
gfc_13_roads_pol <- maskForestLoss(gfc$forest_loss_13, roads_union)

intersect_13 <- roads_union_12[intersect(roads_union_12, gfc_13_roads_pol)]
roads_union_13 <- roads_union_12[is.na(over(roads_union_12, gfc_13_roads_pol))]


## year 2014
gfc_14_roads_pol <- maskForestLoss(gfc$forest_loss_14, roads_union)

intersect_14 <- roads_union_13[intersect(roads_union_13, gfc_14_roads_pol)]
roads_union_14 <- roads_union_13[is.na(over(roads_union_13, gfc_14_roads_pol))]


## year 2015
gfc_15_roads_pol <- maskForestLoss(gfc$forest_loss_15, roads_union)

intersect_15 <- roads_union_14[intersect(roads_union_14, gfc_15_roads_pol)]
roads_union_15 <- roads_union_14[is.na(over(roads_union_14, gfc_15_roads_pol))]


## year 2016
gfc_16_roads_pol <- maskForestLoss(gfc$forest_loss_16, roads_union)

intersect_16 <- roads_union_15[intersect(roads_union_15, gfc_16_roads_pol)]
roads_union_16 <- roads_union_15[is.na(over(roads_union_15, gfc_16_roads_pol))]


## year 2017
gfc_17_roads_pol <- maskForestLoss(gfc$forest_loss_17, roads_union)

intersect_17 <- roads_union_16[intersect(roads_union_16, gfc_17_roads_pol)] # no intersecting roads
roads_union_17 <- roads_union_16[is.na(over(roads_union_16, gfc_17_roads_pol))]


## year 2018
gfc_18_roads_pol <- maskForestLoss(gfc$forest_loss_18, roads_union)

intersect_18 <- roads_union_17[intersect(roads_union_17, gfc_18_roads_pol)]
roads_union_18 <- roads_union_17[is.na(over(roads_union_17, gfc_18_roads_pol))]


## year 2019
gfc_19_roads_pol <- maskForestLoss(gfc$forest_loss_19, roads_union)

intersect_19 <- roads_union_18[intersect(roads_union_18, gfc_19_roads_pol)]
roads_union_19 <- roads_union_18[is.na(over(roads_union_18, gfc_19_roads_pol))]


## year 2020
gfc_20_roads_pol <- maskForestLoss(gfc$forest_loss_20, roads_union)

intersect_20 <- roads_union_19[intersect(roads_union_19, gfc_20_roads_pol)] # no intersecting roads
roads_union_20 <- roads_union_19[is.na(over(roads_union_19, gfc_20_roads_pol))]


# roads not intersecting with forest loss pixels
plot(aoi, main= "Roads not intersecting with Forest Loss Pixels")
plot(roads_union_20, add=T, col = "magenta")



#### VISUALIZATION ####

# plot road development per time period of 5 years
plot(aoi, main= "Road Development per Time Period")
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
intersect_01_sf <- st_as_sf(intersect_01)
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

# convert year to factor (for plotting)
intersect$year <- as.factor(intersect$year)


# basic plot
area <- 
  ggplot(data=intersect) +
  geom_sf(aes(fill = year)) +
  #scale_fill_hue() +
  borders(aoi, colour = "gray85") +
  theme_map() 

plot(area)


# animation - road development per year
  area +
  theme(legend.position = "none") +
  # Here comes the gganimate part
  labs(title = "Road Development 2000-2020", subtitle = "Year: {closest_state}") +
  transition_states(year) +
  ease_aes('linear') 

## rather consecutive road development of previous years

