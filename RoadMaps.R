
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

# set directory to folder where input files lie
dir <- "C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/"

# input area of interest
aoi <- readOGR(paste(dir, "prediction_area.shp", sep = ""))

# input manually digitized roads from IZW
roads_2014 <- readOGR(paste(dir, "2014.shp", sep = ""))


## input OSM - roads data

# convert bb to an overpass query object (API)
roads_osm <- opq(aoi@bbox) %>%
  add_osm_feature(key = "highway") %>%  
  # output as simple features object (or R spatial (sp) - osmdata_sp()) - advantage of sf: use geom_sf() for ggplot2
  osmdata_sp()

# extract roads as lines
roads_osm_lines <- roads_osm$osm_lines


## input Hansen GFC raster files 

# read all files from directory into list
files <- list.files(path= dir, pattern="*.tif", full.names=T)

# stack raster images
gfc_stack <- stack(files)

# convert to raster brick to shorten processing time
gfc <- brick(gfc_stack)




#### PRE-PROCESSING ####

crs(gfc$forest_loss_01)
# CRS("+proj=longlat +datum=WGS84")

# set crs of aoi
crs(aoi) <- crs(gfc$forest_loss_01)

# transform crs of roads
roads_osm_lines <- spTransform(roads_osm_lines, CRS("+proj=longlat +datum=WGS84"))
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

# combine roads from OSM & IZW

# remove roads where OSM & IZW roads intersect
# take OSM data (as better for scaling)
# (only include additional roads from IZW which aren't in OSM data)
roads_2014_new <- roads_2014_clip[is.na(over(roads_2014_clip,roads_osm_clip))]

# only use osm roads, add additional IZW roads
roads_union <- rbind(roads_osm_clip, roads_2014_new)



## buffer around roads (as OSM only provides lines, not the original street width)
# highway width = 22m
# crs needs to be in long/lat for buffer in meters
roads_union@proj4string # WGS84 (EPSG: 4326)

# set crs to EPSG:3405 VN-2000 / UTM zone 48N
# https://spatialreference.org/ref/epsg/vn-2000-utm-zone-48n/
# Proj4js.defs["EPSG:3405"] = "+proj=utm +zone=48 +ellps=WGS84 +units=m +no_defs";
#roads_union <- spTransform(roads_union, CRS("+proj=utm +zone=48 +ellps=WGS84 +units=m +no_defs"))

roads_union_buffer <- buffer(roads_union, width = 0.0005, dissolve=F) #~50m buffer

plot(roads_union)
plot(roads_union_buffer, add=T)

#### question for Martin: why does area(gfc_01_roads_pol) correspond to pixel size (900m²), but width in buffer() only works with degrees?? ####
## buffer() - Unit of width is meter if x has a longitude/latitude CRS


#### ANALYSIS ####

## function for masking forest loss areas in gfc loss data according to road data
maskForestLoss <- function(forest_loss, roads) {
  # masking forest loss areas in gfc loss data according to road data
  forest_loss[forest_loss < 1] <- NA
  forest_loss_roads <- mask(forest_loss, roads)
  forest_loss_roads_pol <- rasterToPolygons(forest_loss_roads, dissolve=F)
  return(forest_loss_roads_pol)
}


## iterate over all raster layers (=years): 
# 1. mask forest loss according to roads 
# 2. polygonize raster pixels
# 3. save polygons into one spatialpolygonsdataframe & add year (as numeric)

# index for iteration: starts at year 2001 (no forest loss in year 2000) = second layer in raster
d <- 1

# create empty spatialpolygons object
pol<- SpatialPolygons(list(), proj4string=crs(gfc$forest_loss_01))
dat <- data.frame(matrix(ncol=1,nrow=0, dimnames=list(NULL, c("year"))))
sp <- SpatialPolygonsDataFrame(pol, data= dat)

while(d < (nlayers(gfc))) {  
  # mask forest loss according to roads & polygonize raster pixels via function maskForestLoss
  # save as spatialpolygons (remove data section which would lead to errors in rbind due to different column naming)
  sp_new <- SpatialPolygons(maskForestLoss(gfc[[d+1]], roads_union_buffer)@polygons, proj4string=crs(gfc$forest_loss_01))
  # add year as column
  sp_new$year <- rep(d, length(sp_new))
  # save all polygons into one spatialpolygonsdataframe (via rbind())
  sp <- rbind(sp, sp_new)
  d=d+1
}

table(sp$year)
plot(sp[sp$year == 1,])

# TODO: convert numeric values to actual years?? ####


# 4. stepwise (cumulative) intersection of forest loss pixels & roads

# via loop
test <- bufferRoads(sp, roads_union_buffer)



## function to buffer around each road, only when x% of forest loss pixels in buffer, then road dev
bufferRoads <- function(gfc_roads_pol, roads_buffer) {
  # calculate area of buffers
  roads_buffer$area <- area(roads_buffer)
  # calculate area of forest loss polygons
  gfc_roads_pol$area <- area(gfc_roads_pol)
  # count number of & calculate sum of area of forest loss pixels which lie in each buffer
  pixelsInBuffer <- over(roads_buffer, gfc_roads_pol, fn = sum)
  # get roads for which x% of forest loss pixels lie in buffer, then road dev
  roads <- roads_buffer[((pixelsInBuffer$area / area(roads_buffer)) >= 0.3) & !is.na(pixelsInBuffer$area),]
  return(roads)
}

# ## function to get corresp. forest loss pixels which lie on roads
# bufferPixels <- function(gfc_roads_pol, roads) {
#   # check if roads is empty
#   if (length(roads) == 0) {
#     # bullshit comparison to pass empty spatialpolygonsdataframe
#     pixels <- gfc_roads_pol[gfc_roads_pol$area == "hello",]
#   } else {
#     # forest loss pixels which lie on roads
#     pixels <- gfc_roads_pol[!is.na(over(gfc_roads_pol, roads)[,1]),]
#   }
#   return(pixels)
# }





## year 2001

gfc_01_roads_pol <- maskForestLoss(gfc$forest_loss_01, roads_union_buffer)

plot(gfc_01_roads_pol)

# roads which correspond to road dev. in year 2001
roads_01 <- bufferRoads(gfc_01_roads_pol, roads_union_buffer)

# forest loss pixels which corresponds to road dev. in year 2001
gfc_roadDev_01 <- bufferPixels(gfc_01_roads_pol, roads_01)


#### TODO: better & consistent naming ! ####


## year 2002
# remove intersecting roads from previous year
if (length(roads_01) != 0) {
  roads_01_new <- erase(roads_union_buffer, roads_01)
} else {
  roads_01_new <- roads_union_buffer
}

gfc_02_roads_pol <- maskForestLoss(gfc$forest_loss_02, roads_union_buffer)

# combine pixels with non-road dev. pixels from previous year
rbind(gfc_02_roads_pol, gfc_01_roads_pol)


gfc_02_roads_pol_new <- rbind(pixelsNotIn_01, gfc_02_roads_pol, make.row.names=F)


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

