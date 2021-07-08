
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
# highway expected to be opened in 2021
# crs needs to be in long/lat for buffer in meters
roads_union@proj4string # WGS84 (EPSG: 4326)

# set crs to EPSG:3405 VN-2000 / UTM zone 48N
# https://spatialreference.org/ref/epsg/vn-2000-utm-zone-48n/
# Proj4js.defs["EPSG:3405"] = "+proj=utm +zone=48 +ellps=WGS84 +units=m +no_defs";
#roads_union <- spTransform(roads_union, CRS("+proj=utm +zone=48 +ellps=WGS84 +units=m +no_defs"))

roads_union_buffer <- buffer(roads_union, width = 0.0001, dissolve=F) #~10m buffer: 0.0001

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
d <- as.numeric(1)

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


## function to buffer around each road, only when x% of forest loss pixels in buffer, then road dev
bufferRoads <- function(gfc_roads_pol, roads_buffer) {
  # calculate area of buffers
  roads_buffer$area <- area(roads_buffer)
  # calculate area of forest loss polygons
  gfc_roads_pol$area <- area(gfc_roads_pol)
  # count number of & calculate sum of area of forest loss pixels which lie in each buffer
  pixelsInBuffer <- over(roads_buffer, gfc_roads_pol, fn = sum)
  # get roads for which x% of forest loss pixels lie in buffer, then road dev
  roads <- roads_buffer[((pixelsInBuffer$area / area(roads_buffer)) >= 0.05) & !is.na(pixelsInBuffer$area),]
  return(roads)
}


# 4. stepwise (cumulative) intersection of forest loss pixels & roads to find year of road development
# a) intersect forest loss pixels of certain year with all roads
# b) if >x% pixels intersecting with road - save year in road spatialpolygonsdataframe
# c) else cumulatively add forest loss pixels to next year & repeat from a)

# index for years
y <- 1
# list for adding cumulative years
year_list <- list(y)

# add year variable to road polygons
roads_union_buffer$year <- rep(NA, length(roads_union_buffer))

while(y < nlayers(gfc)) {
  polygons_current_year <- sp[sp$year %in% year_list,]
  # buffer around each road, only when x% of forest loss pixels in buffer, then road dev
  roads <- bufferRoads(polygons_current_year, roads_union_buffer)
  # if no road dev for that year, add year to list 
  if (length(roads) == 0) {
    year_list <- c(year_list, y)
  } else {
  # check if not already year saved (to not overwrite earlier year)
  if (is.na(mean(roads_union_buffer$year[roads_union_buffer@polygons %in% roads@polygons]))) {
    # note year of intersecting roads
    roads_union_buffer$year[roads_union_buffer@polygons %in% roads@polygons] <- y
  }
  # reset list to current year
  year_list <- list(y)
  }
  y=y+1
}

# drop roads with na year
roads_clean <- roads_union_buffer[!is.na(roads_union_buffer$year),]

lis <- seq(1:15)
plot(sp[sp$year %in% lis,])
#plot(roads,add=T)
plot(roads_clean[roads_clean$year==16,], add=T, col="magenta")

# roads not intersecting with forest loss pixels
roads_na <- roads_union_buffer[is.na(roads_union_buffer$year),]
plot(sp)
plot(roads_na, col="green")


#### TODO: better & consistent naming ! ####

# change numbers to actual years
roads_clean$year <- as.numeric(roads_clean$year) #### TODO: directly specify as double when adding new column ####
roads_clean$year[roads_clean$year <= 9] <- paste("200", roads_clean$year[roads_clean$year <= 9], sep="")
roads_clean$year[!grepl("200", roads_clean$year)] <- paste("20", roads_clean$year[!grepl("200", roads_clean$year)], sep="")


#### VISUALIZATION ####

# first convert data to sf
roads_clean_sf <- st_as_sf(roads_clean)

# convert year to factor (for plotting)
roads_clean_sf$year <- as.factor(roads_clean_sf$year)

# basic plot
area <- 
  ggplot(data=roads_clean_sf) +
  geom_sf(aes(fill = year)) +
  #scale_fill_hue() +
  borders(aoi, colour = "gray85") +
  theme_map() 

plot(area)

## Animation

# animation - road development per year
area +
  theme(legend.position = "none") +
  # Here comes the gganimate part
  transition_states(year) +
  transition_manual(year, cumulative = T) +
  labs(title = "Road Development 2000-2020", subtitle = roads_clean_sf$year) + #"Year: {closest_state}"
  ease_aes('linear')

## rather consecutive road development of previous years

