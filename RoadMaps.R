## Purpose of this script:
# Detect road development in Bach Ma Nationalpark (Vietnam) for the years 2000 untill 2020
# OSM Road data and Hansen Forest Loss data is intersected to detect the year in which the development of each road started
## Data sources: Hansen Global Forest Change Data v1.8, OSM Road Data, Sentinel2 (ESA) Satellite Data
## Author: Caroline Busse
## email: caroline.busse@stud-mail.uni-wuerzburg.de
## git: https://github.com/caro-b
## Date: July, 2021
## R version: 4.0.5
## operating system used for testing: Windows 10



#### SETUP ####

## install required packages (if not installed yet)
packagelist <- c("dplyr","gganimate","ggplot2","ggthemes","maps","osmdata","raster","rgee","rgdal","rgeos","RStoolbox","sf","tidyverse")
new.packages <- packagelist[!(packagelist %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

## load required packages
lapply(packagelist, require, character.only = TRUE)



#### DATA IMPORT ####

## get working directory
getwd()

## set directory to folder where input files lie
dir <- "C:/Users/carob/Dropbox/EAGLE/SS21/Conservation/RoadNetwork/"

## input area of interest (Bach Ma Nationalpark)
aoi <- readOGR(paste(dir, "prediction_area.shp", sep = ""))

## input manually digitized roads from IZW
roads_IZW_2014 <- readOGR(paste(dir, "2014.shp", sep = ""))

## input OSM road data
# convert bounding box to an overpass query object (API)
roads_osm <- opq(aoi@bbox) %>%
  # get roads from OSM
  add_osm_feature(key = "highway") %>%  
  # output as spatial object (sp)
  osmdata_sp()

## extract roads as line objects
roads_osm_lines <- roads_osm$osm_lines


## input Hansen GFC raster files 
# read all forest_loss files from directory into list
gfc_files <- list.files(path = dir, pattern = "forest_loss", full.names=T)

# stack raster layers
gfc_stack <- stack(gfc_files)

# convert to raster brick to shorten processing time
gfc <- brick(gfc_stack)


## input Sentinel-2 RGB as raster stack as background for plotting
#### TODO: change naming ####
sent_2020 <- stack(paste(dir, "Vt_area_2020.tif", sep = ""))



#### PRE-PROCESSING ####

## set to CRS of GFC data
# set CRS of aoi 
crs(gfc[[1]]) # CRS("+proj=longlat +datum=WGS84")
crs(aoi) <- crs(gfc[[1]])

# transform crs of roads
roads_osm_lines <- spTransform(roads_osm_lines, CRS("+proj=longlat +datum=WGS84"))
crs(roads_IZW_2014) <- crs(gfc[[1]])

## clip roads to aoi 
# byid = T to keep single features (roads)
roads_osm_clip <- gIntersection(roads_osm_lines, aoi, byid = T)
roads_IZW_2014_clip <- gIntersection(roads_IZW_2014, aoi, byid = T)

## plot discrepancies between OSM & IZW road data
plot(aoi)
plot(roads_IZW_2014_clip, add = T, col = "blue")
plot(roads_osm_clip, add = T, col = "orange")
legend("topleft", legend = c("Roads IZW (2014)", "Roads OSM (2020)"),
       col = c("blue", "orange"), lty = 1)

## combine roads from OSM & IZW
# take OSM data (as better for scaling to other areas)
# remove roads where OSM & IZW roads intersect
roads_IZW_2014_new <- roads_IZW_2014_clip[is.na(over(roads_IZW_2014_clip, roads_osm_clip))]

# only use osm roads, add additional IZW roads which aren't in OSM data
roads_union <- rbind(roads_osm_clip, roads_IZW_2014_new)

## buffer around roads (as OSM only provides lines, not the original road width)
# highway width = 22m
# highway expected to be opened in 2021
# source: https://vinlove.net/2021/04/05/expressway-through-bach-ma-national-park-before-traffic-opening-date/
roads_union_buffer <- buffer(roads_union, width = 0.0001, dissolve = F) # ~10m buffer



#### ANALYSIS ####

## function for masking forest loss areas in gfc data according to road data
maskForestLoss <- function(forest_loss, roads) {
  # masking forest loss areas in gfc loss data according to road data
  forest_loss[forest_loss < 1] <- NA
  forest_loss_roads <- mask(forest_loss, roads)
  forest_loss_roads_pol <- rasterToPolygons(forest_loss_roads, dissolve = F)
  return(forest_loss_roads_pol)
}


## Procedure: iterate over all raster layers (layers correspond to years): 
## 1. mask forest loss according to roads 
## 2. polygonize forest loss pixels
## 3. save polygons into one spatialpolygonsdataframe & add year column

# index for iteration: starts at year 2001 (no forest loss in year 2000) = second layer in raster
d <- 1

# create empty spatialpolygonsdataframe with year column
polygons_empty <- SpatialPolygons(list(), proj4string = crs(gfc[[1]]))
data_empty <- data.frame(matrix(ncol = 1, nrow = 0, dimnames = list(NULL, c("year"))))
sp_forestloss <- SpatialPolygonsDataFrame(polygons_empty, data = data_empty)

# iterate over rasterstack layers (years)
while(d < (nlayers(gfc))) {  
  # mask forest loss according to roads & polygonize raster pixels via function maskForestLoss()
  # save as spatialpolygons
  sp_forestloss_new <- SpatialPolygons(maskForestLoss(gfc[[d+1]], roads_union_buffer)@polygons, proj4string = crs(gfc[[1]]))
  # add year of forestloss
  sp_forestloss_new$year <- rep(d, length(sp_forestloss_new))
  # save all polygons into one spatialpolygonsdataframe
  sp_forestloss <- rbind(sp_forestloss, sp_forestloss_new)
  # iterate to next year
  d = d+1
}

## function to check if x% of forest loss pixels lie in road buffer
roadDev <- function(forestloss_pol, roads_buffer) {
  # calculate area of road buffers
  roads_buffer$area <- area(roads_buffer)
  # calculate area of forest loss polygons
  forestloss_pol$area <- area(forestloss_pol)
  # calculate area of forest loss polygons which lie in each road buffer
  ForestLossOnRoad <- over(roads_buffer, forestloss_pol, fn = sum)
  # get roads for which x% of forest loss polygons lie in buffer, then road dev
  roads <- roads_buffer[((ForestLossOnRoad$area / area(roads_buffer)) >= 0.05) & !is.na(ForestLossOnRoad$area),]
  return(roads)
}


## 4. stepwise (cumulative) intersection of forest loss polygons & roads to find year of road development
## a) intersect forest loss pixels of certain year with all roads
## b) if >x% pixels intersecting with road - save year in road spatialpolygonsdataframe
## c) else cumulatively add forest loss pixels to next year & repeat from a)

# index for years
y <- 1

# list for adding cumulative years
year_list <- list(y)

# add year variable to road polygons
roads_union_buffer$year <- rep(NA, length(roads_union_buffer))

# iterate over rasterstack layers (years)
while(y < nlayers(gfc)) {
  # get forest loss polygons of current year (from year list)
  polygons_current_year <- sp_forestloss[sp_forestloss$year %in% year_list,]
  # check if x% of forest loss polygons lie on road buffers, then road development
  roads <- roadDev(polygons_current_year, roads_union_buffer)
  # if no road dev for that year, add year to list 
  if (length(roads) == 0) {
    year_list <- c(year_list, y)
  } else {
  # check if already year saved in year column (to not overwrite earlier year)
  if (is.na(mean(roads_union_buffer$year[roads_union_buffer@polygons %in% roads@polygons]))) {
    # if no year saved yet (na value), save year of intersecting roads
    roads_union_buffer$year[roads_union_buffer@polygons %in% roads@polygons] <- y
  }
  # reset list to current year
  year_list <- list(y)
  }
  # iterate to next year
  y = y + 1
}

## remove roads with no corresponding forest loss between 2000 and 2020 (na in year column)
roads_final <- roads_union_buffer[!is.na(roads_union_buffer$year),]

## convert year numbers to actual years
roads_final$year <- as.numeric(roads_final$year)
roads_final$year[roads_final$year <= 9] <- paste("200", roads_final$year[roads_final$year <= 9], sep = "")
roads_final$year[!grepl("200", roads_final$year)] <- paste("20", roads_final$year[!grepl("200", roads_final$year)], sep = "")

# convert to numeric
roads_final$year <- as.numeric(roads_final$year)

# lis <- seq(1:19)
# plot(sp_forestloss[sp_forestloss$year %in% lis,])
# plot(roads,add=T)
# plot(roads_final[roads_final$year==19,], add=T, col="magenta")



#### DATA EXPORT ####

## export polygons to shapefiles
# export polygons of forest loss with corresponding years
writeOGR(sp_forestloss, "GFC_ForestLossOnRoads_20002020.shp", driver = "ESRI Shapefile", dir)

# export final roads with the corresponding years of there development
writeOGR(roads_final, "RoadDevelopmentYears_20002020.shp", driver = "ESRI Shapefile", dir)



#### VISUALIZATION ####

# first convert data to simple feature object (sf) (for easier plotting with ggplot)
roads_final_sf <- st_as_sf(roads_final)

# convert year to factor (for plotting)
roads_final_sf$year <- as.factor(roads_final_sf$year)

## basic plot
basic <- 
  ggplot(data = roads_final_sf) +
  # add RGB plot of Sentinel2 data as background
  # set scale to maximum pixel value of the 3 raster bands
  ggRGB(sent_2020, 1, 2, 3, scale = 0.3, ggLayer = T) +
  # add road data with year of development as coloring factor
  geom_sf(aes(color = year, fill = year)) +
  # remove background grid & axes
  theme_map()


## map of road development
basic +
  # change legend title
  labs(
    #fill = "Year",
    title = "Road Development", 
    subtitle = "Bach Ma Nationalpark (Vietnam), 2000-2020") +
  # remove legend for color parameter
  guides(#color = F,
         # write fill legend into 2 rows
         colour = guide_legend(nrow = 2, byrow = TRUE,
         keyheight = unit(2, units = "mm"))
  ) +
  theme(
    # change title format
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 12),
    # change legend position
    legend.position = "bottom", legend.justification = 'right',
    # change legend title formatting
    legend.title = element_text(face = "bold"),
    legend.key = element_rect(colour = "white"))


## animation of consecutive road development per year
basic +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 14, face = "bold")) +
  # Here comes the gganimate part
  labs(title = "Road Development \nBach Ma Nationalpark (Vietnam), 2000-2020", subtitle = "Year: {current_frame}") +
  transition_manual(year, cumulative = T) +
  ease_aes('linear')

