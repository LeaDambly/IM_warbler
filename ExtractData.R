library(plyr)
library(dplyr)
library(rBBS)
library(USAboundaries)
library(sf)
library(sp)
library(raster)
library(rgdal)
library(RColorBrewer)
library(PointedSDMs)
# library(mapview)
library(censusapi)

library(readr)
library(elevatr)
library(FedData)
proj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# key for US census data
censuskey <- "499974bdafb06dc0d903c9b082e8d37e2f9e08ce"

source("warblerfunctions.R")
# Get outline of PA
PA <- us_states(states = "Pennsylvania")
PA <- PA$geometry[1]
PA <- as(PA, 'Spatial')

# BBA from Miller et al. appendix
load("Data/PA_BBA.RData")

# sum counts (counts are undertaken in 5 time intervals at a single spatial point)
BBA_Wren <- bba %>%
  mutate(total = dplyr::select(., v1:v5) %>% rowSums(na.rm = TRUE)) %>%
  dplyr::select(-c(v1:v5)) %>%
  mutate(present = if_else(total == 0, FALSE, TRUE)) %>%
  rename(X = Longitude, Y = Latitude)

# change to spatial dataframe
BBA_sp <- SpatialPointsDataFrame(BBA_Wren[,c("X", "Y")],
                                 data = BBA_Wren[,c("present", "point")],
                                 proj4string = crs(proj))

# BBS (using rBBS package)
RegionMetaData <- GetRegions()
WeatherMetaData <- GetWeather()
RoutesMetaData <- GetRoutes()

PACode <- RegionMetaData$RegionCode[RegionMetaData$`State/Prov/TerrName` == "PENNSYLVANIA"]
PAYears <- 2005:2009

# fixed getRouteData
RegionsForZipFiles <- GetRegions(ZipFiles = TRUE)


BBS_Wren <- GetRouteData(AOU=6540, countrynum = 840, states = PACode, year = PAYears,
                         weather = WeatherMetaData, routes = RoutesMetaData, TenStops = FALSE,
                         Zeroes = TRUE)

# counts are made along a route
# need to be made into number of presences and number of trials (routes)
BBS_Wren <- BBS_Wren %>%
  mutate(NPres = rowSums(dplyr::select(., starts_with("stop")) > 0)) %>%
  mutate(Ntrials = rowSums(!is.na(dplyr::select(., starts_with("stop"))))) %>%
  group_by(route) %>%
  summarise(Ntrials = sum(Ntrials),
            NPres = sum(NPres),
            Latitude = first(Latitude),
            Longitude = first(Longitude)) %>%
  rename(X = Longitude, Y = Latitude)

# change to spatial points
BBS_Wren <- as.data.frame(BBS_Wren)

BBS_sp <- SpatialPointsDataFrame(BBS_Wren[,c("X", "Y")],
                                 data = BBS_Wren[,c("NPres", "Ntrials")],
                                 proj4string = crs(proj))



# eBird / GBIF
# get GBIF data, filter for 2005-2009: eBird (actually Lab of O) only 
GBIF <- readr::read_delim(unzip("~/Dropbox/DataIntegration/0023215-190415153152247.zip"), "\t", escape_double = FALSE, 
                   trim_ws = TRUE)

# make into spatial points
GBIF <- GBIF %>% 
  dplyr::select(decimalLatitude, decimalLongitude, year) %>%
  filter(year %in% c(2005:2009)) %>%
  rename(X = decimalLongitude, Y = decimalLatitude)

# make into spatial points
GBIF_coords <- cbind(GBIF$X, GBIF$Y)
GBIF_pts <- SpatialPoints(coords = GBIF_coords, proj4string = proj)
# trim to keep only those occuring in PA (with probably unnecessary back-and-forth of data formats)
GBIF_pts <- over(GBIF_pts, PA)
GBIF_ptsC <- Add2010Census(data=GBIF_pts, proj=proj, censuskey=censuskey)
  
GBIF_pts <- data.frame(GBIF_coords[!is.na(GBIF_pts),])
GBIF_sp <- SpatialPoints(coords = GBIF_pts, proj4string = proj)

# Covariates

# elevation data using elevatr (could theoretically also use FedData but get holes in elev raster)
elev <- get_elev_raster(PA, z = 8, clip = "locations") #z = 1 for lowest res, z = 14 for highest (DL time very long)
elevation <- as.data.frame(elev, xy = T, na.rm = T)

# canopy from the NLCD
NLCD_canopy <- get_nlcd(template = PA, year = 2011, dataset = "canopy", label = "PA_lc")
NLCD_canopy <- projectRaster(from = NLCD_canopy, to = elev)
NLCD_canopy <- mask(NLCD_canopy, elev)
canopy <- as.data.frame(NLCD_canopy, xy = T, na.rm = T)

covariates <- full_join(elevation, canopy, by = c("x", "y"))
covariates <- covariates %>%
  rename(elevation = layer, canopy = PA_lc_NLCD_2011_canopy, X = x, Y = y)
covariates <- SpatialPointsDataFrame(covariates[,c("X", "Y")], 
                                     data = covariates[,c("elevation", "canopy")], proj4string = crs(proj))
covariates@data <-data.frame(apply(covariates@data, 2 , scale))  # scale the covariates


censuskey <- "ba4f7dd49b22e58b5e6a5cc99b349b555bbf95a8"
covariates_gbif <- Add2010Census(covariates, proj, censuskey)
covariates_gbif@data$FIPS <- as.numeric(covariates_gbif@data$FIPS)
covariates_gbif@data <-data.frame(apply(covariates_gbif@data, 2 , scale))



save(proj, PA, BBA_sp, BBS_sp, GBIF_sp, covariates, covariates_gbif, file="Data/BTWarblerData.RData")
