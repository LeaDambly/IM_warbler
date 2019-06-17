library(plyr)
library(dplyr)
library(rBBS)
library(USAboundaries)
library(sf)
library(sp)
library(raster)
library(rgdal)
library(RColorBrewer)
library(censusapi)
library(spocc)

library(readr)
library(elevatr)
library(FedData)
proj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Local Location of file with data from Miller et al. (2019)
# Data downloaded from here:
# https://besjournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2F2041-210X.13110&file=mee313110-sup-0001-supplementA.zip
Miller.file <- "Data/mee313110-sup-0001-supplementa.zip"

# Read relevant functions
source("Functions/GetBBSData.R")
source("Functions/Add2010Census.R")

# Get outline of PA
PA <- us_states(states = "Pennsylvania")
PA <- PA$geometry[1]
PA <- as(PA, 'Spatial')


# BBA from Miller et al. appendix

if(!file.exists("Data/BBA.csv")) {
  load(unzip(Miller.file, files="DI_Data.Rdata"))
  
  # sum counts (counts are undertaken in 5 time intervals at a single spatial point)
  BBA_Wren <- bba %>%
    mutate(total = dplyr::select(., v1:v5) %>% rowSums(na.rm = TRUE)) %>%
    dplyr::select(-c(v1:v5)) %>%
    mutate(present = if_else(total == 0, FALSE, TRUE)) %>%
    rename(X = Longitude, Y = Latitude)
  write.csv(BBA_Wren, file="Data/BBA.csv")
} else {
  BBA_Wren <- read.csv(file="Data/BBA.csv")
}

# change to spatial dataframe
BBA_sp <- SpatialPointsDataFrame(BBA_Wren[,c("X", "Y")],
                                 data = BBA_Wren[,c("present", "point")],
                                 proj4string = crs(proj))



# Get BBS data (using rBBS package)

if(!file.exists("Data/BBS.csv")) {
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
  write.csv(BBS_Wren, file="Data/BBS.csv")
} else {
  BBS <- read.csv(file="Data/BBS.csv")
}

BBS_sp <- SpatialPointsDataFrame(BBS_Wren[,c("X", "Y")],
                                 data = BBS_Wren[,c("NPres", "Ntrials")],
                                 proj4string = crs(proj))



# eBird, downloaded from GBIF
 if(!file.exists("Data/eBird.csv")) {
  eBird.raw <- occ("Setophaga caerulescens", from="gbif", date=c("2005-01-01", "2005-12-31"), geometry=PA@bbox)$gbif
  eBird <- eBird.raw$data$Setophaga_caerulescens[grep("EBIRD", eBird.raw$data$Setophaga_caerulescens$collectionCode), 
                                                 c("longitude", "latitude", "year")]
  
  # make into spatial points
  eBird_coords <- cbind(eBird$longitude, eBird$latitude)
  colnames(eBird_coords) <- c("X", "Y")
  write.csv(eBird_coords, file="Data/eBird.csv", row.names = FALSE)
} else {
  eBird_coords <- read.csv(file="Data/eBird.csv")
}

eBird_pts <- SpatialPoints(coords = eBird_coords, proj4string = proj)
# trim to keep only those occuring in PA (with probably unnecessary back-and-forth of data formats)
eBird_pts <- over(eBird_pts, PA)

eBird_pts <- data.frame(eBird_coords[!is.na(eBird_pts),])
eBird_sp <- SpatialPoints(coords = eBird_pts, proj4string = proj)


# Covariates

elev <- get_elev_raster(PA, z = 6, clip = "locations") #z = 1 for lowest res, z = 14 for highest (DL time very long)
elevation <- as.data.frame(elev, xy = T, na.rm = T)

elev <- get_elev_raster(PA, z = 6, clip = "locations") #z = 1 for lowest res, z = 14 for highest (DL time very long)
elevation <- as.data.frame(elev, xy = T, na.rm = T)
elevation$layer[elevation$layer<0] <- 0
# elevation data using elevatr (could theoretically also use FedData but get holes in elev raster)

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


# Add population density

censuskey <- "ba4f7dd49b22e58b5e6a5cc99b349b555bbf95a8"
covariates_eBird <- Add2010Census(covariates, proj, censuskey)
covariates_eBird@data$FIPS <- as.numeric(covariates_eBird@data$FIPS)
covariates_eBird@data <-data.frame(apply(covariates_eBird@data, 2 , scale))


# Save the data
save(proj, PA, BBA_sp, BBS_sp, eBird_sp, covariates, covariates_eBird, file="Data/BTWarblerData.RData")
