# required packages
library(raster)
library(sf)
library(USAboundaries)
library(elevatr)
library(FedData)
library(dplyr)
library(censusapi)

# required function
source("Functions/Add2010Census.R")

# the following packages need to be installed to download and format from
# scratch the data that is already included in this repository:

# devtools::install_github("oharar/rBBS")
# install.packages("plyr")
# install.packages("spocc")

# coordinate reference system to use throughout
proj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Get outline of PA
PA <- USAboundaries::us_states(states = "Pennsylvania")
PA <- PA$geometry[1]
PA <- as(PA, "Spatial")

# BBA from Miller et al. appendix
if (!file.exists("Data/BBA.csv")) {

  # Local Location of file with data from Miller et al. (2019)
  # Data downloaded from here:
  Miller.url <- "https://besjournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2F2041-210X.13110&file=mee313110-sup-0001-supplementA.zip"
  Miller.file <- "Data/mee313110-sup-0001-supplementa.zip"

  # download the Miller file if needed
  if (!file.exists(Miller.file)) {
    download.file(Muller.url, Miller.file)
  }

  load(unzip(Miller.file, files = "DI_Data.Rdata"))

  # sum counts (counts undertaken in 5 time intervals at a single spatial point)
  BBA_Wren <- bba %>%
    mutate(total = dplyr::select(., v1:v5) %>% rowSums(na.rm = TRUE)) %>%
    dplyr::select(-c(v1:v5)) %>%
    mutate(present = if_else(total == 0, FALSE, TRUE)) %>%
    dplyr::rename(X = Longitude, Y = Latitude)

  write.csv(BBA_Wren, file = "Data/BBA.csv")
} else {
  BBA_Wren <- read.csv(file = "Data/BBA.csv")
}

# change to spatial dataframe
BBA_sp <- SpatialPointsDataFrame(
  coords = BBA_Wren[, c("X", "Y")],
  data = BBA_Wren[, c("present", "point")],
  proj4string = crs(proj)
)



# Get BBS data (using rBBS package)
if (!file.exists("Data/BBS.csv")) {
  library(rBBS)
  source("Functions/GetBBSData.R")
  ldply <- plyr::ldply

  RegionMetaData <- GetRegions()
  WeatherMetaData <- GetWeather()
  RoutesMetaData <- GetRoutes()

  idx <- RegionMetaData[["State/Prov/TerrName"]] == "PENNSYLVANIA"
  PACode <- RegionMetaData$RegionCode[idx]
  PAYears <- 2005:2009

  # fixed getRouteData
  RegionsForZipFiles <- GetRegions(ZipFiles = TRUE)

  BBS_Wren <- GetRouteData(
    AOU = 6540,
    countrynum = 840,
    states = PACode,
    year = PAYears,
    weather = WeatherMetaData,
    routes = RoutesMetaData,
    TenStops = FALSE,
    Zeroes = TRUE
  )

  # counts are made along a route
  # need to be made into number of presences and number of trials (routes)
  BBS_Wren <- BBS_Wren %>%
    mutate(NPres = rowSums(dplyr::select(., starts_with("stop")) > 0)) %>%
    mutate(Ntrials = rowSums(!is.na(dplyr::select(., starts_with("stop"))))) %>%
    group_by(route) %>%
    summarise(
      Ntrials = sum(Ntrials),
      NPres = sum(NPres),
      Latitude = first(Latitude),
      Longitude = first(Longitude)
    ) %>%
    dplyr::rename(X = Longitude, Y = Latitude)

  # change to spatial points
  BBS_Wren <- as.data.frame(BBS_Wren)
  write.csv(BBS_Wren, file = "Data/BBS.csv")
} else {
  BBS_Wren <- read.csv(file = "Data/BBS.csv")
}

BBS_sp <- SpatialPointsDataFrame(
  coords = BBS_Wren[, c("X", "Y")],
  data = BBS_Wren[, c("NPres", "Ntrials")],
  proj4string = crs(proj)
)

# eBird, downloaded from GBIF
if (!file.exists("Data/eBird.csv")) {
  eBird.raw <- spocc::occ(
    query = "Setophaga caerulescens",
    from = "gbif",
    date = c("2005-01-01", "2005-12-31"),
    geometry = PA@bbox
  )$gbif

  rows <- grep("EBIRD", eBird.raw$data$Setophaga_caerulescens$collectionCode)
  cols <- c("longitude", "latitude", "year")
  eBird <- eBird.raw$data$Setophaga_caerulescens[rows, cols]

  # make into spatial points
  eBird_coords <- cbind(eBird$longitude, eBird$latitude)
  colnames(eBird_coords) <- c("X", "Y")
  write.csv(eBird_coords, file = "Data/eBird.csv", row.names = FALSE)
} else {
  eBird_coords <- read.csv(file = "Data/eBird.csv")
}

eBird_pts <- SpatialPoints(coords = eBird_coords, proj4string = proj)
# trim to keep only those occuring in PA (with probably unnecessary
# back-and-forth of data formats)
eBird_pts <- over(eBird_pts, PA)

eBird_pts <- data.frame(eBird_coords[!is.na(eBird_pts), ])
eBird_sp <- SpatialPoints(coords = eBird_pts, proj4string = proj)


# Covariates

# elevation data using elevatr (could theoretically also use FedData but get
# holes in elev raster)
elev <- elevatr::get_elev_raster(PA, z = 6, clip = "locations")
# z = 1 for lowest res, z = 14 for highest (DL time very long)
elevation <- as.data.frame(elev, xy = TRUE, na.rm = TRUE)
elevation$layer[elevation$layer < 0] <- 0

# canopy from the NLCD
NLCD_canopy <- get_nlcd(
  template = PA,
  year = 2011,
  dataset = "canopy",
  label = "PA_lc"
)
NLCD_canopy <- projectRaster(from = NLCD_canopy, to = elev)
NLCD_canopy <- mask(NLCD_canopy, elev)
canopy <- as.data.frame(NLCD_canopy, xy = TRUE, na.rm = TRUE)

covariates <- full_join(elevation, canopy, by = c("x", "y"))
covariates <- covariates %>%
  dplyr::rename(
    elevation = layer,
    canopy = PA_lc_NLCD_2011_canopy,
    X = x, Y = y
  )

covariate_coords <- covariates[, c("X", "Y")]
covariate_data <- covariates[, c("elevation", "canopy")]
covariates <- SpatialPointsDataFrame(
  coords = covariate_coords,
  data = covariate_data,
  proj4string = crs(proj)
)
# scale the covariates
covariates@data <- data.frame(apply(covariates@data, 2, scale))


# Add population density

censuskey <- "ba4f7dd49b22e58b5e6a5cc99b349b555bbf95a8"
covariates_eBird <- Add2010Census(covariates, proj, censuskey)
covariates_eBird@data$FIPS <- as.numeric(covariates_eBird@data$FIPS)
covariates_eBird@data <- data.frame(apply(covariates_eBird@data, 2, scale))


# Save the data
save(
  proj, PA, BBA_sp, BBS_sp, eBird_sp, covariates, covariates_eBird,
  file = "Data/BTWarblerData.RData"
)
