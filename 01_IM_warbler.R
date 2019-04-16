# Integrated Model of black-throated warbler in Pennsylvania
# 11.04.19 - Lea Dambly

# Data: BBS, eBird, PA bird atlas
# Cov: Elevation and forest cover

# libraries----
library(plyr)
library(dplyr)
library(rBBS)
library(USAboundaries)
library(sf)
library(sp)
library(readr)
library(raster)
library(rgdal)
library(FedData)
library(RColorBrewer)

# data prep----
# get PA outline
PA = us_states(states = "Pennsylvania")
PA <- PA$geometry[1]
PA <- as(PA, 'Spatial')
proj <- proj4string(PA)

# PA BBA from Miller et al. appendix (can I avoid to hardcode here?)
load("N:/IM_warbler/warbler_data.RData")

# merge count intervals
BBA_Wren <- bba %>%
  mutate(total = dplyr::select(., v1:v5) %>% rowSums(na.rm = TRUE)) %>%
  dplyr::select(-c(v1:v5)) %>%
  mutate(present = if_else(total == 0, FALSE, TRUE))

BBA_sp <- SpatialPointsDataFrame(BBA_Wren[,c("Longitude", "Latitude")],
                               data = BBA_Wren[,c("present", "point")],
                               proj4string = crs(proj))


# BBS
RegionMetaData <- GetRegions()
WeatherMetaData <- GetWeather()
# GetRoutes function has an issue. It calls Routes.zip instead of routes.zip
GetRoutes <- function(Dir="ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/") {
  routes=GetUnzip(ZipName=paste0(Dir, "routes.zip"), FileName="routes.csv")
  routes$routeID=paste(routes$statenum, routes$Route)
  routes
}
RoutesMetaData <- GetRoutes()

PACode <- RegionMetaData$RegionCode[RegionMetaData$`State/Prov/TerrName` == "PENNSYLVANIA"]
PAYears <- 2005:2009

# fixed getRouteData
RegionsForZipFiles <- GetRegions(ZipFiles = TRUE)
GetRouteData <- function(AOU=NULL, countrynum=NULL, states=NULL, year, weather=NULL, routes=NULL, 
                         Zeroes=TRUE, TenStops = TRUE, 
                         Dir="ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/") {
  
  if(TenStops) {
    DirData <- paste0(Dir, "States/")
    CountString <- "^count"
  } else {
    if(any(year<1997)) stop("Data only available from 1997: pre-1997 data not integrated into this function for 50 stop data (yet)")
    DirData <- paste0(Dir, "50-StopData/1997ToPresent_SurveyWide/")
    CountString <- "^stop"
  }
  if(!is.null(countrynum) & any(!(countrynum%in%c(124, 484, 840)))) stop("countrynum should be either 124 (Canada), 484 (Mexico), or 840 (USA)")
  
  GetDat <- function(file, dir, year, AOU, countrynum, states) {
    dat <- GetUnzip(ZipName=paste0(dir, file), FileName=gsub("^Fifty", "fifty", gsub("zip", "csv", file)))
    names(dat) <- tolower(names(dat))
    if(is.null(year)) {  UseYear <- TRUE  } else {  UseYear <- dat$year%in%year  }
    if(is.null(AOU)) {  UseAOU <- TRUE  } else {  UseAOU <- dat$aou%in%AOU  }
    if(is.null(countrynum)) {  UseCountry <- TRUE  } else {  UseCountry <- dat$countrynum%in%countrynum  }
    if(is.null(states)) {  UseState <- TRUE  } else {  UseState <- dat$statenum%in%states  }
    Use <- UseYear & UseAOU & UseCountry & UseState
    if(sum(Use)>0) {
      dat$routeID <- paste(dat$statenum, dat[,grep("^[Rr]oute$", names(dat))])
      dat <- subset(dat, subset=Use)
      return(dat)      
    } else return(NULL)
  }
  
  # Only use the files we want
  CountriesToUse <- if(!is.null(countrynum)) {
    RegionsForZipFiles$countrynum%in%countrynum 
  } else {
    TRUE
  }
  StatesToUse <- if(!is.null(states)) {
    RegionsForZipFiles$RegionCode%in%states 
  } else {
    TRUE
  }
  ToUse <- CountriesToUse & StatesToUse
  if(TenStops) {
    Files <- RegionsForZipFiles$FileName10stop[ToUse]
    Missing <- ToUse & is.na(RegionsForZipFiles$FileName10stop)
  } else { # 50 stop
    Files <- RegionsForZipFiles$FileName50stop[ToUse]
    Missing <- ToUse & is.na(RegionsForZipFiles$FileName50stop)
  }
  
  if(length(Files)==0) stop("No data for the states specified")
  if(any(is.na(Files))) warning(paste0("No data for these states: ", paste(RegionsForZipFiles$'State/Prov/TerrName'[Missing], collapse=", ")))
  
  Data.lst <- sapply(Files[!is.na(Files)], GetDat, dir=DirData, year=year, AOU=AOU, countrynum=countrynum, states=states, simplify=FALSE)
  
  if(all(unlist(lapply(Data.lst, is.null)))) {
    warning("no data, sorry")
    AllData <- NULL
  } else {
    Data <- ldply(Data.lst)
    # Get route data for all routes, and annual data
    if(is.null(weather)) weather <-GetWeather(Dir)
    if(is.null(year)) {  UseYear <- TRUE  } else {  UseYear <- weather$Year%in%year  }
    if(is.null(countrynum)) {  UseCountry <- TRUE  } else {  UseCountry <- weather$CountryNum%in%countrynum  }
    if(is.null(states)) {  UseState <- TRUE  } else {  UseState <- weather$StateNum%in%states  }
    UseWeather <- UseYear & UseCountry & UseState
    
    if(is.null(routes)) routes <- GetRoutes(Dir)
    if(is.null(countrynum)) {  UseCountry <- TRUE  } else {  UseCountry <- routes$CountryNum%in%countrynum  }
    if(is.null(states)) {  UseState <- TRUE  } else {  UseState <- routes$StateNum%in%states  }
    UseRoutes <- UseCountry & UseState
    
    CommonNames <- names(Data)[names(Data)%in%names(weather)]
    CommonNames <- CommonNames[CommonNames%in%names(routes)]
    
    # Subset data
    # First, sites sampled in chosen year(s)
    weather <-subset(weather, subset=UseWeather, 
                     select=c(CommonNames, "Year", "Month", "Day", "RunType", "StateNum"))
    # Route data for sites sampled in chosen years
    routes <- subset(routes, subset=UseRoutes & routes$routeID%in%weather$routeID, 
                     select=c(CommonNames, "Latitude", "Longitude", "StateNum"))
    
    # merge data sets
    dat.routeID.year <- paste(Data$routeID, Data$year, sep=".")
    routes$routeID <- paste0(routes$StateNum, routes$routeID)
    weather.routeID.year <- paste(paste0(weather$StateNum, weather$routeID), weather$Year, sep=".")
    
    WeatherWhiches <- match(dat.routeID.year, weather.routeID.year)
    
    RouteWhiches <- match(Data$routeID, routes$routeID)
    
    AllData <- cbind(Data, weather[WeatherWhiches, !names(weather)%in%names(Data)],
                     routes[RouteWhiches, !names(routes)%in%names(Data)])
    
    #  if(!is.na(weather)) AllData <- merge(Data, weather, all=TRUE) # by=c("routeID", "RPID"), 
    #  if(!is.na(routes))  AllData <- merge(AllData, routes, all=TRUE) # by="routeID", 
    AllData$SumCount <- apply(AllData[,grep(CountString, names(AllData))],1,sum, na.rm=TRUE)
    if(!Zeroes) AllData <- subset(AllData, AllData$SumCount>0)
    AllData <- AllData[,!names(AllData)%in%c(".id", "routedataid", "year")]
  }
  
  AllData
}


BBS_Wren <- GetRouteData(AOU=6540, countrynum = 840, states = PACode, year = PAYears,
                      weather = WeatherMetaData, routes = RoutesMetaData,
                      Zeroes = TRUE)
BBS_sp <- SpatialPoints(BBS_Wren[,c("Longitude", "Latitude")], proj4string = crs(proj))

# eBird / GBIF
# get GBIF data, filter for 2005-2009
GBIF <- read_delim("GBIF.csv", "\t", escape_double = FALSE, 
                   trim_ws = TRUE)
GBIF <- GBIF %>% 
  dplyr::select(decimalLatitude, decimalLongitude, year) %>%
  filter(year %in% c(2005:2009))

# make into spatial points
GBIF_coords <- cbind(GBIF$decimalLongitude, GBIF$decimalLatitude)
GBIF_pts <- SpatialPoints(coords = GBIF_coords, proj4string = CRS(proj))

# trim to keep only those occuring in PA
GBIF_pts <- over(GBIF_pts, PA)
GBIF_pts <- data.frame(GBIF_coords[!is.na(GBIF_pts),])
GBIF_sp <- SpatialPoints(coords = GBIF_pts, proj4string = CRS(proj))

# plot all together to double check
cols <- c(brewer.pal(4, "Paired"))
names(cols) <- c("BBS", "GBIF", "BBA, absent", "BBA, present")
plot(PA)
points(BBS_sp, cex = 0.5, pch = 19, col = cols["BBS"])
points(GBIF_sp, cex = 0.5, pch = 19, col = cols["GBIF"])
points(BBA_sp, cex = 0.7, pch = 19, col = cols[c("BBA, absent", "BBA, present")][1+BBA_sp@data$present])
legend("top", legend = names(cols), fill = cols, cex = 0.8)

# covariates----
# elevation data using elevatr (could theoretically also use FedData but get holes in elev raster)
elev <- get_elev_raster(PA, z = 8, clip = "locations") #z = 1 for lowest res, z = 14 for highest (DL time very long)
plot(elev, alpha = 0.2, add = T) # just to double-check

# canopy from the NLCD
NLCD_canopy <- get_nlcd(template = PA, year = 2011, dataset = "canopy", label = "PA_lc")
NLCD_canopy <- projectRaster(from = NLCD_canopy, to = elev)
NLCD_canopy <- mask(NLCD_canopy, elev)
plot(NLCD_canopy, alpha = 0.2, add = T) # again, just to check



