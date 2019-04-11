# Integrated Model of black-throated warbler in Pennsylvania
# 11.04.19 - Lea Dambly

# Data: BBS, eBird, PA bird atlas
# Cov: Elevation and forest cover

# libraries----
library(rBBS)

# data prep----
# PA BBA from Miller et al. appendix
# BBS
RegionMetaData <- GetRegions()
WeatherMetaData <- GetWeather()
# Bob FYI, your GetRoutes function has an issue. It calls Routes.zip instead of routes.zip 
GetRoutes <- function(Dir="ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/") {
  routes=GetUnzip(ZipName=paste0(Dir, "routes.zip"), FileName="routes.csv")
  routes$routeID=paste(routes$statenum, routes$Route)
  routes
}
RoutesMetaData <- GetRoutes()

PACode <- RegionMetaData$RegionCode[RegionMetaData$`State/Prov/TerrName` == "PENNSYLVANIA"]
PAYears <- 2005:2009

PAWren<- GetRouteData(AOU=6540, countrynum = 840, states = PACode, year = PAYears,
                      weather = WeatherMetaData, routes = RoutesMetaData,
                      Zeroes = TRUE)
# leads to warnings();  In FUN(X[[i]], ...) : no ID, so setting to NA
# some issue with IDs in getID function not matching I assume - will work on it
