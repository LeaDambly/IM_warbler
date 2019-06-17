GetRouteData <- function(AOU=NULL, countrynum=NULL, states=NULL, year, weather=NULL, routes=NULL, 
                         Zeroes=TRUE, TenStops = FALSE, 
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

# GetRoutes function has an issue. It calls Routes.zip instead of routes.zip
GetRoutes <- function(Dir="ftp://ftpext.usgs.gov/pub/er/md/laurel/BBS/DataFiles/") {
  routes=GetUnzip(ZipName=paste0(Dir, "routes.zip"), FileName="routes.csv")
  routes$routeID=paste(routes$statenum, routes$Route)
  routes
}


# population data from US 2010 census
# Function to add population data (from 2010 census) to a data frame # Arguments
#   data: spatial points data frame to add the data to
#   proj: Projection (CRS()))
#   censuskey: Census key to get population data.
Add2010Census <- function(data, proj, censuskey) {
  require(censusapi)
  require(maps)
  require(maptools)
  require(rgeos)
  # Decennial Census sf3, 2010
  #  censuskey <- ....
  # get population data
  data2010 <- getCensus(name="dec/sf1", vintage=2010, key=censuskey,
                        vars="P008001", region="county:*")
  rownames(data2010) <- paste0(data2010$state, data2010$county)
  
  # Get county maps, and make into spatial polygon
  county <- map('county', fill=TRUE)
  data(county.fips)
  county.fips$fips <- as.character(county.fips$fips)
  county.fips$fips[nchar(county.fips$fips)==4] <- paste0("0",
                                                         county.fips$fips[nchar(county.fips$fips)==4])
  
  county$FIPS <- county.fips$fips[match(map("county",
                                            plot=FALSE)$names, county.fips$polyname)]
  CountyMap <- map2SpatialPolygons(map=county, IDs = county$FIPS, proj4string = proj)
  Countydat <- data.frame(area = unlist(lapply(CountyMap@polygons,
                                               function(ply) ply@area)),
                          FIPS = unlist(lapply(CountyMap@polygons,
                                               function(ply) ply@ID)),
                          row.names = unlist(lapply(CountyMap@polygons,
                                                    function(ply) ply@ID)),
                          stringsAsFactors = FALSE)
  
  # Merge population and county data
  Countydat.m <- merge(Countydat, data2010, by="row.names")
  
  rownames(Countydat.m) <- rownames(Countydat)
  Countydat.m$density <- Countydat.m$P008001/Countydat.m$area
  CountyPops <- SpatialPolygonsDataFrame(data=Countydat.m, Sr=CountyMap)
  
  # extract data from correct polygon for data, and merge
  PopData <- over(x=data, y=CountyPops)
  data@data <- cbind(data@data, PopData[,c("FIPS", "density")])
  data
}


# just took PointedSDM's function apart a bit because buffer function used in MakeSpatialRegion doesn't work with unprojected data - so no buffering here
MakeSpatialRegion2 <- function (data = NULL, coords = c("X", "Y"), meshpars, bdry = NULL, 
                                proj = CRS("+proj=utm")) {
  require(rgeos)
  region.bdry <- inla.sp2segment(bdry)
mesh <- inla.mesh.2d(boundary = region.bdry, cutoff = meshpars$cutoff, 
                     max.edge = meshpars$max.edge, offset = meshpars$offset)
spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
dd <- deldir::deldir(mesh$loc[, 1], mesh$loc[, 2])
tiles <- deldir::tile.list(dd)
poly.gpc <- as(bdry@polygons[[1]]@Polygons[[1]]@coords, "gpc.poly")
w <- sapply(tiles, function(p) rgeos::area.poly(rgeos::intersect(as(cbind(p$x, 
                                                                          p$y), "gpc.poly"), poly.gpc)))
return(list(mesh = mesh, spde = spde, w = w))
}



# modified function to take PC priors
FitModel2 <- function (..., formula = NULL, CovNames = NULL, mesh, spat.ind = "i", 
                       predictions = FALSE, tag.pred = "pred", control.fixed = NULL, 
                       waic = FALSE, dic = FALSE, nthreads = NULL) 
{
  stck <- inla.stack(...)
  if (is.null(CovNames)) {
    CovNames <- unlist(stck$effects$names)
    CovNames <- CovNames[!CovNames %in% c(spat.ind)]
  }
  else {
    if (!is.null(formula)) {
      warning("CovNames and formula are both not NULL: CovNames will be ignored")
    }
  }
  mesh <- inla.spde2.matern(mesh)
  if (!is.null(spat.ind)) {
    CovNames <- c(CovNames, paste0("f(", spat.ind, ", model=mesh)"))
  }
  if (is.null(control.fixed)) {
    control.fixed <- list(mean = 0)
  }
  if (is.null(formula)) {
    Formula <- formula(paste(c("resp ~ 0 ", CovNames), collapse = "+"))
  }
  else {
    if (is.null(spat.ind)) {
      Formula <- formula
    }
    else {
      if (any(grepl(paste0("(", spat.ind, ","), formula, 
                    fixed = TRUE))) {
        warning(paste0(spat.ind, " already in formula, so will be ignored"))
        Formula <- formula
      } else {
        Formula <- update(formula, paste0(" ~ . + f(", spat.ind, 
                                          ", model=mesh)"))
      }
    }
  }
  mod <- inla(Formula, family = c("poisson", "binomial"), control.family = list(list(link = "log"), 
                                                                                list(link = "cloglog")), data = inla.stack.data(stck), 
              verbose = FALSE, control.results = list(return.marginals.random = FALSE, 
                                                      return.marginals.predictor = FALSE), control.predictor = list(A = inla.stack.A(stck), 
                                                                                                                    link = NULL, compute = TRUE), control.fixed = control.fixed, 
              Ntrials = inla.stack.data(stck)$Ntrials, E = inla.stack.data(stck)$e, 
              num.threads = nthreads,
              control.compute = list(waic = waic, dic = dic))
  if (predictions) {
    id <- inla.stack.index(stck, tag.pred)$data
    pred <- data.frame(mean = mod$summary.fitted.values$mean[id], 
                       stddev = mod$summary.fitted.values$sd[id])
    res <- list(model = mod, predictions = pred)
  }
  else {
    res <- mod
  }
  res
}
