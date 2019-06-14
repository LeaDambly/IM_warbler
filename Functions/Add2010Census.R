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
  colnames(PopData) <- gsub("density", "pop.density", colnames(PopData))
  data@data <- cbind(data@data, PopData[,c("FIPS", "pop.density")])
  data
}

