library(USAboundaries)
library(sf)
library(sp)
library(raster)
library(rgdal)
library(elevatr)
library(FedData)
library(RColorBrewer)
library(PointedSDMs)
library(mapview)


load("Data/Stacks.RData")
load("Data/WarblerModelOutput.RData")
load("Data/BTWarblerData.RData")

proj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

Pred <- SpatialPixelsDataFrame(points = stk.pred$predcoords, data = warbler_model$predictions, proj4string = crs(proj))
Pred@data$precision <- Pred@data$stddev^-2


# Pennsylvania border
PA <- us_states(states = "Pennsylvania")
PA <- PA$geometry[1]
PA <- as(PA, 'Spatial')

# Plot of data

png("DataPlot.png", height=480, width = 640)
par(mar=rep(0.1,4))
plot(PA)
points(BBA_sp[which(BBA_sp@data$present == FALSE),], cex=0.2, col="pink") # Add BBA
points(BBA_sp[which(BBA_sp@data$present == TRUE),], cex=0.2, col="red3")
points(BBS_sp, cex=0.5, col="blue", pch=16) # Add BBS
points(GBIF_sp, cex=0.5, col="sandybrown", pch=16) # Add eBird
legend(-79, 42.45, c("BBA, absent", "BBA, present", "BBS route", "eBird"), col=c("pink", "red3", "blue", "sandybrown"), pch=c(1,1,16,16), 
       ncol = 2)
dev.off()


# Plot of predictions

png("PredPlot.png", height=360, width = 640)
par(mar=rep(0.1,4))
plot(Pred, col=grey(seq(0,1,length=100)))
lines(PA)
dev.off()


# Plot of uncertainty
png("StdDevPlot.png", height=360, width = 640)
par(mar=rep(0.1,4))
plot(Pred, attr=3, col=0)
plot(Pred, attr = 'stddev', col=grey(seq(0,1,length=100)))
lines(PA)
dev.off()


GreyCol <- function(x) {
  x.p <- (x - min(x))/(max(x) - min(x))
  grey(x.p)
}
plot(Pred@coords, col=GreyCol(Pred@data$precision))


