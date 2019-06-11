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

png("PredPlot.png", height=240)
plot(Pred, col=grey(seq(0,1,length=100)))
dev.off()



Pred <- SpatialPixelsDataFrame(points = stk.pred$predcoords, data = warbler_model$predictions, proj4string = crs(proj))
Pred@data$precision <- Pred@data$stddev^-2
ncolours <- 200
meancols.fn <- colorRampPalette(brewer.pal(9, 'Oranges'))
meancols <- meancols.fn(ncolours)
map.mean <- mapview(Pred, zcol = c("mean"), legend = TRUE, 
                    col.regions=meancols, layer.name = "Mean pred")

sdcols.fn <- colorRampPalette(brewer.pal(9, 'Blues'))
sdcols <- sdcols.fn(ncolours) 
map.stddev <- mapview(Pred, zcol = c("stddev"), legend = TRUE, alpha=0.3, 
                      col.regions = sdcols, layer.name = "SD pred")


elevation <- covariates[,-(2)]
elev.fn <- colorRampPalette(brewer.pal(9, 'Spectral'))
elevcol <- elev.fn(ncolours)
map.elev <- mapview(elevation, legend = TRUE, col.regions = elevcol, layer.name = "Elevation")

canopy <- covariates[,-(1)]
can.fn <- colorRampPalette(brewer.pal(9, 'BuGn'))
cancol <- can.fn(ncolours)
map.can <- mapview(canopy, legend = TRUE, col.regions = cancol, layer.name = "Canopy")


map.bbs <- mapview(BBS_sp, legend = TRUE, layer.name = "BBS", cex = 4, lwd = 0.5, col.regions = "yellow")
map.gbif <- mapview(GBIF_sp, legend = TRUE, layer.name = "GBIF", cex = 4, lwd = 0.5, col.regions = "pink")

BBAabs <- BBA_sp[which(BBA_sp@data$present == F),]
BBApres <- BBA_sp[which(BBA_sp@data$present == T),]
map.bba1 <- mapview(BBAabs, legend = TRUE, layer.name = "BBA absent", cex = 4, lwd = 0.5, col.regions = "grey")
map.bba2 <- mapview(BBApres, legend = TRUE, layer.name = "BBA present", cex = 4, lwd = 0.5, col.regions = "green")

map.mean + map.bbs + map.gbif + map.bba1 + map.bba2
map.mean + map.stddev + map.elev + map.can + map.bbs + map.gbif + map.bba1 + map.bba2
