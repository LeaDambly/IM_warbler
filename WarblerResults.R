library(raster)
library(sf)
library(USAboundaries)

# load results of all other scripts for plotting
load("Data/Stacks.RData")
load("Data/WarblerModelOutput.RData")
load("Data/BTWarblerData.RData")

Pred <- SpatialPixelsDataFrame(
  points = stk.pred$predcoords,
  data = warbler_model$predictions,
  proj4string = crs(proj)
)

Pred@data$precision <- Pred@data$stddev ^ -2


# Pennsylvania border
PA <- us_states(states = "Pennsylvania")
PA <- PA$geometry[1]
PA <- as(PA, "Spatial")

# Results

warbler_model$model$summary.fixed

par(
  mfrow = c(1, 2),
  mar = c(2.1, 2, 3, 1),
  oma = c(2, 2, 0, 0)
)

plot(
  warbler_model$model$marginals.fixed$elevation,
  type = "l",
  main = "Elevation",
  xlab = ""
)

abline(v = 0, col = "red")

plot(
  warbler_model$model$marginals.fixed$canopy,
  type = "l",
  main = "Canopy",
  xlab = ""
)

abline(v = 0, col = "red")

mtext("Estimate", 1, outer = TRUE, line = 0.2)
mtext("Density", 2, outer = TRUE, line = 0.7)


# Plot of data

png("DataPlot.png", height = 480, width = 640)

par(
  mar = rep(0.1, 4),
  mfrow = c(1, 1)
)

plot(PA)

# Add BBA
points(
  BBA_sp[which(BBA_sp@data$present == FALSE), ],
  cex = 0.4,
  col = "pink",
  pch = 16
)
points(
  BBA_sp[which(BBA_sp@data$present == TRUE), ],
  cex = 0.4,
  col = "red3",
  pch = 16
)

# Add BBS
points(
  BBS_sp,
  cex = 0.5,
  col = "blue",
  pch = 16
)

# Add eBird
points(
  eBird_sp,
  cex = 0.5,
  col = "sandybrown",
  pch = 16
)

legend(
  x = -79, y = 42.45,
  legend = c("BBA, absent", "BBA, present", "BBS route", "eBird"),
  col = c("pink", "red3", "blue", "sandybrown"),
  pch = 16,
  ncol = 2
)

dev.off()


# Plot of predictions

png("PredPlot.png", height = 360, width = 640)
par(mar = rep(0.1, 4))
plot(Pred, attr = "mean", col = grey(seq(0, 1, length = 100)))
lines(PA)
dev.off()


# Plot of uncertainty
png("StdDevPlot.png", height = 360, width = 640)
par(mar = rep(0.1, 4))
# plot(Pred, attr=3, col=0)
plot(Pred, attr = "stddev", col = grey(seq(0, 1, length = 100)))
lines(PA)
dev.off()

GreyCol <- function(x) {
  x.p <- (x - min(x)) / (max(x) - min(x))
  grey(x.p)
}

par(mfrow = c(1, 1))
plot(Pred@coords, col = GreyCol(Pred@data$precision), pch = 16)