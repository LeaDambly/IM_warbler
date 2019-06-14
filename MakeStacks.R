library(INLA)

load("Data/BTWarblerData.RData")
source("Functions/MakeSpatialRegion.R")
source("Functions/MakeIntegrationStack.R")
source("Functions/GetNearestCovariate.R")
source("Functions/MakeProjectionGrid.R")

source("Functions/MakeBinomStack.R")
source("Functions/MakePointsStack.R")
# source("warblerfunctions.R")

# Create the mesh to approximate the area and the spatial field
Meshpars <- list(max.edge = c(0.05, 0.4), offset = c(0.1, 0.4), cutoff = 0.1)
Mesh <- MakeSpatialRegion(data = NULL, bdry = PA, meshpars = Meshpars,
                           proj =  proj)
# plot(Mesh$mesh)

# Make stack for background mesh
stk.ip <- MakeIntegrationStack(mesh = Mesh$mesh, data = covariates_eBird, area=Mesh$w, 
                                   tag='ip', InclCoords=TRUE)

# make data for projections
Nxy.scale <- 0.01 # change the resolution of the predictions
Boundary <- Mesh$mesh$loc[Mesh$mesh$segm$int$idx[,2],]
Nxy <- round(c(diff(range(Boundary[,1])), diff(range(Boundary[,2])))/Nxy.scale)

# Make stack for projections
stk.pred <- MakeProjectionGrid(nxy = Nxy, mesh = Mesh$mesh, data = covariates,
                               tag = 'pred', boundary = Boundary)
# Create data stacks
stk.BBS <- MakeBinomStack(observs = BBS_sp,data = covariates, mesh=Mesh$mesh,
                          presname="NPres", trialname="Ntrials", tag='BBS', InclCoords=TRUE)
stk.eBird <- MakePointsStack(presences = eBird_sp, data = covariates_eBird, mesh = Mesh$mesh, 
                            tag = 'eBird', InclCoords = TRUE)
stk.BBA <- MakeBinomStack(observs = BBA_sp, data = covariates, mesh = Mesh$mesh,
                          presname='present', tag='BBA', InclCoords=TRUE)

save(stk.ip, stk.pred, stk.eBird, stk.BBS, stk.BBA, Mesh, file="Data/Stacks.RData")
