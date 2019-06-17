# This is intended to fit a model with a random field for the eBird observation level.At the moment it hasn't got very far.

library(raster)
# library(rgdal)
# library(RColorBrewer)
library(rgeos)
library(INLA)
# library(mapview)

#source("warblerfunctions.R")
source("Functions/FitModel.R")
# load("BTWarblerData.RData")
load("Data/Stacks.RData")


  # create priors 
  C.F. <- list(mean = list(int.BBS = 0, int.eBird = 0, int.BBA = 0),
               mean.intercept=0, prec.intercept = 0.001,
               prec = list(int.BBS = 1, int.eBird = 1, int.BBA = 1))
  
  # The parameters "prior.range" and "prior.sd" control the joint prior on range and standard deviation of the spatial field.
  # see here for details: 
  # https://groups.google.com/d/msg/r-inla-discussion-group/dunoXK_yAco/JhmYb5JoAQAJ
  # e.g.: c(0.01, 0.05) means a 5% that the range will be less than 0.01
  # specific notes: "it is advisable not to favour ranges that are smaller than the resolution of the mesh"
  # "At this point I think one has to do some experimenting with the priors. As far as I know, we still do not have enough experience with the priors to come up with clearer guidelines"
  
  spde <- inla.spde2.pcmatern(mesh = Mesh$mesh, alpha = 2, 
                              prior.range = c(0.02, 0.5),
                              prior.sigma = c(5, 0.1))
  
  form <- formula(resp ~ 0 + elevation + canopy + Intercept + X + Y + int.BBS + density + 
                    int.eBird + int.BBA + f(i, model = spde) + f(e, model = spde))
  
  warbler_model <- FitModel(stk.ip, stk.pred$stk, stk.eBird, stk.BBS, stk.BBA,
                             formula = form, CovNames = NULL, mesh = Mesh$mesh,
                             predictions = TRUE, control.fixed = C.F., waic = TRUE, nthreads=16)
  
  save(C.F., spde, form, warbler_model, file="Data/WarblerModelOutputDensity.RData")
 
 q("no")
