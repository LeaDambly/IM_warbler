library(raster)
# library(rgdal)
# library(RColorBrewer)
library(rgeos)
library(PointedSDMs)
# library(mapview)

source("warblerfunctions.R")
# load("BTWarblerData.RData")
load("Data/Stacks.RData")


  # create priors 
  C.F. <- list(mean = list(int.BBS = 0),
               mean.intercept=0, prec.intercept = 0.001,
               prec = list(int.BBS = 1))
  
  # The parameters "prior.range" and "prior.sd" control the joint prior on range and standard deviation of the spatial field.
  # see here for details: 
  # https://groups.google.com/d/msg/r-inla-discussion-group/dunoXK_yAco/JhmYb5JoAQAJ
  # e.g.: c(0.01, 0.05) means a 5% that the range will be less than 0.01
  # specific notes: "it is advisable not to favour ranges that are smaller than the resolution of the mesh"
  # "At this point I think one has to do some experimenting with the priors. As far as I know, we still do not have enough experience with the priors to come up with clearer guidelines"
  
  spde <- inla.spde2.pcmatern(mesh = Mesh$mesh, alpha = 2, 
                              prior.range = c(0.02, 0.5),
                              prior.sigma = c(5, 0.1))
  
  form <- formula(resp ~ 0 + elevation + canopy + Intercept + X + Y + int.BBS + 
                    f(i, model = spde))
  
  warbler_model_BBS <- FitModel2(stk.ip, stk.pred$stk, stk.BBS,
                                 formula = form, CovNames = NULL, mesh = Mesh$mesh,
                                 predictions = TRUE, control.fixed = C.F., waic = TRUE, nthreads=8)
  
  save(C.F., spde, form, warbler_model_BBS, file="Data/WarblerModelOutputBBS.RData")
  
q("no")
