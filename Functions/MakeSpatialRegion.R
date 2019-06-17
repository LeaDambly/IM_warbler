# just took PointedSDM's function apart a bit because buffer function used in MakeSpatialRegion doesn't work with unprojected data - so no buffering here
MakeSpatialRegion <- function (data = NULL, coords = c("X", "Y"), meshpars, bdry = NULL, 
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
