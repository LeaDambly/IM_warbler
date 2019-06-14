
# modified function to take PC priors
FitModel <- function (..., formula = NULL, CovNames = NULL, mesh, spat.ind = "i", 
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
