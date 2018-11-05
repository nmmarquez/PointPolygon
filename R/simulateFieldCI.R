#' Estimate Confidence Intervals (CI) for Probability Field Via Simulation
#'
#' @description Calculate the mean, 95% confidence intervals, and the standard
#' deviation of an estimated probability field from a model. 
#' 
#' @param field field object which simulated underlying data
#' @param modelFit fitted model object from the field
#' @param draws int, number of draws for simulation to calculate CI's
#' 
#' @return data.frame, mean, CI, and sd of underlying field.
#' 
#' ## Not run:
#' 
#' unitSim <- simField(
#' N = 200, rangeE = .7,
#' offset = c(0.1, 0.2), 
#' max.edge = c(0.1,0.2),
#' beta0 = -2,
#' betaList = list(
#'     list(type="random", value=2),
#'     list(type="spatial", value=-.5),
#'     list(type="cluster", value=-2)
#' ))
#'
#' pointDF <- samplePoints(unitSim, 500, 100)
#'
#' unitFit <- runFieldModel(unitSim, pointDF)
#'
#' head(simulateFieldCI(unitSim, unitFit))
#'
#' ## End(**Not run**)
#' 
#' @export

simulateFieldCI <- function(field, modelFit, draws=1000){
    modelLRP <- buildModelInputs(field, model=FALSE)
    
    parDraws <- t(ar.matrix::sim.AR(draws, modelFit$sd$jointPrecision)) +
        c(modelFit$opt$par, modelFit$sd$par.random)
    
    betaRows <- row.names(modelFit$sd$jointPrecision) == "beta"
    zRows <- row.names(modelFit$sd$jointPrecision) == "z"
    
    betaDraws <- parDraws[betaRows,]
    zDraws <- parDraws[zRows,]
    
    fieldProbs <- arm::invlogit(
        modelLRP$Data$covs %*% betaDraws + modelLRP$Data$AprojObs %*% zDraws)
    
    data.frame(
        mu = apply(fieldProbs, 1, mean),
        sd = apply(fieldProbs, 1, stats::sd),
        lwr = apply(fieldProbs, 1, stats::quantile, probs=.025),
        upr = apply(fieldProbs, 1, stats::quantile, probs=.975),
        id = field$spdf$id
        )
}
