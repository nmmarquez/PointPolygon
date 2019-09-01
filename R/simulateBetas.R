#' Estimate Confidence Intervals (CI) for Beta Values in Linear Model
#'
#' @description Calculate the mean, 95% confidence intervals, and the standard
#' deviation of an estimated beta coefficient. 
#' 
#' @param field field object which simulated underlying data
#' @param modelFit fitted model object from the field
#' @param draws int, number of draws for simulation to calculate CI's
#' 
#' @return data.frame, mean, sd and true value of underlying field.
#' 
#' @examples
#' \dontrun{
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
#' head(simulateBetas(unitSim, unitFit))
#'
#' }
#' 
#' @export

simulateBetas <- function(field, modelFit, draws=1000){
    if(modelFit$moption != 2){
        blocs <- names(modelFit$opt$par) == "beta"
        betaSims <- ar.matrix::sim.AR(draws, modelFit$sd$jointPrecision)
        betaHat <- unname(modelFit$opt$par[blocs])
        betaStErr <- apply(betaSims[,1:sum(blocs)], 2, sd)
        betaDF <- data.frame(
            betaHat = betaHat,
            betaStErr = betaStErr,
            betaTV = field$betas)
    }
    else{
        betaDF <- modelFit$summary.fixed[,c("mean", "sd")]
        betaDF <- betaDF[startsWith(row.names(betaDF), "V"),]
        row.names(betaDF) <- NULL
        names(betaDF) <- c("betaHat", "betaStErr")
        betaDF$betaTV <- field$betas
    }
    return(betaDF)
}
