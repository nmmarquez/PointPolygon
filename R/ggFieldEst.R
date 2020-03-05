#' Run a Model and Estimate Underlying Field From Point and Polygon Data
#'
#' @description Run a model and estimate underlying latent and probability
#' field from point and polygon data using TMB. 
#' 
#' @param field field object which simulated underlying data
#' @param predList list of model predictions
#' @param sd boolean, plot the standard error of field estimates
#' 
#' @return ggplot object of underlying field and model estimates
#'
#' @examples
#' \dontrun{
#' unitSim <- simField(
#' N = 100, rangeE = .5,
#' offset = c(0.1, 0.2), 
#' max.edge = c(0.1,0.2),
#' beta0 = -2)
#' 
#' pointDF1 <- samplePoints(unitSim, N=1200, M=100) 
#' pointDF2 <- samplePoints(unitSim, N=1200, M=100)
#' 
#' modelPredList <- lapply(list(pointDF1, pointDF2), function(x){
#'     simulateFieldCI(unitSim, runFieldModel(unitSim, x))
#' })
#' 
#' ggFieldEst(unitSim, modelPredList)
#' }
#' 
#' @export

ggFieldEst <- function(field, predList, sd=FALSE){
    DF <- dplyr::select(dplyr::as_tibble(field$spdf), -geometry)
    allDF <- DF
    allDF$Type <- "True"
    allDF$sd <- NA
    if(is.null(names(predList))){
        names(predList) <- 1:length(predList)
    }
    for(i in 1:length(predList)){
        DFpred <- DF
        DFpred$Type <- names(predList)[i]
        DFpred$theta <- predList[[i]]$mu
        DFpred$sd <- predList[[i]]$sd
        allDF <- dplyr::bind_rows(allDF, DFpred)
    }
    
    new <- field
    new$spdf <- allDF
    if(sd){
        new$spdf <- allDF[allDF$Type != "True",]
        return(ggField(new, "sd") + ggplot2::facet_grid(Type~tidx))
    }
    ggField(new) +
        ggplot2::facet_grid(Type~tidx)
}
