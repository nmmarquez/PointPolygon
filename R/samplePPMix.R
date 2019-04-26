#' Take Observations from Points and Polygons within the latent field
#'
#' @description Takes binomial trials from a probability field created by 
#' simField randomly from the field where each observation takes place at the 
#' point level and then p% polygons are removed of their point data and
#' replaced with representative polygon data.
#' 
#' @param field field object which simulated underlying data
#' @param N int, Number of points sampled for each polygon
#' @param M int, Number of trials for each point
#' @param times vector of int, times to sample from
#' @param p numeric >0 & <=1, percentage of polygons sampled
#' @param polygonList list, list of polygons to sample from
#' @param rWidth integer, instead of using a polygon list divide the original
#' bounding box into approximately equal sized squares with rWidth number 
#' squares along the x axis.
#' @param replace logical, replace values in the sampling
#' @param ... further argumnets for compatibility 
#' 
#' @return list of data.frames with observation number, number of trials, 
#' and the id of the pixel that the trail took place in.
#' 
#' @examples 
#' require(sp)
#' 
#' set.seed(123)
#' unitSim <- simField(
#'     N = 300, 
#'     offset = c(0.1, 0.2), 
#'     max.edge = c(0.05,0.2),
#'     beta0 = -2,
#'     betaList = list(
#'         list(type="random", value=2),
#'         list(type="spatial", value=-.5),
#'         list(type="cluster", value=-2)
#'     ))
#'
#' mixSample <- samplePPMix(unitSim, 30, 150, p=.5, rWidth=3)
#' head(mixSample$pointDF)
#' head(mixSample$polyDF)
#' 
#' @export

samplePPMix <- function(
    field, N, M, times=NULL, p=.5, polygonList=NULL, rWidth=NULL, replace=TRUE, ...){

    obsDF <- samplePolygons(
        field, N, M, times=times, p=1, polygonList, rWidth, replace, ...)
    polyIDS <- unique(obsDF$polyid)
    
    polySamples <- sample(polyIDS, round(p*length(polyIDS)))
    pointSamples <- setdiff(polyIDS, polySamples)
    polyDF <- obsDF[obsDF$polyid %in% polySamples,]
    pointDF <- obsDF[obsDF$polyid %in% pointSamples,]
    pointDF$id <- NULL
    pointDF$id <- pointDF$trueid
    pointDF$trueid <- NULL

    list(pointDF=pointDF, polyDF=polyDF)
}
