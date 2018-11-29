#' Take Observations from Polygons within the latent field
#'
#' @description Takes binomial trials from a probability field created by 
#' simField randomly from the field where each observation takes place within
#' a given polygon and is representative of that polygon.
#' 
#' @param field field object which simulated underlying data
#' @param M int, Number of trials for each polygon
#' @param p numeric >0 & <=1, percentage of polygons sampled
#' @param polygonList list, list of polygons to sample from
#' @param rWidth integer, instead of using a polygon list divide the original
#' bounding box into approximately equal sized squares with rWidth number 
#' squares along the x axis.
#' @param ... further argumnets for compatibility 
#' 
#' @return data.frame with observation number, number of trials, and the id
#' of the pixel that the trail took place in.
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
#' samplePolygons(unitSim, round(1200*100/9), rWidth=3)
#' 
#' @export

samplePolygons <- function(field, M, p=1., polygonList=NULL, rWidth=NULL, ...){
    if(is.null(polygonList)){
        sectionedSPDF <- dividePolygon(field$bound, rWidth)
        polygonList <- lapply(1:nrow(sectionedSPDF@data), function(i){
            sectionedSPDF[i,]
        })
    }
    polyN <- length(polygonList)
    polySamples <- sort(sample.int(polyN, round(p*polyN), replace=F))
    sampleN <- length(polySamples)

    obsDF <- data.frame(
        id = I(lapply(1:sampleN, function(x) NA)),
        trials = rep(M, sampleN),
        obs = NA, 
        polyid = NA)

    i <- 0
    for(j in polySamples){
        i <- i + 1
        subSPDF <- polygonList[[j]]
        subSPDF$isPresent <- TRUE
        pointsDF <- cbind(sp::over(field$spdf, subSPDF), field$spdf@data)
        pointsDF <- pointsDF[!is.na(pointsDF$isPresent),]
        obsDF$polyid[i] <- subSPDF$polyid
        obsDF$id[[i]] <- list(ids=pointsDF$id)
        obsDF$obs[i] <- sum(stats::rbinom(
            M, 1, sample(pointsDF$theta, M, replace=T)))
    }

    obsDF
}