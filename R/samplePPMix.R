#' Take Observations from Points and Polygons within the latent field
#'
#' @description Takes binomial trials from a probability field created by 
#' simField randomly from the field where each observation takes place at the 
#' point level and then p% polygons are removed of their point data and
#' replaced with representative polygon data.
#' 
#' @param field field object which simulated underlying data
#' @param N int, Number of points sampled
#' @param M int, Number of trials for each point
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
#' mixSample <- samplePPMix(unitSim, 30, 150, .5, rWidth=3)
#' head(mixSample$pointDF)
#' head(mixSample$polyDF)
#' 
#' @export

samplePPMix <- function(
    field, N, M, p=.5, polygonList=NULL, rWidth=NULL, replace=TRUE, ...){
    if(is.null(polygonList)){
        sectionedSPDF <- dividePolygon(field$bound, rWidth)
        polygonList <- lapply(1:nrow(sectionedSPDF@data), function(i){
            sectionedSPDF[i,]
        })
    }

    DF <- dplyr::sample_n(field$spdf@data, N, replace=replace)
    DF$trials <- M
    DF$obs <- stats::rbinom(N, size=DF$trials, prob=DF$theta)
    DF <- DF[,c("id", "trials", "obs")]

    polyN <- length(polygonList)
    polySamples <- sort(sample.int(polyN, round(p*polyN), replace=F))
    sampleN <- length(polySamples)

    obsDF <- data.frame(
        id = I(lapply(1:sampleN, function(x) NA)),
        trials = rep(0, sampleN),
        obs = 0,
        polyid = NA)
    
    i <- 0
    for(j in polySamples){
        i <- i + 1
        subSPDF <- polygonList[[j]]
        subSPDF$isPresent <- TRUE
        pointsDF <- cbind(sp::over(field$spdf, subSPDF), field$spdf@data)
        pointsDF <- pointsDF[!is.na(pointsDF$isPresent),]
        obsDF$id[[i]] <- list(ids=pointsDF$id)
        pRemoveDF <- subset(DF, id %in% pointsDF$id)
        mPolyTotal <- sum(pRemoveDF$trials)
        obsDF$trials[i] <- mPolyTotal
        obsDF$polyid[i] <- subSPDF$polyid
        obsDF$obs[i] <- sum(stats::rbinom(
            mPolyTotal, 1, sample(pointsDF$theta, mPolyTotal, replace=T)))
        DF <- subset(DF, !(id %in% pointsDF$id))
    }

    obsDF <- subset(obsDF, trials != 0)

    list(pointDF=DF, polyDF=obsDF)
}