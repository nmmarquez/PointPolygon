#' Take Observations from Polygons within the latent field
#'
#' @description Takes binomial trials from a probability field created by 
#' simField randomly from the field where each observation takes place within
#' a given polygon and is representative of that polygon.
#' 
#' @param field field object which simulated underlying data
#' @param N int, number of points for each polygon
#' @param M int, Number of trials for each point
#' @param times vector of int, times to sample from
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
#' samplePolygons(unitSim, N=20, M=100, rWidth=3)
#' 
#' @export

samplePolygons <- function(
    field, N, M, times=NULL, p=1., polygonList=NULL, rWidth=NULL, ...){
    if(is.null(polygonList)){
        if(is.null(rWidth)){
            stop("rWidth and polygonList can not both be NULL")
        }
        sectionedSPDF <- dividePolygon(field$bound, rWidth)
        polygonList <- lapply(1:nrow(sectionedSPDF@data), function(i){
            sf::st_as_sf(sectionedSPDF[i,])
        })
    }
    polyN <- length(polygonList)
    polySamples <- sort(sample.int(polyN, round(p*polyN), replace=F))
    sampleN <- length(polySamples)

    obsDF <- as.data.frame(dplyr::bind_rows(lapply(polySamples, function(j){
        subSPDF <- polygonList[[j]]
        subSPDF$isPresent <- TRUE
        pointsDF <- cbind(
            isPresent=sapply(sf::st_intersects(field$spdf,subSPDF), function(z){
                length(z)>0}), 
            field$spdf) %>%
            filter(isPresent)

        miniField <- list(spdf=pointsDF)

        samplePoints(miniField, N, M, times=times) %>%
            dplyr::rename(trueid=id) %>%
            dplyr::mutate(polyid=subSPDF$polyid) %>%
            dplyr::as_tibble() %>%
            dplyr::mutate(id=lapply(1:N, function(x) pointsDF$id)) %>%
            dplyr::select(id, tidx, trials, obs, polyid, trueid)
    })))

    obsDF
}

