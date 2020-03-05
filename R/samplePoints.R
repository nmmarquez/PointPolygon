#' Take Point Observations from a spatial field.
#'
#' @description Takes binomial trials from a probability field created by 
#' simField randomly from points on the field.
#' 
#' @param field field object which simulated underlying data
#' @param N int 1, Number of points to simulate
#' @param M int, Number of trials for each point either 1 or N long 
#' @param replace logical, replace pixel probabilities in sampling?
#' @param times vector of int, times to sample from
#' @param ... further argumnets for compatibility 
#' 
#' @return data.frame with observation number, number of trials, and the id
#' of the pixel that the trail took place in.
#' 
#' @examples
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
#' samplePoints(unitSim, N=20, M=100) 
#' 
#' @export

samplePoints <- function(field, N, M, replace=TRUE, times=NULL, ...){
    if(is.null(times)){
        times <- unique(field$spdf$tidx)
    }
    DF <- field$spdf %>%
        dplyr::as_tibble() %>%
        dplyr::filter(sapply(tidx, function(t) t %in% times)) %>%
        dplyr::sample_n(N, replace=replace)
    DF$trials <- M
    DF$obs <- stats::rbinom(N, size=DF$trials, prob=DF$theta)
    DF[,c("id", "tidx", "trials", "obs")]
}
