#' Calculate the height to width ratio of a spatial polygon
#'
#' @description Calculate the height to width ratio of a spatial polygon
#' 
#' @param SpatPolygon Spatial Polygons object
#'
#' @return numeric, ratio of diff
#'
#' @examples
#' require(ar.matrix)
#' hRatio(US.df)
#'
#' @export

hRatio <- function(SpatPolygon){
    dims <- unname(apply(SpatPolygon@bbox, 1, diff))
    dims[2]/dims[1]
}