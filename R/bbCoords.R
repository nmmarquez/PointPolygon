#' Extract Boundary of Polygon as a Spatial Points Object
#'
#' @description Extracts a bounding box object taken from a spatial polygons
#' data frame and returns a spatial points object with four coordinates
#' representing the bounding box.
#'
#' @param SpatPolygon Spatial Polygons object
#'
#' @return Spatial Points object of length 4
#'
#' @examples
#' require(ar.matrix)
#' bbCoords(US.df)
#'
#' @export

bbCoords <- function(SpatPolygon){
    bb <- SpatPolygon@bbox 
    sp::SpatialPoints(
        as.matrix(expand.grid(x=bb[1,], y=bb[2,])),
        proj4string = SpatPolygon@proj4string)
}