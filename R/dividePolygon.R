#' TDivide field into sections
#'
#' @description Takes a field object and divides it into rWidth square sections.
#' 
#' @param shape spatial polygons data frame
#' @param rWidth integer, instead of using a polygon list divide the original
#' bounding box into approximately equal sized squares with rWidth number 
#' squares along the x axis.
#' 
#' @return spatial polygons data frame section into apprximately square sections
#' 
#' @examples 
#' require(sp)
#' require(ar.matrix)
#' 
#' set.seed(123)
#' usSim <- simField(
#'     N = 300,
#'     shape = US.df,
#'     rangeE = 3,
#'     offset = c(1, 2),
#'     max.edge = c(.25, 1))
#'
#' plot(dividePolygon(usSim$bound, 4))
#' plot(dividePolygon(usSim$bound, 10))
#' plot(dividePolygon(usSim$bound, 16))
#' 
#' @export

dividePolygon <- function(shape, rWidth){
    bb <- shape@bbox
    baseRaster <- raster::rasterToPolygons(raster::raster(
        ncols=rWidth, nrows=round(rWidth*hRatio(shape)), 
        xmn=bb[1,1], xmx=bb[1,2], ymn=bb[2,1], ymx=bb[2,2], 
        crs=shape@proj4string))
    sectionedSPDF <- raster::intersect(shape, baseRaster)
    sectionedSPDF$polyid <- (1:nrow(sectionedSPDF@data))-1
    sectionedSPDF
}