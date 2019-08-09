#' Find index of closest neighbor
#'
#' cpp code for finding closest neighbor
#'
#' @param x1 matrix of long lats for desired match
#' @param x2 matrix of long lats for matching to
#' @param fullDF DF indicating the full id names
#' @return integer vector indicating closest neighbor
#' @export
#' @useDynLib PointPolygon

find_closest <- function(x1, x2, fullDF) {
  
  toRad <- pi / 180
  lat1  <- x1[,2]  * toRad
  long1 <- x1[,1] * toRad
  lat2  <- x2[,2]  * toRad
  long2 <- x2[,1] * toRad
  
  ord1  <- order(lat1)
  rank1 <- match(seq_along(lat1), ord1)
  ord2  <- order(lat2)
  
  ind <- find_closest_point(lat1[ord1], long1[ord1], lat2[ord2], long2[ord2])
  
  fullDF$id[ord2[ind + 1][rank1]]
}