#' Run a Model and Estimate Underlying Field From Point and Polygon Data
#'
#' @description Run a model and estimate underlying latent and probability
#' field from point and polygon data using TMB. 
#' 
#' @param field field object which simulated underlying data
#' @param plotVar character, which variable to plot
#' 
#' @return ggplot object of underlying field
#'
#' @examples 
#' unitSim <- simField(
#' N = 100, rangeE = .5,
#' offset = c(0.1, 0.2), 
#' max.edge = c(0.1,0.2),
#' beta0 = -2)
#' 
#' ggField(unitSim)
#'   
#' @export

ggField <- function(field, plotVar="theta"){
    ggplot2::ggplot(
        field$spdf, 
        ggplot2::aes_string("x", "y", fill=plotVar)) +
        ggplot2::geom_raster() +
        ggplot2::coord_equal() +
        ggplot2::theme_void() +
        ggplot2::scale_fill_distiller(palette = "Spectral") +
        ggplot2::facet_wrap(~tidx)
}
