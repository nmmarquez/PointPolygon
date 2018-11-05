#' Simulate A continous 2D spatial field using SPDE
#'
#' @description Calculates a continous underlying spatial field for a given
#' shape and uses an additive model to combine with covariates to create a final
#' observed field which is a linear combination of covariates their beta values
#' with an optional link function applied.
#' 
#' @param N Number of width pixels for the projection
#' @param sigmaE spatial variance
#' @param rangeE spatial range
#' @param shape polygon object to populate, if NULL unit square is used
#' @param beta0 itercept term added to all points
#' @param betaList list of 2 item lists with a type item which is either random,
#' cluster, or spatial and a value item which is the beta coefficient.
#' @param link link function to apply to the linear combination
#' @param ... other parameters to pass to mesh
#'
#' @return field object, contains 6 items. Spatial points data frame with raster
#' values for transformed linear combination, and beta value. A mesh that was 
#' used to create the latent field and possibly covariates. The latent field
#' itself. A bounding shape where all observations take place. A projection 
#' matrix from the mesh to the entire field of interest. The spde for the matern
#' approximation.
#'
#' @examples
#' ## Not run:
#' require(tidyr)
#' require(gridExtra)
#' require(ar.matrix)
#' require(dplyr)
#' require(ggplot2)
#' 
#' unitSim <- simField(
#'     N = 500, 
#'     offset = c(0.1, 0.2), 
#'     max.edge = c(0.05,0.2),
#'     beta0 = -2,
#'     betaList = list(
#'         list(type="random", value=2),
#'         list(type="spatial", value=-.5),
#'         list(type="cluster", value=-2)
#'     ))
#' 
#' plotList <- lapply(c("V1", "V2", "V3", "theta"), function(eff){
#'     unitSim$spdf@data %>%
#'         dplyr::select(-V0, -z) %>%
#'         gather("Effect", "Value", V1:theta) %>%
#'         filter(Effect==eff) %>%
#'         ggplot(aes(x, y, fill=Value)) +
#'         geom_raster() +
#'         coord_equal() +
#'         theme_void() +
#'         scale_fill_distiller(palette = "Spectral") +
#'         ggtitle(eff)
#' })
#' 
#' do.call(grid.arrange, c(plotList, ncol=2))
#' 
#' usSim <- simField(
#'     N = 600,
#'     shape = US.df,
#'     rangeE = 1.7,
#'     offset = c(1, 2),
#'     max.edge = c(.5, 1))
#' 
#' plot(usSim$mesh)
#' usSim$spdf@data  %>%
#'     ggplot(aes(x, y, fill=z)) +
#'     geom_raster() +
#'     coord_equal() +
#'     theme_void() +
#'     scale_fill_distiller(palette = "Spectral")
#'
#' ## End(**Not run**)
#'
#' @export

simField <- function(N=60, sigmaE=1, rangeE=.3, shape=NULL, 
                     beta0=0, betaList=list(), link=arm::invlogit, ...){
    if(is.null(shape)){
        shape <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(
            matrix(c(0,1,1,0,0,0,1,1), ncol=2))), 1)))
    }
    bb <- shape@bbox
    shape$isPresent <- TRUE
    baseRaster <- raster::raster(
        ncols=N, nrows=round(N*hRatio(shape)), 
        xmn=bb[1,1], xmx=bb[1,2], ymn=bb[2,1], ymx=bb[2,2], 
        crs=shape@proj4string)
    shapeRaster <- raster::rasterize(shape, baseRaster, field=1)
    shapeExtPointsDF <- sp::SpatialPointsDataFrame(
        sp::coordinates(shapeRaster),
        data=as.data.frame(sp::coordinates(shapeRaster)),
        proj4string = shape@proj4string)
    
    validIndex <- !is.na(sp::over(shapeExtPointsDF, shape)$isPresent)
    shapePointsDF <- shapeExtPointsDF[validIndex,]
    mesh <- INLA::inla.mesh.2d(loc=bbCoords(shape), ...)
    AprojField <- INLA::inla.spde.make.A(mesh=mesh, loc=shapePointsDF)
    
    kappaE <- sqrt(8) / rangeE
    tauE <- 1/(sqrt(4*pi)*kappaE*sigmaE)
    spde <- INLA::inla.spde2.matern(mesh)
    Q <- tauE**2 * (kappaE**4 * spde$param.inla$M0 +
                        2 * kappaE**2 *spde$param.inla$M1 + spde$param.inla$M2)
    x_ <- as.vector(ar.matrix::sim.AR(n=1, Q))
    x <- x_ - mean(x_)
    N <- nrow(shapePointsDF@data)
    betaVec <- beta0
    covMat <- matrix(1, nrow=N)
    shapePointsDF$V0 <- 1
    i <- 0
    for(covb in betaList){
        i <- i + 1
        betaVec <- c(betaVec, covb$value)
        if(covb$type == "spatial"){
            cov_ <- as.vector(ar.matrix::sim.AR(n=1, Q))
            cov <- cov_ - mean(cov_)
            newCov <- as.numeric(AprojField %*% cov)
        }
        else if(covb$type == "cluster"){
            suppressWarnings(clusters <- stats::kmeans(
                shapePointsDF@coords, centers = 10, iter.max = 1))
            newCov <- stats::runif(10)[clusters$cluster]
        }
        else{
            newCov <- stats::runif(N)
        }
        shapePointsDF@data[,paste0("V", i)] <- newCov
    }
    shapePointsDF$z <- as.numeric(AprojField %*% x)
    if(i == 0){
        shapePointsDF$theta <- link(beta0 + shapePointsDF$z)
    }
    else{
        shapePointsDF$theta <- link(c(
            as.matrix(shapePointsDF@data[,paste0("V", 0:i)]) %*% betaVec) + 
                shapePointsDF$z)
    }
    
    shapePointsDF$id <- (1:nrow(shapePointsDF@data)) - 1
    
    shapePointsDF$Bound <- 1
    boundShape <- rgeos::gUnaryUnion(shape, id=shape@data$Bound)
    boundShape <- sp::SpatialPolygonsDataFrame(boundShape,
            data.frame(
                Bound=1,
                row.names="1"))
    
    field <- list(
        spdf = shapePointsDF, 
        mesh = mesh, 
        latent = x, 
        bound = boundShape,
        AprojField = AprojField,
        spde = spde)
    class(field) <- "field"
    field
}