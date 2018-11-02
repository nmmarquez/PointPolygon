rm(list=ls())
library(INLA)
library(ar.matrix)
library(dplyr)
library(ggplot2)
library(rgdal)
library(raster)
set.seed(123)

bbCoords <- function(SpatPolygon){
    bb <- SpatPolygon@bbox 
    SpatialPoints(
        as.matrix(expand.grid(x=bb[1,], y=bb[2,])),
        proj4string = SpatPolygon@proj4string)
}

hRatio <- function(SpatPolygon){
    dims <- unname(apply(SpatPolygon@bbox, 1, diff))
    dims[2]/dims[1]
}

simField <- function(N=60, sigmaE=1, rangeE=.3, shape=NULL, ...){
    if(is.null(shape)){
        shape <- SpatialPolygons(list(Polygons(list(Polygon(
            matrix(c(0,1,1,0,0,0,1,1), ncol=2))), 1)))
    }
    bb <- shape@bbox
    shape$isPresent <- TRUE
    baseRaster <- raster(
        ncols=N, nrows=round(N*hRatio(shape)), 
        xmn=bb[1,1], xmx=bb[1,2], ymn=bb[2,1], ymx=bb[2,2], 
        crs=shape@proj4string)
    shapeRaster <- rasterize(shape, baseRaster, field=1)
    shapeExtPointsDF <- SpatialPointsDataFrame(
        coordinates(shapeRaster),
        data=as.data.frame(coordinates(shapeRaster)),
        proj4string = shape@proj4string)
    
    validIndex <- !is.na(over(shapeExtPointsDF, shape)$isPresent)
    shapePointsDF <- shapeExtPointsDF[validIndex,]
    mesh <- inla.mesh.2d(
        loc=bbCoords(shape),
        ...)
    AprojField <- inla.spde.make.A(mesh=mesh, loc=shapePointsDF)
    
    kappaE <- sqrt(8) / rangeE
    tauE <- 1/(sqrt(4*pi)*kappaE*sigmaE)
    spde <- inla.spde2.matern(mesh)
    Q <- tauE**2 * (kappaE**4 * spde$param.inla$M0 +
                        2 * kappaE**2 *spde$param.inla$M1 + spde$param.inla$M2)
    x_ <- as.vector(sim.AR(n=1, Q))
    x <- x_ - mean(x_)
    
    shapePointsDF$z <- as.numeric(AprojField %*% x)
    list(spdf=shapePointsDF, mesh=mesh, latent=x)
}

simField <- function(N=60, sigmaE=1, rangeE=.3,
                     shape=NULL, beta0=0, betaList=list(), link=identity, ...){
    if(is.null(shape)){
        shape <- SpatialPolygons(list(Polygons(list(Polygon(
            matrix(c(0,1,1,0,0,0,1,1), ncol=2))), 1)))
    }
    bb <- shape@bbox
    shape$isPresent <- TRUE
    baseRaster <- raster(
        ncols=N, nrows=round(N*hRatio(shape)), 
        xmn=bb[1,1], xmx=bb[1,2], ymn=bb[2,1], ymx=bb[2,2], 
        crs=shape@proj4string)
    shapeRaster <- rasterize(shape, baseRaster, field=1)
    shapeExtPointsDF <- SpatialPointsDataFrame(
        coordinates(shapeRaster),
        data=as.data.frame(coordinates(shapeRaster)),
        proj4string = shape@proj4string)
    
    validIndex <- !is.na(over(shapeExtPointsDF, shape)$isPresent)
    shapePointsDF <- shapeExtPointsDF[validIndex,]
    mesh <- inla.mesh.2d(
        loc=bbCoords(shape),
        ...)
    AprojField <- inla.spde.make.A(mesh=mesh, loc=shapePointsDF)
    
    kappaE <- sqrt(8) / rangeE
    tauE <- 1/(sqrt(4*pi)*kappaE*sigmaE)
    spde <- inla.spde2.matern(mesh)
    Q <- tauE**2 * (kappaE**4 * spde$param.inla$M0 +
                        2 * kappaE**2 *spde$param.inla$M1 + spde$param.inla$M2)
    x_ <- as.vector(sim.AR(n=1, Q))
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
            cov_ <- as.vector(sim.AR(n=1, Q))
            cov <- cov_ - mean(cov_)
            newCov <- as.numeric(AprojField %*% cov)
        }
        else if(covb$type == "cluster"){
            suppressWarnings(clusters <- kmeans(
                shapePointsDF@coords, centers = 10, iter.max = 1))
            newCov <- runif(10)[clusters$cluster]
        }
        else{
            newCov <- runif(N)
        }
        shapePointsDF@data[,paste0("V", i)] <- newCov
    }
    shapePointsDF$z <- as.numeric(AprojField %*% x)
    if(i == 0){
        shapePointsDF$theta <- beta0 + shapePointsDF$z
    }
    else{
        shapePointsDF$theta <- link(c(
            as.matrix(shapePointsDF@data[,paste0("V", 0:i)]) %*% betaVec) + 
                shapePointsDF$z)
    }
    
    list(spdf=shapePointsDF, mesh=mesh, latent=x)
}

unitSim <- simField(
    N = 500, 
    offset = c(0.1, 0.2), 
    max.edge = c(0.05,0.2),
    beta0 = -2,
    betaList = list(
        list(type="random", value=2),
        list(type="spatial", value=-3),
        list(type="cluster", value=1)
    ))

unitSim$spdf@data %>%
    ggplot(aes(x, y, fill=theta)) +
    geom_raster() +
    coord_equal() +
    theme_void() +
    scale_fill_distiller(palette = "Spectral")

unitSim$spdf@data %>%
    ggplot(aes(x, y, fill=V1)) +
    geom_raster() +
    coord_equal() +
    theme_void() +
    scale_fill_distiller(palette = "Spectral")

unitSim$spdf@data %>%
    ggplot(aes(x, y, fill=V2)) +
    geom_raster() +
    coord_equal() +
    theme_void() +
    scale_fill_distiller(palette = "Spectral")

unitSim$spdf@data %>%
    ggplot(aes(x, y, fill=V3)) +
    geom_raster() +
    coord_equal() +
    theme_void() +
    scale_fill_distiller(palette = "Spectral")

unitSim$spdf@data %>%
    ggplot(aes(x, y, fill=z)) +
    geom_raster() +
    coord_equal() +
    theme_void() +
    scale_fill_distiller(palette = "Spectral")

usSim <- simField(
    N = 600,
    shape = US.df,
    rangeE = 3,
    offset = c(1, 2),
    max.edge = c(.25, 1))


plot(usSim$mesh)
usSim$spdf@data  %>%
    ggplot(aes(x, y, fill=z)) +
    geom_raster() +
    coord_equal() +
    theme_void() +
    scale_fill_distiller(palette = "Spectral")
