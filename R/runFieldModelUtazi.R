#' Run a Model and Estimate Underlying Field From Point and Polygon Data
#'
#' @description Run a model and estimate underlying latent and probability
#' field from point and polygon data using INLA and Utazi method.
#' 
#' @param field field object which simulated underlying data
#' @param pointDF data simulated from samplePoints
#' @param polyDF data simulated from samplePolygns
#' @param moption int, intger indicating how polygon data should be estimated 0
#' is by Reimann sum approximation, 1 is by redistribution, and 2 is by Utazi
#' approach.
#' @param verbose logical, print model fitting information
#' @param symbolic logical, use metas reordering in model fitting
#' @param control list, control list passed to nlminb
#' @param rWidth int, only used in moption 2 to build besag prior
#' 
#' @return List of fitted model objects.
#' 
#' @import sp
#' @import dplyr
#' @import INLA
#' 
#' @examples
#' \dontrun{
#' unitSim <- simField(
#' N = 500, rangeE = .7,
#' offset = c(0.1, 0.2), 
#' max.edge = c(0.1,0.2),
#' beta0 = -2,
#' betaList = list(
#'     list(type="random", value=2),
#'     list(type="spatial", value=-.5),
#'     list(type="cluster", value=-2)
#' ))
#'
#' polyDF <- samplePolygons(unitSim, round(1200*100/9), rWidth=3)
#'
#' runFieldModelUtazi(unitSim, polyDF=polyDF, moption=2, rWidth=3)
#'
#' }
#' 
#' @export

runFieldModelUtazi <- function(
    field, 
    pointDF = NULL,
    polyDF = NULL,
    moption = 0,
    verbose = FALSE,
    symbolic = TRUE,
    control = list(eval.max=1e4, iter.max=1e4),
    rWidth = NULL
    ){
    library(sp)
    library(dplyr)
    library(INLA)

    if(is.null(pointDF)){
        empty <- vector("integer")
        pointDF <- data.frame(obs=empty, trials=empty, id=empty)
    }
    
    polyDF <- polyDF %>% 
        select(-id, -trueid) %>%
        group_by(polyid) %>% 
        summarize_all(sum) %>%
        as.data.frame

    shape3 <- dividePolygon(field$bound, rWidth)
    shape3@data <- dplyr::left_join(shape3@data, polyDF, by="polyid")
    pointA <- field$AprojField[pointDF$id + 1,]
    
    grid.poly.no <- sp::over(as(field$spdf, "Spatial"), shape3)$polyid
    ran.pr <- as.numeric(summary(dist(sf::st_coordinates(field$spdf)))[3])
    kap.pr <- sqrt(8)/ran.pr
    spde.model <- INLA::inla.spde2.matern(
        mesh=field$mesh, alpha=2,
        B.tau=matrix(c(0, 1, 0),nrow=1,ncol=3),
        B.kappa=matrix(c(0, 0, 1),nrow=1,ncol=3),
        theta.prior.mean=c(0, log(kap.pr)),
        theta.prior.prec=c(1, 1))
    
    area.poly.nb <- spdep::poly2nb(shape3, snap = 1)
    area.poly.adj <- as(spdep::nb2mat(area.poly.nb, style = "B"), "dgTMatrix")
    
    # extract covariate values for polygons
    covPoly <- field$spdf %>%
        dplyr::as_tibble() %>%
        select(matches("V[0-9]+")) %>%
        mutate(polyid=grid.poly.no) %>%
        group_by(polyid) %>%
        summarise_all(mean, na.rm=T)
    
    covPoint <- field$spdf %>%
        dplyr::as_tibble() %>%
        dplyr::select(-geometry) %>%
        dplyr::right_join(pointDF) %>%
        dplyr::select(dplyr::matches("V[0-9]+")) %>%
        as.data.frame
    
    mesh.coord.in <-
        field$mesh$loc[as.vector(which(!is.na(over(
            SpatialPoints(field$mesh$loc), shape3)$polyid))),]
    
    mesh.coord.poly.no <- over(SpatialPoints(mesh.coord.in), shape3)$polyid
    
    areaA <- inla.spde.make.A(
        mesh=field$mesh,
        loc=mesh.coord.in,
        block=mesh.coord.poly.no+1, block.rescale="sum")
    
    stack.area <- inla.stack(
        tag='areal',
        data=list(y=shape3$obs, n=shape3$trials),
        # A maps how yo get from the effects estimated to Areal data observation
        # its a 1 for both non SPDE random effects and covariates.
        A=list(areaA,1,1),
        effects=list(
            s=1:field$spde$n.spde,
            sa=1:nrow(shape3@data),
            as.data.frame(select(covPoly, -polyid))))
    
    stack.pred <- inla.stack(
        tag='pred',
        data=list(
            y=rep(NA, nrow(field$spdf)),
            n=rep(NA, nrow(field$spdf))),
        A=list(field$AprojField,1,1),
        effects=list(
            s=1:field$spde$n.spde,
            sa=grid.poly.no+1,
            field$spdf %>%
                dplyr::as_tibble() %>%
                dplyr::select(dplyr::matches("V[0-9]+")) %>%
                as.data.frame))
    if(nrow(pointDF) != 0){
        pointAreas <- pointDF %>%
            dplyr::left_join(dplyr::as_tibble(field$spdf)) %>%
            select(x, y) %>%
            SpatialPoints %>%
            over(shape3) %>%
            select(polyid) %>%
            unlist %>%
            unname
        
        stack.points <- inla.stack(
            tag = 'points',
            data=list(y=pointDF$obs, n=pointDF$trials),
            A=list(pointA,1,1),
            effects=list(
                s=1:field$spde$n.spde,
                sa=pointAreas+1,
                covPoint))
        
        stack.all <- inla.stack(stack.points, stack.area, stack.pred)
    }
    
    else{
        stack.all <- inla.stack(stack.area, stack.pred)
    }
    f <- y ~ -1 + V0 + V1 + f(s, model=field$spde) + f(
        sa, 
        model = "besagproper2", 
        graph = area.poly.adj,
        hyper=list(prec=list(prior="loggamma", param=c(1,0.01))))
    
    startTime <- Sys.time()
    
    modfit <- inla(
        f, 
        data = inla.stack.data(stack.all),
        family = "binomial",
        Ntrials = stack.all$data$data$n,
        control.predictor=list(compute=TRUE, A=inla.stack.A(stack.all), link=1),
        control.compute=list(dic=TRUE, config=TRUE),
        verbose = verbose)
    
    runtime <- Sys.time() - startTime
    
    modfit$opt <- list(convergence=as.numeric(!modfit$ok))
    modfit$moption <- 2
    modfit$stack <- stack.all
    modfit$runtime <- runtime
    
    modfit
}
