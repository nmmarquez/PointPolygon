#' Run a Model and Estimate Underlying Field From Point and Polygon Data
#'
#' @description Run a model and estimate underlying latent and probability
#' field from point and polygon data using INLA and Utazi method.
#' 
#' @param field field object which simulated underlying data
#' @param pointDF data simulated from samplePoints
#' @param polyDF data simulated from samplePolygns
#' @param moption int, intger indicating how polygon data should be estimated 0
# is by Reimann sum approximation, 1 is by redistribution, and 2 is by Utazi
#' approach.
#' @param verbose logical, print model fitting information
#' @param symbolic logical, use metas reordering in model fitting
#' @param control list, control list passed to nlminb
#' @param rWidth int, only used in moption 2 to build besag prior
#' @param aggData bool, whether to aggregate data as passed in
#' @param shape3 the shape to build the adjaceny matrix with
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
    rWidth = NULL,
    aggData = TRUE,
    shape3 = NULL
    ){
    library(sp)
    library(dplyr)
    library(INLA)

    if(is.null(pointDF)){
        empty <- vector("integer")
        pointDF <- data.frame(obs=empty, trials=empty, id=empty)
    }
    
    if(aggData){
        # presence of strat indicates this is DR sim model
        if(("strat" %in% names(polyDF))){
            polyDF_ <- polyDF %>% 
                select(polyid, tidx, obs, trials, strat) %>%
                group_by(polyid, tidx, strat) %>% 
                summarize_all(sum) %>%
                as.data.frame %>%
                mutate(urban=(str_split(strat, "_", n=2, simplify=T)[,2])) %>%
                mutate(urban=as.numeric(as.numeric(urban)==1))
        }
        else{
            polyDF_ <- polyDF %>% 
                select(polyid, tidx, obs, trials) %>%
                group_by(polyid, tidx) %>% 
                summarize_all(sum) %>%
                mutate(urban = TRUE) %>%
                as.data.frame
            
            field$spdf$urban <- TRUE 
        }
    }
    
    else{
        polyDF_ <- polyDF
    }

    if(!is.null(rWidth) & is.null(shape3)){
        shape3 <- dividePolygon(field$bound, rWidth)
    }
    
    if(!is.null(rWidth)){
        # this is a direct copy from the Utazi paper
        area.poly.nb <- spdep::poly2nb(shape3, snap = 1)
        area.poly.adj <- as(spdep::nb2mat(area.poly.nb, style = "B"), "dgTMatrix")
    }
    
    if(is.null(rWidth)){
        # this modifys the adjacency structure
        area.poly.nb <- spdep::poly2nb(shape3)
        area.poly.adj <- as(spdep::nb2mat(area.poly.nb, style = "B"), "dgTMatrix")
    }
    
    pointA <- field$AprojField[pointDF$id + 1,]
    
    tmpSPDF <- as(field$spdf, "Spatial")
    tmpSPDF@proj4string <- CRS(st_crs(field$spdf)$proj4string)
    
    grid.poly.no <- sp::over(tmpSPDF, shape3)$polyid
    ran.pr <- as.numeric(summary(dist(sf::st_coordinates(field$spdf)))[3])
    kap.pr <- sqrt(8)/ran.pr
    spde.model <- INLA::inla.spde2.matern(
        mesh=field$mesh, alpha=2,
        B.tau=matrix(c(0, 1, 0),nrow=1,ncol=3),
        B.kappa=matrix(c(0, 0, 1),nrow=1,ncol=3),
        theta.prior.mean=c(0, log(kap.pr)),
        theta.prior.prec=c(1, 1))
    
    # extract covariate values for polygons
    covPoly <- field$spdf %>%
        dplyr::as_tibble() %>%
        select(matches("V[0-9]+"), tidx, urban) %>%
        dplyr::mutate(polyid=grid.poly.no) %>%
        group_by(polyid, tidx, urban) %>%
        summarise_all(mean, na.rm=T)
    
    covPoint <- field$spdf %>%
        dplyr::as_tibble() %>%
        dplyr::select(-geometry) %>%
        dplyr::right_join(pointDF) %>%
        dplyr::select(dplyr::matches("V[0-9]+"), urban) %>%
        as.data.frame

    mesh.coord.in <-
        field$mesh$loc[as.vector(which(!is.na(over(
            SpatialPoints(field$mesh$loc, shape3@proj4string), shape3)$polyid))),]

    mesh.coord.poly.no <- over(SpatialPoints(
        mesh.coord.in, shape3@proj4string), shape3)$polyid

    if(field$nTimes == 1){

        areaA <- inla.spde.make.A(
            mesh=field$mesh,
            loc=mesh.coord.in,
            block=mesh.coord.poly.no+1, block.rescale="sum")

        stack.area <- inla.stack(
            tag='areal',
            data=list(y=polyDF_$obs, n=polyDF_$trials),
            # A maps how you get from the effects estimated to Areal data observation
            # its a 1 for both non SPDE random effects and covariates.
            A=list(areaA[polyDF_$polyid+1,],1,1),
            effects=list(
                s=1:ncol(areaA),
                sa=polyDF_$polyid+1,
                polyDF_ %>%
                    left_join(covPoly, by = c("polyid", "tidx", "urban")) %>%
                    select(V0, V1, urban) %>%
                    as.data.frame()))

        stack.pred <- inla.stack(
            tag='pred',
            data=list(
                y=rep(NA, nrow(field$spdf)),
                n=rep(NA, nrow(field$spdf))),
            A=list(
                do.call(rbind, lapply(1:field$nTimes, function(x) field$AprojField)),
                1,
                1),
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
    }
    
    if(field$nTimes > 1){
        knots <- 1:field$nTimes
        mesh_time <- inla.mesh.1d(loc = knots, degree = 2, boundary = "free")
        index <- inla.spde.make.index(
            "s", n.spde = spde.model$n.spde, n.group = mesh_time$m)

        # create a  mesh data frame that tells you which mesh points 
        # correspond to which polygons
        meshDF <- tibble::tibble(
            mesh_coord_poly_no = mesh.coord.poly.no,
            og_position = 1:length(mesh.coord.poly.no)) %>%
            arrange(mesh_coord_poly_no)

        # take the original matrix that was created and repeat it for the
        # number of obseravtions that we have
        mesh.coord.in.new <- do.call(rbind, lapply(polyDF_$polyid, function(i){
            mesh.coord.in[meshDF$og_position[meshDF$mesh_coord_poly_no == i],]
        }))

        # this counts for which obs each row in the above matrix corresponds to
        # and 
        mesh.coord.poly.no.new <- unlist(lapply(1:nrow(polyDF_), function(i){
            rep(i, sum(meshDF$mesh_coord_poly_no == polyDF_$polyid[i]))
        }))

        # this tells us the time each 
        mesh.coord.time <- unlist(lapply(1:nrow(polyDF_), function(i){
            polyDF_ %>%
                filter()
            repCount <- sum(meshDF$mesh_coord_poly_no == polyDF_$polyid[i])
            rep(polyDF_$tidx[i], repCount)
        }))

        # this first projection Matrix maps from the mesh nodes to a SxT 
        # of all polygons observed in the data set
        areaA <- inla.spde.make.A(
            mesh=field$mesh,
            loc=mesh.coord.in.new,
            # blocks need to be independently defined per time periods
            block=mesh.coord.poly.no.new,
            block.rescale="sum",
            group = mesh.coord.time + 1,
            group.mesh = mesh_time
        )

        stack.area <- inla.stack(
            tag='areal',
            data=list(y=polyDF_$obs, n=polyDF_$trials),
            # A maps how you get from the effects estimated to Areal data observation
            # its a 1 for both non SPDE random effects and covariates.
            A=list(areaA,1,1),
            effects=list(
                s=index,
                sa=polyDF_$polyid+1,
                polyDF_ %>%
                    left_join(covPoly, by = c("polyid", "tidx", "urban")) %>%
                    select(V0, V1, urban) %>%
                    as.data.frame()))
        
        predSpaceA <- inla.spde.make.A(
            mesh=field$mesh,
            loc=field$spdf %>%
                as.data.frame() %>%
                select(x, y) %>%
                as.matrix(),
            group = field$spdf$tidx + 1,
            group.mesh = mesh_time
        )

        stack.pred <- inla.stack(
            tag='pred',
            data=list(
                y=rep(NA, nrow(field$spdf)),
                n=rep(NA, nrow(field$spdf))),
            A=list(predSpaceA,1,1),
            effects=list(
                s=index,
                sa=grid.poly.no+1,
                field$spdf %>%
                    dplyr::as_tibble() %>%
                    dplyr::select(dplyr::matches("V[0-9]+"), urban) %>%
                    as.data.frame))

        if(nrow(pointDF) != 0){
            pointAreas <- pointDF %>%
                dplyr::left_join(dplyr::as_tibble(field$spdf)) %>%
                select(x, y) %>%
                SpatialPoints(proj4string = shape3@proj4string) %>%
                over(shape3) %>%
                select(polyid) %>%
                unlist %>%
                unname
            
            pointSpaceA <- inla.spde.make.A(
                mesh=field$mesh,
                loc=pointDF %>% 
                    left_join(select(as.data.frame(field$spdf), x, y, id, tidx)) %>%
                    select(x, y) %>%
                    as.matrix(),
                group = pointDF$tidx + 1,
                group.mesh = mesh_time
            )
            
            stack.points <- inla.stack(
                tag = 'points',
                data=list(y=pointDF$obs, n=pointDF$trials),
                A=list(pointSpaceA,1,1),
                effects=list(
                    s=index,
                    sa=pointAreas+1,
                    covPoint))
            
            stack.all <- inla.stack(stack.points, stack.area, stack.pred)
        }
        
        else{
            stack.all <- inla.stack(stack.area, stack.pred)
        }
        f <- y ~ -1 + V0 + V1 + urban + 
            f(
                s,
                model=spde.model,
                group=s.group,
                control.group = list(
                    model = "ar1",
                    hyper = list(
                        theta1 = list(prior="pc.cor0", param=c(0.9, 0.9))))) + 
            f(
                sa,
                model = "besagproper2", 
                graph = area.poly.adj,
                hyper=list(prec=list(prior="loggamma", param=c(1,0.01))))
    }
    
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
