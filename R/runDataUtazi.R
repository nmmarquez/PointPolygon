#' Run a Model and Estimate Underlying Field From Point and Polygon Data for the
#' Utazi model.
#'
#' @description Run the utazi model on the actual DR data
#'
#' @param field field object which simulated underlying data
#' @param pointDF data simulated from samplePoints
#' @param polyDF_ data simulated from samplePolygns
#' @param shape3 shape to build adjaceny matrix
#' @param fullDF id to xy location and urbanicity status
#' @param timeStructured bool add time structured effects
#' @return Model fit object
#' @export

runDataUtazi <- function(
    field, pointDF, polyDF_, shape3, fullDF, timeStructured = FALSE){
  area.poly.nb <- spdep::poly2nb(shape3)
  area.poly.adj <- as(spdep::nb2mat(area.poly.nb, style = "B"), "dgTMatrix")
  
  tmpSPDF <- as(field$spdf, "Spatial")
  tmpSPDF@proj4string <- CRS(st_crs(field$spdf)$proj4string)
  
  dummify <- function(y){
      data.frame(sapply(1:7, function(x) as.numeric(x==y)))
  }
  
  timeify <- function(y){
      df_ <- data.frame(sapply(1:7, function(x) y))
      names(df_) <- paste0("T", 1:7)
      df_
  }
  
  grid.poly.no <- sp::over(tmpSPDF, shape3)$polyid
  ran.pr <- 1.152895 # as.numeric(summary(dist(sf::st_coordinates(field$spdf)))[3])
  kap.pr <- sqrt(8)/ran.pr
  spde.model <- INLA::inla.spde2.matern(
    mesh=field$mesh, alpha=2,
    B.tau=matrix(c(0, 1, 0),nrow=1,ncol=3),
    B.kappa=matrix(c(0, 0, 1),nrow=1,ncol=3),
    theta.prior.mean=c(0, log(kap.pr)),
    theta.prior.prec=c(1, 1))
  
  predDF <- bind_rows(lapply(0:6, function(a){
    field$spdf %>%
      as_tibble() %>%
      select(id, tidx, urban, x, y) %>%
      mutate(aid=a)})) %>%
    mutate(mu=NA, sd=NA, lwr=NA, upr=NA, sid=4)
  
  area.poly.nb <- spdep::poly2nb(shape3)
  area.poly.adj <- as(spdep::nb2mat(area.poly.nb, style = "B"), "dgTMatrix")
  
  pointA <- field$AprojField[pointDF$id + 1,]
  
  covPoly <- mutate(select(polyDF_, urban, yid, aid), INT=1) %>%
      mutate(yid=yid+1, sid=3, aid=aid+1)

  covPoint <- pointDF %>%
      select(id, yid, sid, aid) %>%
      left_join(fullDF, by="id") %>%
      select(urban, yid, sid, aid) %>%
      mutate(yid=yid+1, INT=1, sid=sid+1, aid=aid+1)

  mesh.coord.in <-
      field$mesh$loc[as.vector(which(!is.na(over(
          SpatialPoints(field$mesh$loc, shape3@proj4string), shape3)$polyid))),]
  
  mesh.coord.poly.no <- over(SpatialPoints(
      mesh.coord.in, shape3@proj4string), shape3)$polyid
  
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
      repCount <- sum(meshDF$mesh_coord_poly_no == polyDF_$polyid[i])
      rep(polyDF_$yid[i], repCount)
  }))
  
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
      data=list(y=polyDF_$obs, n=polyDF_$denom),
      # A maps how you get from the effects estimated to Areal data observation
      # its a 1 for both non SPDE random effects and covariates.
      A=list(areaA,1,1),
      effects=list(
          s=index,
          sa=polyDF_$polyid+1,
          covPoly %>%
              cbind(dummify(.$aid)) %>%
              cbind(timeify(.$yid))))
  
  predSpaceA <- inla.spde.make.A(
      mesh=field$mesh,
      loc=predDF %>%
          select(x, y) %>%
          as.matrix(),
      group = predDF$tidx + 1,
      group.mesh = mesh_time
  )
  
  grid.poly.no <- cbind(predDF$x, predDF$y) %>%
      sp::SpatialPoints(proj4string = shape3@proj4string) %>%
      sp::over(shape3) %>%
      .$polyid
  
  stack.pred <- inla.stack(
      tag='pred',
      data=list(
          y=rep(NA, nrow(predDF)),
          n=rep(NA, nrow(predDF))),
      A=list(predSpaceA,1,1),
      effects=list(
          s=index,
          sa=grid.poly.no+1,
          predDF %>%
              mutate(yid=tidx+1, sid=4, INT=1, aid=aid+1) %>%
              select(urban, yid, INT, sid, aid) %>%
              cbind(dummify(.$aid)) %>%
              cbind(timeify(.$yid))))
  
  pointAreas <- pointDF %>%
      dplyr::left_join(fullDF) %>%
      select(x, y) %>%
      SpatialPoints(proj4string = shape3@proj4string) %>%
      over(shape3) %>%
      select(polyid) %>%
      unlist %>%
      unname
  
  pointSpaceA <- inla.spde.make.A(
      mesh=field$mesh,
      loc=pointDF %>% 
          left_join(select(fullDF, id, x, y), by="id") %>%
          select(x, y) %>%
          as.matrix(),
      group = pointDF$yid + 1,
      group.mesh = mesh_time
  )
  
  stack.points <- inla.stack(
      tag = 'points',
      data=list(y=pointDF$obs, n=pointDF$denom),
      A=list(pointSpaceA,1,1),
      effects=list(
          s=index,
          sa=pointAreas+1,
          covPoint %>%
              cbind(dummify(.$aid)) %>%
              cbind(timeify(.$yid))))
  
  stack.all <- inla.stack(stack.points, stack.area, stack.pred)
  
  f_ <- y ~ -1 + INT + urban + as.factor(aid) +
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
          hyper=list(prec=list(prior="loggamma", param=c(1,0.01)))) +
      f(sid, model="iid")
  
  if(timeStructured){
      f_ <- y ~ -1 + INT + urban + as.factor(aid) +
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
              hyper=list(prec=list(prior="loggamma", param=c(1,0.01)))) +
          f(sid, model="iid") +
          f(T1, X1, model="rw2") +
          f(T2, X2, model="rw2") +
          f(T3, X3, model="rw2") +
          f(T4, X4, model="rw2") +
          f(T5, X5, model="rw2") +
          f(T6, X6, model="rw2") +
          f(T7, X7, model="rw2")
  }
  
  startTime <- Sys.time()
  
  modfit <- inla(
      f_, 
      data = inla.stack.data(stack.all),
      family = "binomial",
      Ntrials = stack.all$data$data$n,
      control.predictor=list(compute=TRUE, A=inla.stack.A(stack.all), link=1),
      verbose = TRUE)
  
  runtime <- Sys.time() - startTime
  
  modfit$opt <- list(convergence=as.numeric(!modfit$ok))
  modfit$moption <- 2
  modfit$stack <- stack.all
  modfit$runtime <- runtime

  modfit
}