rm(list=ls())
library(sp)
library(PointPolygon)
library(INLA)
library(dplyr)
library(spdep)

unitSim <- simField(
    N = 60, rangeE = .7,
    offset = c(0.1, 0.2), 
    max.edge = c(0.1,0.2),
    beta0 = -2,
    betaList = list(list(type="random", value=2)))

mixSample <- samplePPMix(unitSim, 30, 150, .5, rWidth=3)
pointDF <- mixSample$pointDF
polyDF <- mixSample$polyDF

shape3 <- dividePolygon(unitSim$bound, 3)
shape3@data <- left_join(shape3@data, polyDF)

pointAreas <- pointDF %>%
    left_join(unitSim$spdf@data) %>%
    select(x, y) %>%
    SpatialPoints %>%
    over(shape3) %>%
    select(polyid) %>%
    unlist %>%
    unname

pointA <- unitSim$AprojField[pointDF$id + 1,]

grid.poly.no <- over(unitSim$spdf, shape3)$polyid
ran.pr <- as.numeric(summary(dist(unitSim$spdf@coords))[3])
kap.pr <- sqrt(8)/ran.pr
spde.model <- inla.spde2.matern(
    mesh=unitSim$mesh, alpha=2,
    B.tau=matrix(c(0, 1, 0),nrow=1,ncol=3),
    B.kappa=matrix(c(0, 0, 1),nrow=1,ncol=3),
    theta.prior.mean=c(0, log(kap.pr)),
    theta.prior.prec=c(1, 1))

area.poly.nb <- poly2nb(shape3, snap = 1)
area.poly.adj <- as(nb2mat(area.poly.nb, style = "B"), "dgTMatrix")

# extract covariate values for polygons
covPoly <- unitSim$spdf@data %>%
    select(matches("V[0-9]+")) %>%
    mutate(polyid=grid.poly.no) %>%
    group_by(polyid) %>%
    summarise_all(mean, na.rm=T)

covPoint <- unitSim$spdf@data %>%
    right_join(pointDF) %>%
    select(matches("V[0-9]+")) %>%
    as.data.frame
    
mesh.coord.in <-
    unitSim$mesh$loc[as.vector(which(!is.na(over(SpatialPoints(unitSim$mesh$loc), shape3)$polyid))),]

mesh.coord.poly.no <- over(SpatialPoints(mesh.coord.in), shape3)$polyid

areaA <- inla.spde.make.A(
    mesh=unitSim$mesh,
    loc=mesh.coord.in,
    block=mesh.coord.poly.no, block.rescale="sum")

stack.area <- inla.stack(
    tag='areal',
    data=list(y=shape3$obs, n=shape3$trials),
    A=list(areaA,1,1),
    effects=list(
        s=1:unitSim$spde$n.spde,
        sa=1:nrow(shape3@data),
        as.data.frame(select(covPoly, -polyid))))

stack.pred <- inla.stack(
    tag='pred',
    data=list(
        y=rep(NA, nrow(unitSim$spdf@data)),
        n=rep(NA, nrow(unitSim$spdf@data))),
    A=list(unitSim$AprojField,1,1),
    effects=list(
        s=1:unitSim$spde$n.spde,
        sa=grid.poly.no,
        unitSim$spdf@data %>%
            select(matches("V[0-9]+"))))

stack.points <- inla.stack(
    tag = 'points',
    data=list(y=pointDF$obs, n=pointDF$trials),
    A=list(pointA,1,1),
    effects=list(
        s=1:unitSim$spde$n.spde,
        sa=pointAreas,
        covPoint))

stack.all <- inla.stack(stack.points, stack.area, stack.pred)

f <- y ~ -1 + V0 + V1 + f(s, model=unitSim$spde) + f(
    sa, 
    model = "besagproper2", 
    graph = area.poly.adj,
    hyper=list(prec=list(prior="loggamma", param=c(1,0.01))))

modfit <- inla(
    f, 
    data = inla.stack.data(stack.all),
    family = "binomial",
    Ntrials = stack.all$data$data$n,
    control.predictor=list(compute=TRUE, A=inla.stack.A(stack.all), link=1),
    control.compute=list(dic=TRUE))

summary(modfit)
