# qsub -e ~/errors/ -o ~/outputs/ -l mem_free=999G -l m_mem_free=999G \
# -P proj_geo_nodes_u5m -l fthread=10 -q geospatial.q -N datarun \
# -l h_rt='07:00:00:00' \
# /share/singularity-images/lbd/shells/singR.sh -m 10 -o 5 \
# -e s ~/Documents/PointPolygon/demo/dataRun.R
.libPaths(c("~/R3.6/", .libPaths()))
rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(PointPolygon)
library(sp)
library(Rcpp)
library(Matrix)
library(TMB)
library(sparseMVN)
setwd("~/Documents/PointPolygon/")
sourceCpp("./demo/dist.cpp")
load("./demo/prepData.rda")
setwd("./src/")

args <- commandArgs(trailingOnly=TRUE)
yearHO <- as.numeric(args[1]) 
regHO <- as.numeric(args[2]) 

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

maxYear <- 2014
minYear <- 2000
nYears <- length(minYear:maxYear)
ageVec <- c(NN=0, PNN1=1, PNN2=2, `1yr`=3, `2yr`=4, `3yr`=5, `4yr`=6)

field <- simField(
    N = as.matrix(select(fullDF, long, lat)), rangeE = 1,
    offset = c(0.1, 0.2),
    max.edge = c(0.18,0.2),
    beta0 = -4,
    sigmaE = .3,
    nTimes = nYears,
    rho = .9, shape=spDF)

fullDF <- fullDF %>%
    rename(x=long, y=lat, oldid=id) %>%
    mutate(strat=paste(sprintf("%02d", reg), (2-urban), sep="_")) %>%
    left_join(
        select(filter(as_tibble(field$spdf), tidx==0), x, y, id),
        by = c("x", "y"))

syearWDF <- yearWDF %>%
    rename(oldid=id) %>%
    mutate(tidx=year-maxYear-1+nYears) %>%
    left_join(select(fullDF, oldid, id), by="oldid") %>%
    filter(year <= maxYear) %>%
    select(-year, -oldid) %>%
    arrange(id) %>%
    filter(tidx==max(tidx))

covArray <- array(1, dim=c(nrow(fullDF), nYears, 2))
covArray[,,2] <- sapply(1:nYears, function(i) fullDF$urban)

pointDF <- pointDF %>%
    filter((year >= minYear) & (year <= maxYear)) %>%
    mutate(denom=N, obs=died, yid=year-minYear, aid=ageVec[age_group]) %>%
    mutate(sid=as.numeric(as.factor(source))-1) %>%
    mutate(cid=group_indices(., sid, psu)-1) %>%
    rename(oldid=id) %>%
    left_join(select(fullDF, id, oldid), by="oldid") %>%
    select(denom, obs, id, yid, aid, sid, cid)

polyDF <- polyDF %>%
    filter((year >= minYear) & (year <= maxYear)) %>%
    mutate(denom=N, obs=died, yid=year-minYear, aid=ageVec[age_group]) %>%
    mutate(sid=max(pointDF$sid) + 1) %>%
    mutate(cid=group_indices(., psu) + max(pointDF$cid)) %>%
    mutate(polyid=as.numeric(as.factor(strat))-1) %>%
    select(denom, obs, id, yid, aid, sid, cid, polyid, strat, weight)

# Turns out IHME only has one point or each large admin level :/
ihmePolyDF <- read_csv("../demo/ihmeResampleDF.csv") %>%
    select(longitude, latitude) %>%
    unique %>%
    as.matrix %>%
    find_closest(
        syearWDF %>%
            left_join(
                select(as_tibble(field$spdf), -geometry), 
                by = c("tidx", "id")) %>%
            select(x, y) %>%
            as.matrix(),
        fullDF
    ) %>%
    tibble(id=.) %>%
    mutate(id=id-1) %>%
    left_join(syearWDF, by="id") %>%
    rename(trueid=id) %>%
    mutate(location_code=as.numeric(str_split(strat, "_", simplify=T)[,1])) %>%
    select(trueid, location_code) %>%
    right_join(
        polyDF %>%
            mutate(location_code=
                       as.numeric(str_split(strat, "_", simplify=T)[,1])), 
        by="location_code") %>%
    select(-id, -polyid, -location_code) %>%
    rename(id=trueid)

utaziPolyDF <- polyDF %>%
    group_by(yid, aid, sid, strat) %>%
    summarize(
        obs = sum(weight * obs),
        denom = sum(weight * denom)) %>%
    mutate(urban=(str_split(strat, "_", simplify = T)[,2] == "1")+0) %>%
    mutate(polyid=(as.numeric(str_split(strat, "_", simplify = T)[,1]))-1) %>%
    ungroup()

shape3 <- rgeos::gUnaryUnion(spDF, id = spDF$REG)
shape3 <- SpatialPolygonsDataFrame(
    shape3, data.frame(reg=as.numeric(names(shape3)), row.names=names(shape3)))
shape3 <- spTransform(shape3, "+proj=longlat +ellps=WGS84 +no_defs")
shape3$polyid <- shape3$reg - 1

stratOrder <- arrange(
    as.data.frame(unique(select(polyDF, polyid, strat))), polyid)
AprojPoly <- syearWDF %>%
    left_join(stratOrder, by="strat") %>%
    mutate(col=polyid+1, row=id+1) %>%
    {sparseMatrix(
        i = .$row,
        j = .$col,
        x = .$popW,
        dims = c(max(.$row), max(.$col)))}

field$spdf <- field$spdf %>%
    left_join(select(fullDF, id, urban), by="id")

if(!is.na(yearHO)){
    pointDF <- filter(pointDF, yid != yearHO)
    polyDF <- filter(polyDF, yid != yearHO)
    ihmePolyDF <- filter(ihmePolyDF, yid != yearHO)
}
    
if(!is.na(regHO)){
    pointDF <- pointDF %>%
        left_join(select(fullDF, id, strat)) %>%
        mutate(location_code=
                   as.numeric(str_split(strat, "_", simplify=T)[,1]) - 1) %>%
        filter(location_code != regHO) 
    polyDF <- polyDF %>%
        mutate(location_code=
                   as.numeric(str_split(strat, "_", simplify=T)[,1]) - 1) %>%
        filter(location_code != regHO)
    ihmePolyDF <- ihmePolyDF %>%
        mutate(location_code=
                   as.numeric(str_split(strat, "_", simplify=T)[,1]) - 1) %>%
        filter(location_code != regHO)
}

modelRun <- function(
    pointDF, polyDF=NULL, moption = 0, priors = 0, nugget = TRUE, model = "u5m",
    time_structured = TRUE, time_unstructured = TRUE, survey_effect=FALSE,
    verbose=TRUE, symbolic=TRUE, control=list(eval.max=1e4, iter.max=1e4)){
    if(is.null(polyDF)){
        moption <- 0
        empty <- vector("integer")
        polyDF <- data.frame(
            denom=empty, obs=empty, id=empty, yid=empty, aid=empty,
            sid=empty, cid=empty, polyid=empty, strat=empty)
    }

    Data <- list(
        yPoint=pointDF$obs,
        yPoly=polyDF$obs,
        denomPoint=pointDF$denom,
        denomPoly=polyDF$denom,
        idPoint=pointDF$id,
        idPoly=polyDF$polyid,
        idtPoint=pointDF$yid,
        idtPoly=polyDF$yid,
        idaPoint=pointDF$aid,
        idaPoly=polyDF$aid,
        idnPoint=pointDF$sid,
        idnPoly=polyDF$sid,
        idcPoint=pointDF$cid,
        idcPoly=polyDF$cid,
        covs=covArray,
        AprojObs=field$AprojField,
        AprojPoly=AprojPoly,
        M0=field$spde$param.inla$M0,
        M1=field$spde$param.inla$M1,
        M2=field$spde$param.inla$M2,
        moption=moption,
        priors=priors
    )

    Params <- list(
        beta = c(0, 0),
        beta_age = rep(0, length(ageVec) - 1),
        log_tau = 0,
        log_kappa = 0,
        logit_rho = 0,
        log_sigma_phi = rep(0, length(ageVec)),
        log_sigma_nu = 0,
        log_sigma_epsilon = 0,
        log_sigma_eta = 0,
        z = array(0, dim=c(field$mesh$n, field$nTimes)),
        epsilon = rep(0, field$nTimes),
        phi = array(0,dim=c(length(ageVec), field$nTimes)),
        nu = rep(0, max(c(pointDF$sid, polyDF$sid)) + 1),
        eta = rep(0, max(c(pointDF$cid, polyDF$cid)) + 1)
    )

    Map <- NULL
    random <- c("z", "epsilon", "phi", "nu", "eta")

    if(!nugget){
        Map <- c(
            Map,
            list(
                log_sigma_eta = factor(NA),
                eta = factor(rep(NA, length(Params$eta)))
            )
        )
        random <- random[random != "eta"]
    }

    if((length(Params$nu) < 3) | !survey_effect){
        Map <- c(
            Map,
            list(
                log_sigma_nu = factor(NA),
                nu = factor(rep(NA, length(Params$nu)))
            )
        )
        random <- random[random != "nu"]
    }
    
    if(!time_structured){
        oldDim <- dim(Map$phi)
        Map <- c(
            Map,
            list(
                phi = factor(rep(NA, length(c(Params$phi)))),
                log_sigma_phi = factor(rep(NA, length(Params$log_sigma_phi)))
            )
        )
        dim(Map$phi) <- oldDim
        random <- random[random != "phi"]
    }
    
    if(!time_unstructured){
        Map <- c(
            Map,
            list(
                log_sigma_epsilon = factor(NA),
                epsilon = factor(rep(NA, length(Params$epsilon)))
            )
        )
        random <- random[random != "epsilon"]
    }

    print(str(Map))
    model <- "u5m"
    compile(paste0(model, ".cpp"))
    dyn.load(dynlib(model))
    config(tape.parallel=0, DLL=model)
    
    Obj <- MakeADFun(
        data=Data, parameters=Params, DLL=model, random=random, map = Map,
        silent=!verbose)
    Obj$env$tracemgc <- verbose
    Obj$env$inner.control$trace <- verbose

    startTime <- Sys.time()

    if(symbolic){
        TMB::runSymbolicAnalysis(Obj)
    }
    Opt <- stats::nlminb(
        start = Obj$par,
        objective = Obj$fn,
        gradient = Obj$gr,
        control = control)
    sdrep <- TMB::sdreport(Obj, getJointPrecision=TRUE)
    runtime <- Sys.time() - startTime

    dyn.unload(dynlib(model))

    list(
        obj = Obj,
        opt = Opt,
        runtime = runtime,
        sd = sdrep,
        moption = moption,
        stack=NULL)
}

test <- runDataUtazi(
    field, pointDF = pointDF, polyDF_ = utaziPolyDF, shape3 = shape3,
    fullDF = fullDF, timeStructured = FALSE)
# 
# modelList <- list()
# modelList$noRE <- list()
# modelList$temporal <- list()
# modelList$full <- list()
# 
# modelList$noRE$point <- modelRun(
#     pointDF, polyDF = NULL, nugget = FALSE, time_structured = FALSE, 
#     time_unstructured = FALSE, survey_effect = FALSE, priors = 1)
# 
# modelList$noRE$mixture <- modelRun(
#     pointDF, polyDF = polyDF, nugget = FALSE, time_structured=FALSE, 
#     time_unstructured = FALSE, survey_effect = TRUE, priors = 1)
# 
# modelList$noRE$resample <- modelRun(
#     bind_rows(pointDF, ihmePolyDF), nugget = FALSE, time_structured=FALSE, 
#     time_unstructured = FALSE, survey_effect = TRUE, priors = 1)
# 
# modelList$temporal$point <- modelRun(
#     pointDF, polyDF = NULL, nugget = FALSE, time_structured = TRUE, 
#     time_unstructured = FALSE, survey_effect = FALSE, priors = 1)
# 
# modelList$temporal$mixture <- modelRun(
#     pointDF, polyDF = polyDF, nugget = FALSE, time_structured=TRUE, 
#     time_unstructured = FALSE, survey_effect = TRUE, priors = 1)
# 
# modelList$temporal$resample <- modelRun(
#     bind_rows(pointDF, ihmePolyDF), nugget = FALSE, time_structured=TRUE, 
#     time_unstructured = FALSE, survey_effect = TRUE, priors = 1)
# 
# modelList$full$point <- modelRun(
#     pointDF, polyDF = NULL, nugget = TRUE, time_structured = TRUE, 
#     time_unstructured = FALSE, survey_effect = FALSE, priors = 1)
# 
# modelList$full$mixture <- modelRun(
#     pointDF, polyDF = polyDF, nugget = TRUE, time_structured=TRUE, 
#     time_unstructured = FALSE, survey_effect = TRUE, priors = 1)
# 
# modelList$full$resample <- modelRun(
#     bind_rows(pointDF, ihmePolyDF), nugget = TRUE, time_structured=TRUE, 
#     time_unstructured = FALSE, survey_effect = TRUE, priors = 1)
# 
# saveRDS(
#     list(
#         pointDF=pointDF, polyDF=polyDF, ihmePolyDF=ihmePolyDF,
#         field=field, modelList=modelList),
#     sprintf(
#         "~/Documents/PointPolygon/demo/data_run/re_nonug_model_y%s_r%s_.RDS",
#         yearHO,
#         regHO))
# 
# library(ggplot2)
# predDF <- simulateDataFieldCI(modelList$mixture, field)
# 
# predDF %>%
#     {left_join(field$spdf, .)} %>%
#     ggplot(aes(x, y, fill = mu)) +
#     geom_raster() +
#     coord_equal() +
#     theme_void() +
#     scale_fill_distiller(palette = "Spectral") +
#     facet_grid(aid~tidx)
# 
# predDF %>%
#     group_by(id, tidx) %>%
#     # This is not correct look into this later
#     summarize(
#         mu=1-prod(1-mu),
#         lwr=1-prod(1-lwr),
#         uprr=1-prod(1-upr)) %>%
#     mutate(Year=tidx+2000) %>%
#     {left_join(field$spdf, .)} %>%
#     ggplot(aes(x, y, fill = mu)) +
#     geom_raster() +
#     coord_equal() +
#     theme_void() +
#     scale_fill_distiller(palette = "Spectral") +
#     facet_wrap(~Year)
# 
# yearWDF %>%
#     rename(oldid=id) %>%
#     mutate(tidx=year-maxYear-1+nYears) %>%
#     left_join(select(fullDF, oldid, id), by="oldid") %>%
#     filter(year <= maxYear) %>%
#     select(-year, -oldid) %>%
#     arrange(id) %>%
#     filter(tidx>=0) %>%
#     right_join(predDF) %>%
#     group_by(tidx, aid) %>%
#     summarize(
#         mu=weighted.mean(mu, Population),
#         lwr=weighted.mean(lwr, Population),
#         upr=weighted.mean(upr, Population)) %>%
#     mutate(aid=as.factor(aid)) %>%
#     ggplot(aes(x=tidx, y=mu, ymin=lwr, ymax=upr, group=aid)) +
#     geom_line(aes(color=aid)) +
#     geom_ribbon(aes(fill=aid), alpha=.3) +
#     coord_trans(y="log") +
#     theme_classic()
# 
# yearWDF %>%
#     rename(oldid=id) %>%
#     mutate(tidx=year-maxYear-1+nYears) %>%
#     left_join(select(fullDF, oldid, id), by="oldid") %>%
#     filter(year <= maxYear) %>%
#     select(-year, -oldid) %>%
#     arrange(id) %>%
#     filter(tidx>=0) %>%
#     group_by(tidx, strat) %>%
#     summarize(Population=sum(Population)) %>%
#     mutate(popW=Population/sum(Population)) %>%
#     ungroup %>%
#     right_join(
#         polyDF %>%
#             select(-id) %>%
#             rename(id=polyid, tidx=yid) %>%
#             group_by(tidx, strat, aid) %>%
#             summarize(mu=sum(obs)/sum(denom))) %>%
#     group_by(tidx, aid) %>%
#     summarize(mu=weighted.mean(mu, Population)) %>%
#     mutate(aid=as.factor(aid)) %>%
#     ggplot(aes(x=tidx, y=log(mu), color=aid, group=aid)) +
#     geom_line() +
#     theme_classic()
