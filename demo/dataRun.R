.libPaths(c("~/R3.5/", .libPaths()))
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
setwd("~/Documents/PointPolygon/")
sourceCpp("./demo/dist.cpp")
load("./demo/prepData.rda")
setwd("./src/")

maxYear <- 2015
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
    mutate(denom=N, obs=died, yid=year-minYear, aid=ageVec[age_group]) %>%
    mutate(sid=as.numeric(as.factor(source))-1) %>%
    mutate(cid=group_indices(., sid, psu)-1) %>%
    rename(oldid=id) %>%
    left_join(select(fullDF, id, oldid), by="oldid") %>%
    select(denom, obs, id, yid, aid, sid, cid)
    
polyDF <- polyDF %>%
    mutate(denom=N, obs=died, yid=year-minYear, aid=ageVec[age_group]) %>%
    mutate(sid=2, cid=group_indices(., psu) + max(pointDF$cid) + 1) %>%
    mutate(polyid=as.numeric(as.factor(strat))-1) %>%
    select(denom, obs, id, yid, aid, sid, cid, polyid, strat) 

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

modelRun <- function(
    pointDF, polyDF=NULL, moption=0, priors=0, nugget=TRUE, 
    verbose=TRUE, symbolic=TRUE, control=list(eval.max=1e4, iter.max=1e4)){
    snu <- c(0, 0, 0)
    if(is.null(polyDF)){
        moption <- 0
        snu <- c(0, 0)
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
        nu = snu,
        eta = rep(0, length(unique(pointDF$cid)) + length(unique(polyDF$cid)))
    )
    
    Map <- list(
        log_sigma_eta = factor(NA),
        eta = factor(rep(NA, length(Params$eta)))
    )
    
    random <- c("z", "epsilon", "phi", "nu")
    
    if(nugget){
        Map <- NULL
        random <- c("z", "epsilon", "phi", "nu", "eta")
    }
    
    model <- "u5m"
    compile(paste0(model, ".cpp"))
    dyn.load(dynlib(model))
    config(tape.parallel=0, DLL=model)
    
    startTime <- Sys.time()
    Obj <- MakeADFun(
        data=Data, parameters=Params, DLL=model, random=random, map = Map,
        silent=!verbose)
    Obj$env$tracemgc <- verbose
    Obj$env$inner.control$trace <- verbose
    
    if(symbolic){
        nah <- utils::capture.output(TMB::runSymbolicAnalysis(Obj))
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

test <- modelRun(pointDF, polyDF=NULL, moption=0, nugget = FALSE)

