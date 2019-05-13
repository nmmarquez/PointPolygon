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
setwd("~/Documents/PointPolygon/")
sourceCpp("./demo/dist.cpp")
load("./demo/prepData.rda")

nT <- 6

args <- commandArgs(trailingOnly=TRUE)

# for testing
# rangeE <- .3
# covVal <- 2
# covType <- "random"
# seed <- as.integer(123)

rangeE <- as.numeric(args[1]) # range of spatial prces varies from {.2, .4, .6}
covVal <- as.numeric(args[2]) # covariate effect in set {.2, .4, -.5, .2, -2}
covType <- args[3] # either random spatial or cluster 
seed <- as.integer(args[4]) # RNG

modelname <- paste0(
    "range=", rangeE,
    ",cov=", covVal,
    ",covtype=", covType,
    ",seed=", seed, ".Rds"
)

set.seed(seed)

find_closest <- function(x1, x2) {
    
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

maxYear <- max(c(pointDF$year, polyDF$year))

field <- simField(
    N = as.matrix(select(fullDF, long, lat)), rangeE = rangeE,
    offset = c(0.1, 0.2),
    max.edge = c(0.18,0.2),
    beta0 = -2,
    betaList = list(list(type=covType, value=covVal)),
    nTimes = nT,
    rho = .85, shape=spDF)

ggField(field)

fullDF <- fullDF %>%
    rename(x=long, y=lat, oldid=id) %>%
    mutate(strat=paste(sprintf("%02d", reg), (2-urban), sep="_")) %>%
    left_join(
        select(filter(as_tibble(field$spdf), tidx==0), x, y, id),
        by = c("x", "y"))

yearWDF <- yearWDF %>%
    rename(oldid=id) %>%
    mutate(tidx=year-maxYear-1+nT) %>%
    left_join(select(fullDF, oldid, id), by="oldid") %>%
    filter(year <= maxYear) %>%
    select(-year, -oldid) %>%
    arrange(id) %>%
    filter(tidx>=0)

syearWDF <- filter(yearWDF, tidx == max(tidx))

locDF <- syearWDF %>%
    group_by(strat) %>%
    summarize(id=list(id))

# for each psu resample points based off of the population weights of locations
polyDF <- polyDF %>%
    mutate(tidx=year-maxYear-1+nT) %>%
    group_by(strat, tidx) %>%
    summarize(samples=sum(N)) %>%
    ungroup %>%
    filter(tidx >= 0) %>%
    right_join(
        polyDF %>%
            select(psu, strat) %>%
            unique() %>%
            mutate(trueid=sapply(strat, function(s){
                subDF <- filter(syearWDF, strat==s)
                sample(subDF$id, 1, prob=subDF$popW)
            })) %>%
            left_join(tibble(
                psu=rep(.$psu, nT),
                tidx=rep((1:nT)-1, each=nrow(.))
            ), by="psu"),
        by=c("strat", "tidx")) %>%
    rename(id=trueid) %>%
    # We dont need to weight the trials by the weight of the sample of the 
    # pixel as that has already been done in the PSU selection process
    # left_join(select(syearWDF, id, Population), by="id") %>%
    # group_by(strat, tidx) %>%
    # mutate(w=Population/sum(Population)) %>%
    group_by(strat, tidx) %>%
    mutate(w=1/n()) %>%
    ungroup %>%
    mutate(trials=round(w*samples)) %>%
    left_join(
        select(as_tibble(field$spdf), theta, id, tidx),
        by = c("id", "tidx")) %>%
    mutate(obs=rbinom(n(), trials, prob=theta)) %>%
    rename(trueid=id) %>%
    left_join(locDF, by="strat") %>%
    mutate(polyid=as.numeric(as.factor(strat))-1) %>%
    select(id, tidx, trials, obs, polyid, trueid, strat) %>%
    filter(trials != 0)

stratOrder <- arrange(
    as.data.frame(unique(select(polyDF, polyid, strat))), polyid)

pointDF <- pointDF %>%
    rename(oldid=id) %>%
    mutate(tidx=year-maxYear-1+nT) %>%
    filter(tidx>=0) %>%
    left_join(select(fullDF, oldid, id), by="oldid") %>%
    left_join(
        select(as_tibble(field$spdf), theta, id, tidx),
        by = c("id", "tidx")) %>%
    rename(trials=N) %>%
    mutate(obs=rbinom(n(), trials, prob=theta)) %>%
    arrange(tidx, id) %>%
    select(id, tidx, trials, obs)

# Turns out IHME only has one point or each large admin level :/
ihmePolyDF <- read_csv("./demo/ihmeResampleDF.csv") %>%
    select(longitude, latitude) %>%
    unique %>%
    as.matrix %>%
    find_closest(
        syearWDF %>%
            left_join(
                select(as_tibble(field$spdf), -geometry), 
                by = c("tidx", "id")) %>%
            select(x, y) %>%
            as.matrix()
    ) %>%
    tibble(id=.) %>%
    mutate(id=id-1) %>%
    left_join(syearWDF, by="id") %>%
    rename(trueid=id) %>%
    mutate(location_code=as.numeric(str_split(strat, "_", simplify=T)[,1])) %>%
    select(trueid, location_code) %>%
    right_join(
        polyDF_ %>%
            mutate(location_code=
                       as.numeric(str_split(strat, "_", simplify=T)[,1])) %>%
            select(-trueid) , 
        by="location_code") %>%
    group_by(trueid, tidx, location_code) %>%
    summarise(
        trials=sum(trials), obs=sum(obs), id=id[1]) %>%
    ungroup() %>%
    mutate(polyid=as.numeric(as.factor(trueid))-1) %>%
    select(id, tidx, trials, obs, polyid, trueid) %>%
    arrange(tidx, polyid)

AprojPoly <- syearWDF %>%
    left_join(stratOrder, by="strat") %>%
    mutate(col=polyid+1, row=id+1) %>%
    {sparseMatrix(
        i = .$row,
        j = .$col,
        x = .$popW,
        dims = c(max(.$row), max(.$col)))}

### The approaches that we want to test are as follows
# 1 Mixture model
# 2 IHME Resample
# 3 Riemman
# 4 Ignore
# 5 Known

modelList <- list(
    `IHME Resample` = runFieldModel(
        field, pointDF, ihmePolyDF, moption=5, verbose=T
    ),
    `Riemann` = runFieldModel(
        field, pointDF, select(polyDF, -strat), moption=3, verbose=T,
        AprojPoly=AprojPoly
    ),
    `Ignore` = runFieldModel(
        field, pointDF, polyDF, moption=4, verbose=T, AprojPoly=AprojPoly
    ),
    `Known` = runFieldModel(
        field, pointDF, polyDF, moption=5, verbose=T, AprojPoly=AprojPoly
    )
)

modelList[["Mixture Model"]] <- runFieldModel(
    field, pointDF, polyDF, moption=0, verbose=T, AprojPoly=AprojPoly,
    start = list(
        z=matrix(unname(modelList$`IHME Resample`$sd$par.random), ncol=nT),
        beta=modelList$`IHME Resample`$opt$par[names(
            modelList$`IHME Resample`$opt$par) == "beta"])
)

unitFitList <- lapply(modelList, function(model){
    simulateFieldCI(field, model)
})

unitBetaList <- lapply(modelList, function(model){
    simulateBetas(field, model)
})

convergeList <- lapply(modelList, function(model){
    model$opt$convergence
})

timeList <- lapply(modelList, function(model){
    as.numeric(model$runtime, units="mins")
})

aggShapes <- readRDS("demo/aggShapes.Rds")

provPredList <- lapply(modelList, function(model){
    arealCI(
        field, 
        model, 
        polygonList = lapply(1:nrow(aggShapes$provShape@data), function(i){
            sf::st_as_sf(aggShapes$provShape[i,])}),
        popDF = mutate(select(syearWDF, -tidx), w=Population))
})

regPredList <- lapply(modelList, function(model){
    arealCI(
        field, 
        model, 
        polygonList = lapply(1:nrow(aggShapes$regShape@data), function(i){
            sf::st_as_sf(aggShapes$regShape[i,])}),
        popDF = mutate(select(syearWDF, -tidx), w=Population))
})


unitResults <- list(
    sim = field,
    pred = unitFitList,
    betas = unitBetaList,
    model = modelList, # this takes up a lot of space so ignore for now
    converge = convergeList,
    covType = covType,
    rangeE = rangeE,
    covVal = covVal,
    seed = seed,
    runtime = timeList,
    provPred = provPredList,
    regPred = regPredList
)

saveRDS(unitResults, file=paste0("~/Data/spaceTimeTest/", modelname))
