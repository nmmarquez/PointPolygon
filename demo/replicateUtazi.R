.libPaths(c("~/R3.5/", .libPaths()))
rm(list=ls())
library(rgeos)
library(sp)
library(PointPolygon)

args <- commandArgs(trailingOnly=TRUE)

# for testing
# rangeE <- .3 
# covVal <- 2 
# covType <- "random"  
# M <- as.integer(150)
# seed <- as.integer(123) 

rangeE <- as.numeric(args[1]) # range of spatial prces varies from {.3, .5, .7}
covVal <- as.numeric(args[2]) # covariate effect in set {.2, .4, -.5, .2, -2}
covType <- args[3] # either random spatial or cluster 
M <- as.integer(args[4]) # number of samples in Polygons chosen from U(50, 300)
seed <- as.integer(args[5]) # RNG

modelname <- paste0(
    "range=", rangeE,
    ",cov=", covVal,
    ",covtype=", covType,
    ",M=", M,
    ",seed=", seed, ".Rds"
)

set.seed(seed)

unitSim <- simField(
    N = 60, # use 60 squares as base field just as in example
    sigmaE = 1, # unit variance for spde process
    rangeE = rangeE,
    shape = NULL, # null shape which by default creates unit square
    beta0 = -2, # intercept
    betaList = list(list(type=covType, value=covVal)), 
    link = arm::invlogit,
    offset = c(0.1, 0.2),
    max.edge = c(0.05,0.2))


rWidthSamples <- c("3"=3, "5"=5, "10"=10)
pointSamples <- 250

# sample from only polygons
polyDFList <- lapply(rWidthSamples, function(x){
    samplePolygons(
        unitSim,
        M = round(M*100/(x^2)),
        p = 1, # sample all polygons
        rWidth = x)})

# sample with points only
pointDF <- samplePoints(unitSim, pointSamples, round(M*100/pointSamples))

# sample with mix of nonoverlapping
mixDFList <- lapply(rWidthSamples, function(x){
    samplePPMix(unitSim, pointSamples, round(M*100/pointSamples), rWidth=x)
})

# sample with mix overlapping units
ovDFList <- lapply(rWidthSamples, function(x){
    list(
        pointDF = samplePoints(
            unitSim, 
            pointSamples,
            round(M*100/pointSamples/2)),
        polyDF = samplePolygons(
            unitSim,
            M = round(M*100/(x^2)/2),
            p = 1, # sample all polygons
            rWidth = x)
    )})

# run all the models!!!
unitModelList <- list(
    riemann = c(
        lapply(1:length(polyDFList), function(i){
            runFieldModel(
                unitSim,
                polyDF=polyDFList[[i]], 
                verbose=T,
                control=list())}),
        lapply(1:length(mixDFList), function(i){
            runFieldModel(
                unitSim,
                pointDF = mixDFList[[i]]$pointDF,
                polyDF = mixDFList[[i]]$polyDF, 
                verbose=T,
                control=list())}),
        lapply(1:length(ovDFList), function(i){
            runFieldModel(
                unitSim,
                pointDF = ovDFList[[i]]$pointDF,
                polyDF = ovDFList[[i]]$polyDF, 
                verbose=T,
                control=list())})
        ),
    resample = c(
        lapply(1:length(polyDFList), function(i){
            runFieldModel(
                unitSim,
                polyDF=polyDFList[[i]], 
                verbose=T,
                control=list(),
                moption=1)}),
        lapply(1:length(mixDFList), function(i){
            runFieldModel(
                unitSim,
                pointDF = mixDFList[[i]]$pointDF,
                polyDF = mixDFList[[i]]$polyDF, 
                verbose=T,
                control=list(),
                moption=1)}),
        lapply(1:length(ovDFList), function(i){
            runFieldModel(
                unitSim,
                pointDF = ovDFList[[i]]$pointDF,
                polyDF = ovDFList[[i]]$polyDF, 
                verbose=T,
                control=list(),
                moption=1)})
        ),
    utazi = c(
        lapply(1:length(polyDFList), function(i){
            runFieldModel(
                unitSim,
                polyDF=polyDFList[[i]], 
                verbose=T,
                control=list(),
                moption=2,
                rWidth=rWidthSamples[i])}),
        lapply(1:length(mixDFList), function(i){
            runFieldModel(
                unitSim,
                pointDF = mixDFList[[i]]$pointDF,
                polyDF = mixDFList[[i]]$polyDF, 
                verbose=T,
                control=list(),
                moption=2,
                rWidth=rWidthSamples[i])}),
        lapply(1:length(ovDFList), function(i){
            runFieldModel(
                unitSim,
                pointDF = ovDFList[[i]]$pointDF,
                polyDF = ovDFList[[i]]$polyDF, 
                verbose=T,
                control=list(),
                moption=2,
                rWidth=rWidthSamples[i])})
        ),
    point=list(point=runFieldModel(unitSim, pointDF, verbose=T, control=list()))
)

names(unitModelList$riemann) <- paste0(
    rep(paste0("rwidth_", rWidthSamples), 3),
    rep(c(" poly", " mix", " ov"), each=length(rWidthSamples)))
names(unitModelList$resample) <- paste0(
    rep(paste0("rwidth_", rWidthSamples), 3),
    rep(c(" poly", " mix", " ov"), each=length(rWidthSamples)))
names(unitModelList$utazi) <- paste0(
    rep(paste0("rwidth_", rWidthSamples), 3),
    rep(c(" poly", " mix", " ov"), each=length(rWidthSamples)))

unitFitList <- lapply(unitModelList, function(type){
    lapply(type, function(ffit){
        simulateFieldCI(unitSim, ffit)})
})

unitBetaList <- lapply(unitModelList, function(type){
    lapply(type, function(ffit){
        simulateBetas(unitSim, ffit)})
})

unitResults <- list(
    sim = unitSim,
    pred = unitFitList,
    betas = unitBetaList,
    # model = unitModelList, # this takes up a lot of space so ignore for now
    covType = covType,
    rangeE = rangeE,
    covVal = covVal,
    M = M,
    seed = seed
)

saveRDS(unitResults, file=paste0("~/Data/utaziTest/", modelname))
