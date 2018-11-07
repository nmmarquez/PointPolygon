rm(list=ls())
library(rgeos)
library(sp)
library(PointPolygon)

args <- commandArgs(trailingOnly=TRUE)

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
unitModelList <- c(
    lapply(polyDFList, function(pdf){
        runFieldModel(
            unitSim,
            polyDF=pdf, 
            verbose=T,
            control=list())}),
    lapply(mixDFList, function(mix){
        runFieldModel(
            unitSim,
            pointDF = mix$pointDF,
            polyDF = mix$polyDF, 
            verbose=T,
            control=list())}),
    lapply(ovDFList, function(ov){
        runFieldModel(
            unitSim,
            pointDF = ov$pointDF,
            polyDF = ov$polyDF, 
            verbose=T,
            control=list())}),
    list(point=runFieldModel(unitSim, pointDF, verbose=T, control=list()))
)

names(unitModelList) <- c(paste0(
    rep(paste0("rwidth_", rWidthSamples), 3),
    rep(c(" poly", " mix", " ov"), each=length(rWidthSamples))), "point")

unitFitList <- lapply(unitModelList, function(ffit){
    simulateFieldCI(unitSim, ffit)})

unitResults <- list(
    sim = unitSim,
    pred = unitFitList,
    model = unitModelList,
    covType = covType,
    rangeE = rangeE,
    covVal = covVal,
    M = M,
    seed = seed
)

saveRDS(unitResults, file=paste0("~/Data/utaziTest/", modelname))