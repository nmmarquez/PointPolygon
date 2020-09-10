.libPaths(c("~/R3.5/", .libPaths()))
#rm(list=ls())
library(rgeos)
library(sp)
library(PointPolygon)
library(dplyr)

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

# note that we may potentially want to change these tests so that it reflects a
# pixel count and mesh node number similar to IHME u5m tests. For example in
# brazil we have the following dimensions
# full proj raster: 1e+06
# mesh nodes: 988

rWidthSamples <- c("3"=3, "5"=5, "10"=10)
pointSamples <- 250
totalTrials <- pointSamples*M

# sample from only polygons
polyDFList <- lapply(rWidthSamples, function(x){
    samplePerPolygon <- totalTrials / x^2
    sampleTrailClusterRatio <- pointSamples / M
    clustersPerPolygon <- round(sqrt(samplePerPolygon/sampleTrailClusterRatio))
    trialPerCluster <- round(sqrt(samplePerPolygon*sampleTrailClusterRatio))
    list(
        polyDF=samplePolygons(
            unitSim,
            N = clustersPerPolygon,
            M = trialPerCluster,
            p = 1, # sample all polygons
            rWidth = x),
        pointDF=NULL)})

# sample with points only
pointDF <- samplePoints(unitSim, pointSamples, M)

# sample with mix of nonoverlapping
mixDFList <- lapply(rWidthSamples, function(x){
    samplePerPolygon <- totalTrials / x^2
    sampleTrailClusterRatio <- pointSamples / M
    clustersPerPolygon <- round(sqrt(samplePerPolygon/sampleTrailClusterRatio))
    trialPerCluster <- round(sqrt(samplePerPolygon*sampleTrailClusterRatio))
    samplePPMix(
        unitSim,
        N = clustersPerPolygon,
        M = trialPerCluster,
        p = .5, # sample all polygons
        rWidth = x)})

# sample with mix overlapping units
ovDFList <- lapply(polyDFList, function(x){
    DF <- x$polyDF
    DF$known <- as.logical(rbinom(nrow(DF), 1, .5))
    list(
        pointDF = DF %>%
            filter(known) %>%
            select(-id) %>%
            rename(id=trueid) %>%
            select(-known),
        polyDF = DF %>%
            filter(!known) %>%
            select(-known)
    )
})

allSamplesList <- list(
    mixed=mixDFList,
    overlap=ovDFList
)

models <- c(
    "Mixed Model" = 0,
    "Resample"    = 1,
    "Utazi"       = 2,
     #  "Reimann"     = 3,
    "Ignore"      = 4,
    "Known"       = 5)

# run all the models!!!
unitModelList <- lapply(models, function(i){
    lapply(allSamplesList, function(p){
        lapply(rWidthSamples, function(np){
            ps <- p[[as.character(np)]]
            runFieldModel(
                unitSim,
                pointDF=ps$pointDF,
                polyDF=ps$polyDF,
                verbose=T,
                control=list(),
                moption=i,
                rWidth=as.numeric(np)
            )
        })
    })
})

unitFitList <- lapply(unitModelList, function(model){
    lapply(model, function(sampleType){
        lapply(sampleType, function(ffit){
            simulateFieldCI(unitSim, ffit)
        })
    })
})

unitBetaList <- lapply(unitModelList, function(model){
    lapply(model, function(sampleType){
        lapply(sampleType, function(ffit){
            simulateBetas(unitSim, ffit)
        })
    })
})

convergeList <- lapply(unitModelList, function(model){
    lapply(model, function(sampleType){
        lapply(sampleType, function(ffit){
            ffit$opt$convergence
        })
    })
})

timeList <- lapply(unitModelList, function(model){
    lapply(model, function(sampleType){
        lapply(sampleType, function(ffit){
            as.numeric(ffit$runtime, units="mins")
        })
    })
})

unitResults <- list(
    sim = unitSim,
    pred = unitFitList,
    betas = unitBetaList,
    # model = unitModelList, # this takes up a lot of space so ignore for now
    converge = convergeList,
    covType = covType,
    rangeE = rangeE,
    covVal = covVal,
    M = M,
    seed = seed,
    runtime = timeList
)

saveRDS(unitResults, file=paste0("/gscratch/csde/nmarquez/Data/utaziTest3/", modelname))
