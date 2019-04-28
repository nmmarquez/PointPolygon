.libPaths(c("~/R3.5/", .libPaths()))
rm(list=ls())
library(PointPolygon)
library(dplyr)
library(Rcpp)
set.seed(1234)

field <- simField(
    N = 100, rangeE = .7,
    offset = c(0.1, 0.2), 
    max.edge = c(0.3,0.4),
    beta0 = -2,
    betaList = list(list(type="spatial", value=1)))

ggField(field)

mixSample <- samplePPMix(field, 40, 100, p=.5, rWidth=5)

modelFits <- lapply(0:5, function(m){
    runFieldModel(
        field, mixSample$pointDF, mixSample$polyDF, moption=m, verbose=T,
        rWidth = 5)
})

modelrez <- lapply(modelFits, function(x){
    simulateFieldCI(field, x)
})

sapply(modelrez, function(df){
    sqrt(mean((df$mu - df$trueValue)^2))
})

mcmcmodel5 <- runFieldModel(
    field, mixSample$pointDF, mixSample$polyDF, moption=5, mcmc=T, chains=1)
mcmcmodel0 <- runFieldModel(
    field, mixSample$pointDF, mixSample$polyDF, moption=0, mcmc=T, chains=1)

# this model runs in 7 seconds
# model2 <- runFieldModel(unitSim, mixSample$pointDF, mixSample$polyDF, moption=0)
# 


ggFieldEst(field, modelPredList)
ggFieldEst(field, modelPredList, sd = T)