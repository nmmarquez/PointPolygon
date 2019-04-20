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

mixSample <- samplePPMix(field, 20, 100, .5, rWidth=3)

model0 <- runFieldModel(
    field, mixSample$pointDF, mixSample$polyDF, moption=0, verbose=T)
model1 <- runFieldModel(
    field, mixSample$pointDF, mixSample$polyDF, moption=1, verbose=T)
model2 <- runFieldModel(
    field, mixSample$pointDF, mixSample$polyDF, moption=2, verbose=T, rWidth = 3)
model3 <- runFieldModel(
    field, mixSample$pointDF, mixSample$polyDF, moption=3, verbose=T)
model4 <- runFieldModel(
    field, mixSample$pointDF, mixSample$polyDF, moption=4, verbose=T)
model5 <- runFieldModel(
    field, mixSample$pointDF, mixSample$polyDF, moption=5, verbose=T)

mcmcmodel5 <- runFieldModel(
    field, mixSample$pointDF, mixSample$polyDF, moption=5, mcmc=T, chains=1)
mcmcmodel3 <- runFieldModel(
    field, mixSample$pointDF, mixSample$polyDF, moption=0, mcmc=T, chains=1)

# this model runs in 7 seconds
# model2 <- runFieldModel(unitSim, mixSample$pointDF, mixSample$polyDF, moption=0)
# 
modelPredList <- lapply(list(model1, model2, model3), function(x){
    simulateFieldCI(field, x)
})

ggFieldEst(field, modelPredList)
ggFieldEst(field, modelPredList, sd = T)