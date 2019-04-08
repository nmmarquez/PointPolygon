.libPaths(c("~/R3.5/", .libPaths()))
rm(list=ls())
library(PointPolygon)
library(dplyr)
set.seed(1234)

field <- simField(
    N = 500, rangeE = .7,
    offset = c(0.1, 0.1), 
    max.edge = c(0.05,0.2),
    beta0 = -2,
    betaList = list())

ggField(field)

mixSample <- samplePPMix(field, 20, 100, .5, rWidth=3)
fullDF <- mixSample$polyDF %>%
    select(-id) %>%
    rename(id=trueid) %>%
    bind_rows(mixSample$pointDF)

model1 <- runFieldModel(field, fullDF, verbose = T)
model2 <- runFieldModel(
    field, mixSample$pointDF, mixSample$polyDF, moption=3, verbose=T)
model3 <- runFieldModel(
    field, mixSample$pointDF, mixSample$polyDF, moption=0, verbose=T)

# this model runs in 7 seconds
# model2 <- runFieldModel(unitSim, mixSample$pointDF, mixSample$polyDF, moption=0)
# 
modelPredList <- lapply(list(model1, model2, model3), function(x){
    simulateFieldCI(field, x)
})

ggFieldEst(field, modelPredList)
ggFieldEst(field, modelPredList, sd = T)