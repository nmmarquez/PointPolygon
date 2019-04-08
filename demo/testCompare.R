rm(list=ls())
library(PointPolygon)

unitSim <- simField(
    N = 500, rangeE = .7,
    offset = c(0.1, 0.2), 
    max.edge = c(0.03,0.06),
    beta0 = -2,
    betaList = list())

pointDF <- samplePoints(unitSim, 500, 100)
mixSample <- samplePPMix(unitSim, 30, 150, .5, rWidth=10)
model1 <- runFieldModel(unitSim, pointDF)
# this model runs in 7 seconds
model2 <- runFieldModel(unitSim, mixSample$pointDF, mixSample$polyDF, moption=1)

modelPredList <- lapply(list(model1, model2), function(x){
    simulateFieldCI(unitSim, x)
})

ggFieldEst(unitSim, modelPredList)
