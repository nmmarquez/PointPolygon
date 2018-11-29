rm(list=ls())
library(PointPolygon)
library(INLA)

unitSim <- simField(
    N = 60, rangeE = .7,
    offset = c(0.1, 0.2), 
    max.edge = c(0.1,0.2),
    beta0 = -2,
    betaList = list(list(type="random", value=2)))

mixSample <- samplePPMix(unitSim, 30, 150, .5, rWidth=3)

shape3 <- dividePolygon(unitSim$bound, 3)
