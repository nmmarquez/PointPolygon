rm(list=ls())
library(TMB)
library(PointPolygon)
library(dplyr)

field <- simField(
    N = 500, rangeE = .7,
    offset = c(0.1, 0.2), 
    max.edge = c(0.1,0.2),
    beta0 = -2,
    betaList = list())

mixSample <- samplePPMix(field, 20, 100, .5, rWidth=5)
fullDF <- mixSample$polyDF %>%
    select(-id) %>%
    rename(id=trueid) %>%
    bind_rows(mixSample$pointDF)
model1 <- runFieldModel(field, fullDF, verbose = T)
model2 <- runFieldModel(
    field, mixSample$pointDF, mixSample$polyDF, moption=0, verbose=T)

# this model runs in 7 seconds
# model2 <- runFieldModel(unitSim, mixSample$pointDF, mixSample$polyDF, moption=0)
# 
# modelPredList <- lapply(list(model1, model2), function(x){
#     simulateFieldCI(unitSim, x)
# })

# model <- "PointPolygon"
# dyn.load("PointPolygon.so")
# 
# inputs <- buildModelInputs(
#     field,
#     pointDF = mixSample$pointDF,
#     polyDF = mixSample$polyDF,
#     moption = 0)
# 
# Obj <- TMB::MakeADFun(
#     data = inputs$Data,
#     parameters = inputs$Params,
#     DLL = model,
#     random = "z")
# 
# Opt <- stats::nlminb(
#     start = Obj$par,
#     objective = Obj$fn,
#     gradient = Obj$gr)
# sdrep <- TMB::sdreport(Obj, getJointPrecision=TRUE)
