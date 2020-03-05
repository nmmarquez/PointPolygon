.libPaths(c("~/R3.5/", .libPaths()))
rm(list=ls())
library(PointPolygon)
library(dplyr)
library(Rcpp)
set.seed(1234)

field <- simField(
    N = 100, rangeE = .4,
    offset = c(0.1, 0.2),
    max.edge = c(0.1,0.2),
    beta0 = -2,
    #betaList = list(list(type="spatial", value=1)),
    nTimes = 6,
    rho = .85)

ggField(field)
#mixSample <- samplePPMix(field, 100, 100, p=.5, rWidth=5)


mixSample <- list(
    polyDF  = samplePolygons(field, 150, 50, p=1., rWidth=5),
    pointDF = samplePoints(field, 200, 50))

modelTypes <- c(
    "Mixture Model" = 0,
    "Resample"      = 1,
    # "Utazi"         = 2, # Currently cant run utazi on more than one time
    "Riemann"       = 3,
    "Ignore"        = 4,
    "Known"         = 5
)

modelFits <- lapply(modelTypes, function(m){
    runFieldModel(
        field, mixSample$pointDF, mixSample$polyDF, moption=m, verbose=T,
        rWidth = 5)
})

modelFits$Best <- modelFits$Resample
modelFits$Best$sd$par.random <- field$latent
modelFits$Best$opt$par <- c(beta=-2, field$fieldPars)

modelrez <- lapply(modelFits, function(x){
    simulateFieldCI(field, x)
})


sapply(modelrez, function(df){
    sqrt(mean((df$mu - df$trueValue)^2))
})

ggFieldEst(field, modelrez)




arealrez <- lapply(modelFits, function(x){
    arealCI(field, x, rWidth=5)
})

# Confidence intervals for pixels
sapply(modelrez, function(x){
    mean(x$trueValue > x$lwr & x$trueValue < x$upr)
})

sapply(modelrez, function(x){
    mean(abs(x$upr - x$lwr))
})

# Confidence intervals for areal units
sapply(arealrez, function(x){
    mean(x$trueValue > x$lwr & x$trueValue < x$upr)
})

sapply(arealrez, function(x){
    mean(abs(x$upr - x$lwr))
})


arealrez$Best %>%
    as_tibble() %>%
    select(-geometry) %>%
    left_join(as_tibble(reimDF)) %>%
    select(trueValue, trials, obs, tidx, polyid) %>%
    mutate(crude=obs/trials, absdiff=abs(trueValue-crude)) %>%
    arrange(-absdiff)


ggplot() + 
    geom_sf(
        data = arealrez$Riemann %>%
            mutate(covered=lwr < trueValue & upr > trueValue),
        aes(fill = covered)) +
    facet_wrap(~tidx) +
    #ggplot2::scale_fill_distiller(palette = "Spectral") +
    theme_classic()

ggplot() + 
    geom_sf(data = arealrez$Known, aes(fill = trueValue)) +
    facet_wrap(~tidx) +
    ggplot2::scale_fill_distiller(palette = "Spectral") +
    theme_void() +
    

ggFieldEst(field, modelrez)
ggFieldEst(field, modelrez, sd = T)

mcmcmodel5 <- runFieldModel(
    field, mixSample$pointDF, mixSample$polyDF, moption=5, mcmc=T, chains=1)
mcmcmodel0 <- runFieldModel(
    field, mixSample$pointDF, mixSample$polyDF, moption=0, mcmc=T, chains=1)

# this model runs in 7 seconds
# model2 <- runFieldModel(unitSim, mixSample$pointDF, mixSample$polyDF, moption=0)
# 


ggFieldEst(field, modelPredList)
ggFieldEst(field, modelPredList, sd = T)