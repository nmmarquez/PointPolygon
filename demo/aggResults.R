.libPaths(c("~/R3.5/", .libPaths()))
rm(list=ls())
library(tibble)
library(dplyr)
library(parallel)
library(readr)
library(PointPolygon)
library(stringr)

rdsPathList <- list.files("~/Data/utaziTest", full.names=TRUE)

smallResults <- mclapply(rdsPathList, function(f_){
    x <- readRDS(f_)
    results <- list(
        covType = x$covType,
        covVal = x$covVal,
        M = x$M,
        seed = x$seed
    )
    results$model <- lapply(x$model, function(m){
        list(hess=m$sd$cov.fixed, est=m$sd$par.fixed)
    })
    results
}, mc.cores=5)

saveRDS(smallResults, "~/Data/utaziResults/smallResults.Rds")

resultsDF <- bind_rows(mclapply(rdsPathList, function(f_){
    x <- readRDS(f_)
    pList <- unlist(x$pred, recursive=FALSE)
    bList <- unlist(x$betas, recursive=FALSE)
    tibble(
        covType = x$covType,
        covVal = x$covVal,
        rangeE = x$rangeE,
        M = x$M,
        seed = x$seed,
        rmse = sapply(pList, function(y) sqrt(mean((y$trueValue - y$mu)^2))),
        coverage = sapply(pList, function(y){
            mean(y$trueValue >= y$lwr & y$trueValue <= y$upr)}),
        correlation = sapply(pList, function(y) cor(y$mu, y$trueValue)),
        b0Cov = sapply(bList, function(b){
            bhat <- b$betaHat[1]
            sder <- b$betaStErr[1]
            as.numeric(((bhat - sder) <= -2) & ((bhat + sder) >= -2))
        }),
        b1Cov = sapply(bList, function(b){
          bhat <- b$betaHat[2]
          sder <- b$betaStErr[2]
          as.numeric(((bhat - sder) <= x$covVal) & ((bhat + sder) >= x$covVal))
        }),
        b0Bias = sapply(bList, function(b){
          bhat <- b$betaHat[1]
          bhat + 2
        }),
        b1Bias = sapply(bList, function(b){
          bhat <- b$betaHat[2]
          bhat - x$covVal
        }),
        model = str_split(names(pList), "\\.", simplify=TRUE)[,1],
        sampling = str_split(names(pList), "\\.", simplify=TRUE)[,2])
    }, mc.cores=5))

resultsDF %>%
  filter(model %in% c("riemann", "utazi")) %>%
  arrange(model) %>%
  group_by(covType, covVal, rangeE, M, seed, sampling) %>%
  summarize(
    rmseDiff=diff(rmse),
    covDiff=diff(coverage),
    b1Bias=diff(b1Bias)) %>%
  summary


table(resultsDF$covType)
aggResDF <- resultsDF %>%
    filter(convergence) %>%
    group_by(type, covType, rangeE) %>%
    summarize(
      mean(coverage), 
      min(coverage), 
      mean(rmse), 
      max(rmse),
      mean(correlation),
      min(correlation),
      `mean(b0)`=mean(b0, na.rm=T),
      `mean(b1)`=mean(b1, na.rm=T)) %>%
    as.data.frame %>%
    arrange(covType, rangeE, -`mean(coverage)`)

write_csv(aggResDF, "~/Data/utaziResults/aggRes.csv")
write_csv(resultsDF, "~/Data/utaziResults/results.csv")
