rm(list=ls())
library(tibble)
library(dplyr)
library(parallel)
library(readr)

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
    tibble(
        covType = x$covType,
        covVal = x$covVal,
        rangeE = x$rangeE,
        M = x$M,
        seed = x$seed,
        rmse = sapply(x$pred, function(y) sqrt(mean((y$trueValue - y$mu)^2))),
        coverage = sapply(x$pred, function(y){
            mean(y$trueValue >= y$lwr & y$trueValue <= y$upr)}),
        correlation = sapply(x$pred, function(y) cor(y$mu, y$trueValue)),
        convergence = sapply(x$model, function(y) y$opt$convergence == 0),
        b0 = sapply(x$model, function(m){
          sder <- 1.96 * sqrt(m$sd$cov.fixed[1,1])
          bhat <- m$sd$par.fixed[1]
          as.numeric(((bhat - sder) <= -2) & ((bhat + sder) >= -2))
        }),
        b1 = sapply(x$model, function(m){
          sder <- 1.96 * sqrt(m$sd$cov.fixed[2,2])
          bhat <- m$sd$par.fixed[2]
          as.numeric(((bhat - sder) <= x$covVal) & ((bhat + sder) >= x$covVal))
        }),
        type = names(x$pred)
    )}, mc.cores=5))

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
