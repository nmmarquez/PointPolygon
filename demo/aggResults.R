rm(list=ls())
library(tibble)
library(dplyr)
library(parallel)
library(readr)

rdsPathList <- list.files("~/Data/utaziTest", full.names=TRUE)

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
      min(correlation)) %>%
    as.data.frame %>%
    arrange(covType, rangeE, -`mean(coverage)`)

write_csv(aggResDF, "~/Data/utaziResults/aggRes.csv")
write_csv(resultsDF, "~/Data/utaziResults/results.csv")
