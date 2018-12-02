.libPaths(c("~/R3.5/", .libPaths()))
rm(list=ls())
library(tibble)
library(dplyr)
library(parallel)
library(readr)
library(PointPolygon)
library(stringr)
library(tidyr)
library(ggplot2)

rdsPathList <- list.files("~/Data/utaziTest", full.names=TRUE)

isnt_out_z <- function(x, thres = 3, na.rm = TRUE) {
  abs(x - mean(x, na.rm = na.rm)) <= thres * sd(x, na.rm = na.rm)
}

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
        sampling = str_split(names(pList), "\\.", simplify=TRUE)[,2],
        converge = unlist(x$converge))
    }, mc.cores=5))

resultsDF %>%
  mutate(b1Bias = abs(b1Bias)) %>%
  filter(model != "point") %>%
  filter(model %in% c("riemann", "utazi")) %>%
  filter(converge == 0) %>%
  group_by(covType, model, rangeE) %>%
  summarize(b1Cov=mean(b1Cov)) %>%
  arrange(covType, rangeE, model)

diagnosticDF <- bind_rows(
    resultsDF %>%
        mutate(b1Bias = abs(b1Bias)) %>%
        filter(model != "point") %>%
        filter(model %in% c("riemann", "utazi")) %>%
        arrange(model) %>%
        group_by(covType, covVal, rangeE, M, seed, sampling) %>%
        mutate(any.fail=sum(converge)) %>%
        filter(b1Bias == min(b1Bias)) %>%
        ungroup %>%
        filter(any.fail == 0) %>%
        mutate(riemann1 = model == "riemann") %>%
        group_by(covType, rangeE) %>%
        summarize(
            mu = mean(riemann1),
            glm = list(glm(riemann1 ~ 1, family=binomial))
        ) %>%
        ungroup %>%
        mutate(lwr=sapply(glm, function(g) arm::invlogit(confint(g))[1])) %>%
        mutate(upr=sapply(glm, function(g) arm::invlogit(confint(g))[2])) %>%
        select(-glm) %>%
        mutate(diagnostic="Beta Bias"),

    resultsDF %>%
        mutate(cov.off = abs(.95 - coverage)) %>%
        filter(model != "point") %>%
        filter(model %in% c("riemann", "utazi")) %>%
        arrange(model) %>%
        group_by(covType, covVal, rangeE, M, seed, sampling) %>%
        mutate(any.fail=sum(converge)) %>%
        filter(cov.off == min(cov.off)) %>%
        ungroup %>%
        filter(any.fail == 0) %>%
        mutate(riemann1 = model == "riemann") %>%
        group_by(covType, rangeE) %>%
        summarize(
            mu = mean(riemann1),
            glm = list(glm(riemann1 ~ 1, family=binomial))
        ) %>%
        ungroup %>%
        mutate(lwr=sapply(glm, function(g) arm::invlogit(confint(g))[1])) %>%
        mutate(upr=sapply(glm, function(g) arm::invlogit(confint(g))[2])) %>%
        select(-glm) %>%
        mutate(diagnostic="Coverage"),

    resultsDF %>%
        mutate(cov.off = abs(.95 - coverage)) %>%
        filter(model != "point") %>%
        filter(model %in% c("riemann", "utazi")) %>%
        arrange(model) %>%
        group_by(covType, covVal, rangeE, M, seed, sampling) %>%
        mutate(any.fail=sum(converge)) %>%
        filter(rmse == min(rmse)) %>%
        ungroup %>%
        filter(any.fail == 0) %>%
        mutate(riemann1 = model == "riemann") %>%
        group_by(covType, rangeE) %>%
        summarize(
            mu = mean(riemann1),
            glm = list(glm(riemann1 ~ 1, family=binomial))
        ) %>%
        ungroup %>%
        mutate(lwr=sapply(glm, function(g) arm::invlogit(confint(g))[1])) %>%
        mutate(upr=sapply(glm, function(g) arm::invlogit(confint(g))[2])) %>%
        select(-glm) %>%
        mutate(diagnostic="RMSE"))

diagnosticDF %>%
    mutate(rangeE = as.character(rangeE)) %>%
    ggplot(aes(x=rangeE, y=mu, ymin=lwr, ymax=upr)) +
    facet_grid(covType ~ diagnostic) +
    geom_point() +
    geom_errorbar() +
    coord_flip() +
    theme_classic() +
    geom_hline(yintercept=.5, linetype=2) +
    labs(x="Spatial Range", y="") +
    ggtitle("Probability of Improved Result Using Riemann") +
    theme(panel.spacing.y = unit(0, "lines"))

aggRes <- resultsDF %>%
    filter(model %in% c("riemann", "utazi")) %>%
    arrange(model) %>%
    group_by(covType, covVal, rangeE, M, seed, sampling) %>%
    summarize(
        rmseDiff = diff(rmse),
        covDiff = diff(coverage),
        b1Bias = diff(abs(b1Bias)),
        converge = sum(converge)) %>%
    filter(converge == 0)

aggRes %>%
    ungroup %>%
    select(rmseDiff, covDiff) %>%
    gather("key", "value") %>%
    group_by(key) %>%
    mutate(not_out=isnt_out_z(value, thres = 2.2)) %>%
    mutate(value.zoom=ifelse(not_out, value, NA)) %>%
    ggplot(aes(value.zoom)) +
    geom_density() +
    theme_classic() +
    facet_wrap(~key, scales = "free")

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
