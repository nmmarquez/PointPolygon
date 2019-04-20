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

rdsPathList <- list.files("~/Data/utaziTest2/", full.names=TRUE)

resultsDF <- bind_rows(mclapply(rdsPathList, function(f_){
    x <- readRDS(f_)
    pList <- unlist(unlist(x$pred, recursive=FALSE), recursive=FALSE)
    bList <- unlist(unlist(x$betas, recursive=FALSE), recursive=FALSE)
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
        b1ERR = sapply(bList, function(b) b$betaStErr[2]),
        model = str_split(names(pList), "\\.", simplify=TRUE)[,1],
        sampling = str_split(names(pList), "\\.", simplify=TRUE)[,2],
        polysize = str_split(names(pList), "\\.", simplify=TRUE)[,3],
        converge = unlist(x$converge))
    }, mc.cores=5))

resultsDF %>%
    mutate(b1Bias = abs(b1Bias)) %>%
    #filter(model != "point") %>%
    #filter(model %in% c("riemann", "utazi")) %>%
    filter(converge == 0) %>%
    group_by(covType, model, rangeE) %>%
    summarize(b1Cov=mean(b1Cov)) %>%
    arrange(covType, rangeE, -b1Cov) %>%
    group_by(covType, rangeE) %>%
    mutate(rank=1:n()) %>%
    as.data.frame

resultsDF %>%
    filter(model != "Known") %>%
    mutate(Model=str_to_title(model)) %>%
    filter(converge == 0) %>%
    group_by(covType, rangeE, Model) %>%
    summarize(
        mu = mean(coverage),
        lwr = quantile(coverage, probs=.025),
        upr = quantile(coverage, probs=.975)
    ) %>%
    ggplot(aes(x=Model, ymin=lwr, y=mu, ymax=upr, color=Model)) +
    geom_point() +
    geom_errorbar() +
    theme_classic() +
    facet_grid(rangeE~covType) +
    coord_flip() +
    geom_hline(yintercept=.95, linetype=2) +
    labs(x="Model", y="") +
    ggtitle("95% Coverage of Underlying Probability Field") +
    theme(panel.spacing.y = unit(0, "lines")) +
    guides(color=FALSE)

resultsDF %>%
    filter(model != "Known") %>%
    mutate(Model=str_to_title(model)) %>%
    filter(converge == 0) %>%
    group_by(covType, rangeE, Model) %>%
    summarize(
        mu = mean(rmse),
        lwr = quantile(rmse, probs=.025),
        upr = quantile(rmse, probs=.975)
    ) %>%
    ggplot(aes(x=Model, ymin=lwr, y=mu, ymax=upr, color=Model)) +
    geom_point() +
    geom_errorbar() +
    theme_classic() +
    facet_grid(rangeE~covType) +
    coord_flip() +
    geom_hline(yintercept=0, linetype=2) +
    labs(x="Model", y="") +
    ggtitle("RMSE plots") +
    theme(panel.spacing.y = unit(0, "lines")) +
    guides(color=FALSE)

resultsDF %>%
    filter(model=="Utazi") %>%
    select(covType:rmse, sampling, polysize) %>%
    rename(rmseUtazi=rmse) %>%
    right_join(select(resultsDF, covType:rmse, model, sampling, converge, polysize)) %>%
    filter(converge == 0 & model != "Known") %>%
    mutate(improveRatio=rmseUtazi/rmse) %>%
    group_by(covType, model, rangeE) %>%
    summarize(
        mu = mean(improveRatio),
        lwr = mean(improveRatio) - 1.96*(sd(improveRatio)/sqrt(n())),
        upr = mean(improveRatio) + 1.96*(sd(improveRatio)/sqrt(n()))) %>%
    ungroup %>%
    rename(Model=model) %>%
    mutate(txt=round(mu, 2)) %>%
    ggplot(aes(x=Model, ymin=lwr, y=mu, ymax=upr, color=Model, label=txt)) +
    geom_point() +
    geom_errorbar() +
    theme_classic() +
    facet_grid(rangeE~covType) +
    coord_flip() +
    geom_hline(yintercept=1, linetype=2) +
    labs(x="Model", y="Relative Improvement") +
    ggtitle("RMSE: Margin of Improvement Over Utazi Model") +
    theme(panel.spacing.y = unit(0, "lines")) +
    guides(color=FALSE) +
    geom_text(nudge_y = .11)

resultsDF %>%
    mutate(cov.off = abs(.95 - coverage)) %>%
    #filter(model != "point") %>%
    #filter(model %in% c("riemann", "utazi")) %>%
    arrange(model) %>%
    group_by(covType, covVal, rangeE, M, seed, sampling) %>%
    summarize(
        any.fail = sum(converge),
        cov.off = abs(nth(cov.off, 2)) - abs(nth(cov.off, 1))
    ) %>%
    ungroup %>%
    filter(any.fail == 0) %>%
    select(cov.off) %>%
    unlist %>%
    (function(x) sd(x) / sqrt(length(x)))


diagnosticFilterDF <- bind_rows(
    resultsDF %>%
        #filter(sampling != "rwidth_3 poly") %>%
        mutate(b1Bias = abs(b1Bias)) %>%
        #filter(model != "point") %>%
        #filter(model %in% c("riemann", "utazi")) %>%
        arrange(model) %>%
        group_by(covType, covVal, rangeE, M, seed, sampling) %>%
        summarize(
            any.fail = sum(converge),
            b1Bias = abs(nth(b1Bias, 2)) - abs(nth(b1Bias, 1))
        ) %>%
        ungroup %>%
        filter(any.fail == 0) %>%
        mutate(riemann1 = b1Bias > 0) %>%
        group_by(covType, rangeE) %>%
        summarize(
            mu = mean(b1Bias),
            st.err = sd(b1Bias) / sqrt(n()),
            pmu = mean(riemann1),
            glm = list(glm(riemann1 ~ 1, family=binomial))) %>%
        ungroup %>%
        mutate(plwr=sapply(glm, function(g) arm::invlogit(confint(g))[1])) %>%
        mutate(pupr=sapply(glm, function(g) arm::invlogit(confint(g))[2])) %>%
        mutate(diagnostic="Beta Bias"),

    resultsDF %>%
        filter(sampling != "rwidth_3 poly") %>%
        mutate(cov.off = abs(.95 - coverage)) %>%
        filter(model != "point") %>%
        filter(model %in% c("riemann", "utazi")) %>%
        arrange(model) %>%
        group_by(covType, covVal, rangeE, M, seed, sampling) %>%
        summarize(
            any.fail = sum(converge),
            cov.off = abs(nth(cov.off, 2)) - abs(nth(cov.off, 1))
        ) %>%
        ungroup %>%
        filter(any.fail == 0) %>%
        mutate(riemann1 = cov.off > 0) %>%
        group_by(covType, rangeE) %>%
        summarize(
            mu = mean(cov.off),
            st.err = sd(cov.off) / sqrt(n()),
            pmu = mean(riemann1),
            glm = list(glm(riemann1 ~ 1, family=binomial))) %>%
        ungroup %>%
        mutate(plwr=sapply(glm, function(g) arm::invlogit(confint(g))[1])) %>%
        mutate(pupr=sapply(glm, function(g) arm::invlogit(confint(g))[2])) %>%
        mutate(diagnostic="Coverage"),

    resultsDF %>%
        filter(sampling != "rwidth_3 poly") %>%
        filter(model != "point") %>%
        filter(model %in% c("riemann", "utazi")) %>%
        arrange(model) %>%
        group_by(covType, covVal, rangeE, M, seed, sampling) %>%
        summarize(
            any.fail = sum(converge),
            rmse = nth(rmse, 2) - nth(rmse, 1)
        ) %>%
        ungroup %>%
        filter(any.fail == 0) %>%
        mutate(riemann1 = rmse > 0) %>%
        group_by(covType, rangeE) %>%
        summarize(
            mu = mean(rmse),
            st.err = sd(rmse) / sqrt(n()),
            pmu = mean(riemann1),
            glm = list(glm(riemann1 ~ 1, family=binomial))) %>%
        ungroup %>%
        mutate(diagnostic="RMSE")) %>%
    ungroup %>%
    mutate(plwr=sapply(glm, function(g) arm::invlogit(confint(g))[1])) %>%
    mutate(pupr=sapply(glm, function(g) arm::invlogit(confint(g))[2])) %>%
    select(-glm) %>%
    mutate(lwr=mu - 1.98*`st.err`) %>%
    mutate(upr=mu + 1.98*`st.err`)

diagnosticDF <- bind_rows(
    resultsDF %>%
        mutate(b1Bias = abs(b1Bias)) %>%
        filter(model != "point") %>%
        filter(model %in% c("riemann", "utazi")) %>%
        arrange(model) %>%
        group_by(covType, covVal, rangeE, M, seed, sampling) %>%
        summarize(
            any.fail = sum(converge),
            b1Bias = abs(nth(b1Bias, 2)) - abs(nth(b1Bias, 1))
        ) %>%
        ungroup %>%
        filter(any.fail == 0) %>%
        mutate(riemann1 = b1Bias > 0) %>%
        group_by(covType, rangeE) %>%
        summarize(
            mu = mean(b1Bias),
            st.err = sd(b1Bias) / sqrt(n()),
            pmu = mean(riemann1),
            glm = list(glm(riemann1 ~ 1, family=binomial))) %>%
        ungroup %>%
        mutate(plwr=sapply(glm, function(g) arm::invlogit(confint(g))[1])) %>%
        mutate(pupr=sapply(glm, function(g) arm::invlogit(confint(g))[2])) %>%
        mutate(diagnostic="Beta Bias"),
    
    resultsDF %>%
        mutate(cov.off = abs(.95 - coverage)) %>%
        filter(model != "point") %>%
        filter(model %in% c("riemann", "utazi")) %>%
        arrange(model) %>%
        group_by(covType, covVal, rangeE, M, seed, sampling) %>%
        summarize(
            any.fail = sum(converge),
            cov.off = abs(nth(cov.off, 2)) - abs(nth(cov.off, 1))
        ) %>%
        ungroup %>%
        filter(any.fail == 0) %>%
        mutate(riemann1 = cov.off > 0) %>%
        group_by(covType, rangeE) %>%
        summarize(
            mu = mean(cov.off),
            st.err = sd(cov.off) / sqrt(n()),
            pmu = mean(riemann1),
            glm = list(glm(riemann1 ~ 1, family=binomial))) %>%
        ungroup %>%
        mutate(plwr=sapply(glm, function(g) arm::invlogit(confint(g))[1])) %>%
        mutate(pupr=sapply(glm, function(g) arm::invlogit(confint(g))[2])) %>%
        mutate(diagnostic="Coverage"),
    
    resultsDF %>%
        filter(model != "point") %>%
        filter(model %in% c("riemann", "utazi")) %>%
        arrange(model) %>%
        group_by(covType, covVal, rangeE, M, seed, sampling) %>%
        summarize(
            any.fail = sum(converge),
            rmse = nth(rmse, 2) - nth(rmse, 1)
        ) %>%
        ungroup %>%
        filter(any.fail == 0) %>%
        mutate(riemann1 = rmse > 0) %>%
        group_by(covType, rangeE) %>%
        summarize(
            mu = mean(rmse),
            st.err = sd(rmse) / sqrt(n()),
            pmu = mean(riemann1),
            glm = list(glm(riemann1 ~ 1, family=binomial))) %>%
        ungroup %>%
        mutate(diagnostic="RMSE")) %>%
    ungroup %>%
    mutate(plwr=sapply(glm, function(g) arm::invlogit(confint(g))[1])) %>%
    mutate(pupr=sapply(glm, function(g) arm::invlogit(confint(g))[2])) %>%
    select(-glm) %>%
    mutate(lwr=mu - 1.98*`st.err`) %>%
    mutate(upr=mu + 1.98*`st.err`)

diagnosticDF %>%
    filter(diagnostic != "Coverage") %>%
    mutate(diagnostic=paste0(diagnostic, " Improvement")) %>%
    mutate(rangeE = as.character(rangeE)) %>%
    ggplot(aes(x=rangeE, y=mu, ymin=lwr, ymax=upr)) +
    facet_grid(covType ~ diagnostic, scales="free_x") +
    geom_point() +
    geom_errorbar() +
    coord_flip() +
    theme_classic() +
    geom_hline(yintercept=0, linetype=2) +
    labs(x="Spatial Range", y="") +
    ggtitle("Margin of Improved Result Using Riemann") +
    theme(panel.spacing.y = unit(0, "lines"))

diagnosticDF %>%
    mutate(rangeE = as.character(rangeE)) %>%
    ggplot(aes(x=rangeE, y=pmu, ymin=plwr, ymax=pupr)) +
    facet_grid(covType ~ diagnostic, scales="free_x") +
    geom_point() +
    geom_errorbar() +
    coord_flip() +
    theme_classic() +
    geom_hline(yintercept=0.5, linetype=2) +
    labs(x="Spatial Range", y="") +
    ggtitle("Probability of Improved Result Using Riemann") +
    theme(panel.spacing.y = unit(0, "lines"))

diagnosticFilterDF %>%
    filter(diagnostic != "Coverage") %>%
    mutate(diagnostic=paste0(diagnostic, " Improvement")) %>%
    mutate(rangeE = as.character(rangeE)) %>%
    ggplot(aes(x=rangeE, y=mu, ymin=lwr, ymax=upr)) +
    facet_grid(covType ~ diagnostic, scales="free_x") +
    geom_point() +
    geom_errorbar() +
    coord_flip() +
    theme_classic() +
    geom_hline(yintercept=0, linetype=2) +
    labs(x="Spatial Range", y="") +
    ggtitle("Margin of Improved Result Using Riemann: Polygon only 3x3 removed") +
    theme(panel.spacing.y = unit(0, "lines"))

diagnosticFilterDF %>%
    mutate(rangeE = as.character(rangeE)) %>%
    ggplot(aes(x=rangeE, y=pmu, ymin=plwr, ymax=pupr)) +
    facet_grid(covType ~ diagnostic, scales="free_x") +
    geom_point() +
    geom_errorbar() +
    coord_flip() +
    theme_classic() +
    geom_hline(yintercept=0.5, linetype=2) +
    labs(x="Spatial Range", y="") +
    ggtitle("Probability of Improved Result Using Riemann: Polygon only 3x3 removed") +
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
    ggplot(aes(value)) +
    geom_density() +
    theme_classic() +
    facet_wrap(~key, scales = "free")

table(resultsDF$covType)
aggResDF <- resultsDF %>%
    filter(converge == 0 & model != "point") %>%
    group_by(sampling, covType, rangeE, model) %>%
    summarize(
      mean(coverage), 
      min(coverage), 
      mean(rmse), 
      max(rmse),
      mean(correlation),
      min(correlation),
      `mean(b0)`=mean(b0Cov, na.rm=T),
      `mean(b1)`=mean(b1Cov, na.rm=T)) %>%
    as.data.frame %>%
    arrange(covType, rangeE, sampling, -`mean(coverage)`)

write_csv(aggResDF, "~/Data/utaziResults/aggRes.csv")
write_csv(resultsDF, "~/Data/utaziResults/results.csv")
