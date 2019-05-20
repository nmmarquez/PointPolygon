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
library(sf)

rdsPathList <- list.files("~/Data/spaceTimeTest/", full.names=TRUE)

resultsDF <- bind_rows(mclapply(rdsPathList, function(f_){
    print(f_)
    dz <- tryCatch({
        x <- readRDS(f_)
        pList <- x$pred
        bList <- x$betas
        tibble(
            covType = x$covType,
            covVal = x$covVal,
            rangeE = x$rangeE,
            seed = x$seed,
            rmse = sapply(pList, function(y){
                sqrt(mean((y$trueValue - y$mu)^2))}),
            provrmse = sapply(x$provPred, function(y){
                sqrt(mean((y$trueValue - y$mu)^2))}),
            bias = sapply(pList, function(y) mean(y$mu - y$trueValue)),
            dissDiff = sapply(pList, function(y){
                (.5*mean(abs(y$mu / mean(y$mu) - (1-y$mu) / mean(1-y$mu)))) - 
                    (.5*mean(abs(y$trueValue / mean(y$trueValue) - 
                                     (1-y$trueValue) / mean(1-y$trueValue))))
            }),
            coverage = sapply(pList, function(y){
                mean(y$trueValue >= y$lwr & y$trueValue <= y$upr)}),
            provcoverage = sapply(x$provPred, function(y){
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
            model = names(pList),
            converge = unlist(x$converge),
            runtime = unlist(x$runtime))
    },
    error= function(cond){
        tibble()
    })
    dz}, mc.cores=8)) %>%
    mutate(model=gsub("Reimann", "Riemann", model))

aggPlotsDR <- list()

(aggPlotsDR$coverage <- resultsDF %>%
    mutate(Model=model) %>%
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
    guides(color=FALSE))

(aggPlotsDR$coveragePaper <- resultsDF %>%
        mutate(Model=model) %>%
        filter(converge == 0 & model != "Riemann") %>%
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
        guides(color=FALSE))

(aggPlotsDR$provcoverage <- resultsDF %>%
        mutate(Model=model) %>%
        filter(converge == 0) %>%
        group_by(covType, rangeE, Model) %>%
        summarize(
            mu = mean(provcoverage),
            lwr = quantile(provcoverage, probs=.025),
            upr = quantile(provcoverage, probs=.975)
        ) %>%
        ggplot(aes(x=Model, ymin=lwr, y=mu, ymax=upr, color=Model)) +
        geom_point() +
        geom_errorbar() +
        theme_classic() +
        facet_grid(rangeE~covType) +
        coord_flip() +
        geom_hline(yintercept=.95, linetype=2) +
        labs(x="Model", y="") +
        ggtitle("95% Coverage of Province Probability(N=32)") +
        theme(panel.spacing.y = unit(0, "lines")) +
        guides(color=FALSE))

(aggPlotsDR$provcoveragePaper <- resultsDF %>%
        mutate(Model=model) %>%
        filter(converge == 0 & model != "Riemann") %>%
        group_by(covType, rangeE, Model) %>%
        summarize(
            mu = mean(provcoverage),
            lwr = quantile(provcoverage, probs=.025),
            upr = quantile(provcoverage, probs=.975)
        ) %>%
        ggplot(aes(x=Model, ymin=lwr, y=mu, ymax=upr, color=Model)) +
        geom_point() +
        geom_errorbar() +
        theme_classic() +
        facet_grid(rangeE~covType) +
        coord_flip() +
        geom_hline(yintercept=.95, linetype=2) +
        labs(x="Model", y="") +
        ggtitle("95% Coverage of Province Probability(N=32)") +
        theme(panel.spacing.y = unit(0, "lines")) +
        guides(color=FALSE))

(aggPlotsDR$rmseRelative <- resultsDF %>%
    filter(model=="IHME Resample") %>%
    select(covType:rmse) %>%
    rename(rmseUtazi=rmse) %>%
    right_join(select(resultsDF, covType:rmse, model, converge)) %>%
    filter(converge == 0) %>%
    mutate(improveRatio=(rmseUtazi-rmse)/rmseUtazi) %>%
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
    geom_hline(yintercept=0, linetype=2) +
    labs(x="Model", y="Relative Improvement") +
    ggtitle("RMSE: Margin of Improvement Over IHME Resample Model") +
    theme(panel.spacing.y = unit(0, "lines")) +
    guides(color=FALSE) +
    geom_text(aes(y=upr), nudge_y = .04))

(aggPlotsDR$rmseRelativePaper <- resultsDF %>%
        filter(model=="IHME Resample") %>%
        select(covType:rmse) %>%
        rename(rmseUtazi=rmse) %>%
        right_join(select(resultsDF, covType:rmse, model, converge)) %>%
        filter(converge == 0 & model != "Riemann") %>%
        mutate(improveRatio=(rmseUtazi-rmse)/rmseUtazi) %>%
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
        geom_hline(yintercept=0, linetype=2) +
        labs(x="Model", y="Relative Improvement") +
        ggtitle("RMSE: Margin of Improvement Over IHME Resample Model") +
        theme(panel.spacing.y = unit(0, "lines")) +
        guides(color=FALSE) +
        geom_text(aes(y=upr), nudge_y = .04))

(aggPlotsDR$rmseProvRelative <- resultsDF %>%
        filter(model=="IHME Resample") %>%
        select(covType:seed, provrmse) %>%
        rename(rmseUtazi=provrmse) %>%
        right_join(select(resultsDF, covType:seed, provrmse, model, converge)) %>%
        filter(converge == 0) %>%
        mutate(improveRatio=(rmseUtazi-provrmse)/rmseUtazi) %>%
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
        geom_hline(yintercept=0, linetype=2) +
        labs(x="Model", y="Relative Improvement") +
        ggtitle("Province RMSE: Margin of Improvement Over IHME Resample Model") +
        theme(panel.spacing.y = unit(0, "lines")) +
        guides(color=FALSE) +
        geom_text(aes(y=upr), nudge_y = .06))

(aggPlotsDR$rmseProvRelativePaper <- resultsDF %>%
        filter(model=="IHME Resample") %>%
        select(covType:seed, provrmse) %>%
        rename(rmseUtazi=provrmse) %>%
        right_join(select(resultsDF, covType:seed, provrmse, model, converge)) %>%
        filter(converge == 0 & model != "Riemann") %>%
        mutate(improveRatio=(rmseUtazi-provrmse)/rmseUtazi) %>%
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
        geom_hline(yintercept=0, linetype=2) +
        labs(x="Model", y="Relative Improvement") +
        ggtitle("Province RMSE: Margin of Improvement Over IHME Resample Model") +
        theme(panel.spacing.y = unit(0, "lines")) +
        guides(color=FALSE) +
        geom_text(aes(y=upr), nudge_y = .06))

(aggPlotsDR$bias <- resultsDF %>%
    filter(converge == 0) %>%
    group_by(covType, model, rangeE) %>%
    summarize(
        mu = mean(bias),
        lwr = quantile(bias, probs=.025),
        upr = quantile(bias, probs=.975)) %>%
    ungroup %>%
    rename(Model=model) %>%
    mutate(txt=round(mu, 2)) %>%
    ggplot(aes(x=Model, ymin=lwr, y=mu, ymax=upr, color=Model, label=txt)) +
    geom_point() +
    geom_errorbar() +
    theme_classic() +
    facet_grid(rangeE~covType) +
    coord_flip() +
    geom_hline(yintercept=0, linetype=2) +
    labs(x="Model", y="Bias") +
    ggtitle("RMSE: Average Bias of Models") +
    theme(panel.spacing.y = unit(0, "lines")) +
    guides(color=FALSE))

(aggPlotsDR$biasPaper <- resultsDF %>%
        filter(converge == 0 & model != "Riemann") %>%
        group_by(covType, model, rangeE) %>%
        summarize(
            mu = mean(bias),
            lwr = quantile(bias, probs=.025),
            upr = quantile(bias, probs=.975)) %>%
        ungroup %>%
        rename(Model=model) %>%
        mutate(txt=round(mu, 2)) %>%
        ggplot(aes(x=Model, ymin=lwr, y=mu, ymax=upr, color=Model, label=txt)) +
        geom_point() +
        geom_errorbar() +
        theme_classic() +
        facet_grid(rangeE~covType) +
        coord_flip() +
        geom_hline(yintercept=0, linetype=2) +
        labs(x="Model", y="Bias") +
        ggtitle("RMSE: Average Bias of Models") +
        theme(panel.spacing.y = unit(0, "lines")) +
        guides(color=FALSE))

(aggPlotsDR$dissDiff <- resultsDF %>%
    filter(converge == 0) %>%
    group_by(covType, model, rangeE) %>%
    summarize(
        mu = mean(dissDiff),
        lwr = quantile(dissDiff, probs=.025),
        upr = quantile(dissDiff, probs=.975)) %>%
    ungroup %>%
    rename(Model=model) %>%
    mutate(txt=round(mu, 2)) %>%
    ggplot(aes(x=Model, ymin=lwr, y=mu, ymax=upr, color=Model, label=txt)) +
    geom_point() +
    geom_errorbar() +
    theme_classic() +
    facet_grid(rangeE~covType) +
    coord_flip() +
    geom_hline(yintercept=0, linetype=2) +
    labs(x="Model", y="Bias") +
    ggtitle("Dissimilarity Difference") +
    theme(panel.spacing.y = unit(0, "lines")) +
    guides(color=FALSE))

(aggPlotsDR$dissDiffPaper <- resultsDF %>%
        filter(converge == 0 & model != "Riemann") %>%
        group_by(covType, model, rangeE) %>%
        summarize(
            mu = mean(dissDiff),
            lwr = quantile(dissDiff, probs=.025),
            upr = quantile(dissDiff, probs=.975)) %>%
        ungroup %>%
        rename(Model=model) %>%
        mutate(txt=round(mu, 2)) %>%
        ggplot(aes(x=Model, ymin=lwr, y=mu, ymax=upr, color=Model, label=txt)) +
        geom_point() +
        geom_errorbar() +
        theme_classic() +
        facet_grid(rangeE~covType) +
        coord_flip() +
        geom_hline(yintercept=0, linetype=2) +
        labs(x="Model", y="Bias") +
        ggtitle("Dissimilarity Difference") +
        theme(panel.spacing.y = unit(0, "lines")) +
        guides(color=FALSE))

(aggPlotsDR$runtime <- resultsDF %>%
    select(covType:seed, model, runtime) %>%
    filter(model %in% c("Mixture Model", "Riemann")) %>%
    left_join(
        resultsDF %>%
            filter(model %in% c("IHME Resample")) %>%
            select(covType:seed, runtime) %>%
            rename(IHME=runtime), 
        by=c("covType", "covVal", "rangeE", "seed")) %>%
    ggplot(aes(x=IHME, y=runtime)) +
    geom_point() +
    geom_abline() +
    theme_classic() +
    facet_wrap(~model, scales="free_y") +
    labs(x="IHME Resample Runtime", y="Runtime"))

write_rds(aggPlotsDR, "~/Documents/PointPolygon/demo/aggplotsDR.Rds")
