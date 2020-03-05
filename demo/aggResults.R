.libPaths(c("~/R3.5/", .libPaths()))
rm(list=ls())
library(GGally)
library(tibble)
library(dplyr)
library(parallel)
library(readr)
library(PointPolygon)
library(stringr)
library(tidyr)
library(ggplot2)
library(forcats)

rdsPathList <- list.files("~/Data/utaziTest3/", full.names=TRUE)

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
        bias = sapply(pList, function(y) mean(y$mu - y$trueValue)),
        dissDiff = sapply(pList, function(y){
            (.5*mean(abs(y$mu / mean(y$mu) - (1-y$mu) / mean(1-y$mu)))) - 
                (.5*mean(abs(y$trueValue / mean(y$trueValue) - 
                                 (1-y$trueValue) / mean(1-y$trueValue))))
        }),
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
    }, mc.cores=5)) %>%
    mutate(model=gsub("Reimann", "Riemann", model)) %>%
    mutate(model=gsub("Known", "Unmasked", model)) %>%
    mutate(model=gsub("Mixed Model", "Mixture", model)) %>%
    mutate(model=gsub("Utazi", "Ecological", model)) %>%
    mutate(model = fct_relevel(
        model,
        "Ignore", "Resample", "Ecological", "Riemann", "Mixture", "Unmasked"))

aggPlots <- list()

(aggPlots$coverage <- resultsDF %>%
    mutate(Model=fct_rev(model)) %>%
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

(aggPlots$coveragePaper <- resultsDF %>%
        mutate(Model=fct_rev(model)) %>%
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
        theme_bw() +
        facet_grid(rangeE~covType) +
        coord_flip() +
        geom_hline(yintercept=.95, linetype=2) +
        labs(x="Model", y="") +
        ggtitle("95% Coverage of Underlying Probability Field") +
        theme(panel.spacing.y = unit(0, "lines")) +
        guides(color=FALSE) +
        theme(
            strip.text = element_text(size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            title = element_text(size=25),
            axis.text.x = element_text(size=12),
            axis.text.y = element_text(size=15),
            axis.title.x = element_text(size=20)))

ggsave(
    "demo/figures/covSim1.png", aggPlots$coveragePaper,
    width=400, height=280, units = "mm")

(aggPlots$rmseRelative <- resultsDF %>%
    filter(model=="Ecological") %>%
    select(covType:rmse, sampling, polysize) %>%
    rename(rmseUtazi=rmse) %>%
    right_join(select(resultsDF, covType:rmse, model, sampling, converge, polysize)) %>%
    filter(converge == 0) %>%
    mutate(improveRatio=(rmseUtazi-rmse)/rmseUtazi) %>%
    group_by(covType, model, rangeE) %>%
    summarize(
        mu = mean(improveRatio),
        lwr = mean(improveRatio) - 1.96*(sd(improveRatio)/sqrt(n())),
        upr = mean(improveRatio) + 1.96*(sd(improveRatio)/sqrt(n()))) %>%
    ungroup %>%
    mutate(Model=fct_rev(model)) %>%
    mutate(txt=round(mu, 2)) %>%
    ggplot(aes(x=Model, ymin=lwr, y=mu, ymax=upr, color=Model, label=txt)) +
    geom_point() +
    geom_errorbar() +
    theme_classic() +
    facet_grid(rangeE~covType) +
    coord_flip() +
    geom_hline(yintercept=0, linetype=2) +
    labs(x="Model", y="Relative Improvement") +
    ggtitle("RMSE: Margin of Improvement Over Ecological Model") +
    theme(panel.spacing.y = unit(0, "lines")) +
    guides(color=FALSE) +
    geom_text(nudge_y = .42))

(aggPlots$rmseRelativePaper <- resultsDF %>%
        filter(model=="Ecological") %>%
        select(covType:rmse, sampling, polysize) %>%
        rename(rmseUtazi=rmse) %>%
        right_join(select(resultsDF, covType:rmse, model, sampling, converge, polysize)) %>%
        filter(converge == 0 & model !="Riemann") %>%
        mutate(improveRatio=(rmseUtazi-rmse)/rmseUtazi) %>%
        group_by(covType, model, rangeE) %>%
        summarize(
            mu = mean(improveRatio),
            lwr = mean(improveRatio) - 1.96*(sd(improveRatio)/sqrt(n())),
            upr = mean(improveRatio) + 1.96*(sd(improveRatio)/sqrt(n()))) %>%
        ungroup %>%
        mutate(Model=fct_rev(model)) %>%
        mutate(txt=round(mu, 2)) %>%
        ggplot(aes(x=Model, ymin=lwr, y=mu, ymax=upr, color=Model, label=txt)) +
        geom_point() +
        geom_errorbar() +
        theme_bw() +
        facet_grid(rangeE~covType) +
        coord_flip() +
        geom_hline(yintercept=0, linetype=2) +
        labs(x="", y="Relative Improvement") +
        ggtitle("RMSE: Margin of Improvement Over Ecological Model") +
        theme(panel.spacing.y = unit(0, "lines")) +
        guides(color=FALSE) +
        geom_text(nudge_y = .5) +
        theme(
            strip.text = element_text(size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            title = element_text(size=25),
            axis.text.x = element_text(size=12),
            axis.text.y = element_text(size=15),
            axis.title.x = element_text(size=20))
        )

ggsave(
    "demo/figures/rmseSim1.png", aggPlots$rmseRelativePaper,
    width=400, height=280, units = "mm")

(aggPlots$rmseRelativePaperb <- resultsDF %>%
        filter(model=="Ecological") %>%
        select(covType:rmse, sampling, polysize) %>%
        rename(rmseUtazi=rmse) %>%
        right_join(select(resultsDF, covType:rmse, model, sampling, converge, polysize)) %>%
        filter(converge == 0 & model !="Riemann" & model != "Ignore") %>%
        mutate(improveRatio=(rmseUtazi-rmse)/rmseUtazi) %>%
        group_by(covType, model, rangeE) %>%
        summarize(
            mu = mean(improveRatio),
            lwr = mean(improveRatio) - 1.96*(sd(improveRatio)/sqrt(n())),
            upr = mean(improveRatio) + 1.96*(sd(improveRatio)/sqrt(n()))) %>%
        ungroup %>%
        mutate(Model=fct_rev(model)) %>%
        mutate(txt=round(mu, 2)) %>%
        ggplot(aes(x=Model, ymin=lwr, y=mu, ymax=upr, color=Model, label=txt)) +
        geom_point() +
        geom_errorbar() +
        theme_bw() +
        facet_grid(rangeE~covType) +
        coord_flip() +
        geom_hline(yintercept=0, linetype=2) +
        labs(x="", y="Relative Improvement") +
        ggtitle("RMSE: Margin of Improvement Over Ecological Model") +
        theme(panel.spacing.y = unit(0, "lines")) +
        guides(color=FALSE) +
        geom_text(nudge_y = .5) +
        theme(
            strip.text = element_text(size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            title = element_text(size=25),
            axis.text.x = element_text(size=12),
            axis.text.y = element_text(size=15),
            axis.title.x = element_text(size=20)))

ggsave(
    "demo/figures/rmseSim1b.png", aggPlots$rmseRelativePaperb,
    width=400, height=280, units = "mm")

(aggPlots$rmseSingleRelativePaper <- resultsDF %>%
        filter(model=="Ecological") %>%
        select(covType:rmse, sampling, polysize) %>%
        rename(rmseUtazi=rmse) %>%
        right_join(select(resultsDF, covType:rmse, model, sampling, converge, polysize)) %>%
        filter(converge == 0) %>%  #& model !="Riemann" & model != "Utazi") %>%
        mutate(improveRatio=(rmseUtazi-rmse)/rmseUtazi) %>%
        group_by(model) %>%
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
        coord_flip() +
        geom_hline(yintercept=0, linetype=2) +
        labs(x="Model", y="Relative Improvement") +
        ggtitle("RMSE: Margin of Improvement Over Utazi Model") +
        theme(panel.spacing.y = unit(0, "lines")) +
        guides(color=FALSE) +
        geom_text(nudge_y = .22))

(aggPlots$rmseRelativeZoom <- resultsDF %>%
    filter(model=="Ecological") %>%
    select(covType:rmse, sampling, polysize) %>%
    rename(rmseUtazi=rmse) %>%
    right_join(select(resultsDF, covType:rmse, model, sampling, converge, polysize)) %>%
    filter(converge == 0) %>%
    filter(model != "Ignore") %>%
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
    ggtitle("RMSE: Margin of Improvement Over Utazi Model") +
    theme(panel.spacing.y = unit(0, "lines")) +
    guides(color=FALSE) +
    geom_text(aes(y=upr), nudge_y = .1))

(aggPlots$bias <- resultsDF %>%
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
    guides(color=FALSE) +
    theme(
        strip.text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        title = element_text(size=25),
        axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=20)))

(aggPlots$biasPaper <- resultsDF %>%
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
        theme_bw() +
        facet_grid(rangeE~covType) +
        coord_flip() +
        geom_hline(yintercept=0, linetype=2) +
        labs(x="Model", y="Bias") +
        ggtitle("RMSE: Average Bias of Models") +
        theme(panel.spacing.y = unit(0, "lines")) +
        guides(color=FALSE) +
        theme(
            strip.text = element_text(size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            title = element_text(size=25),
            axis.text.x = element_text(size=12),
            axis.text.y = element_text(size=15),
            axis.title.x = element_text(size=20)))

ggsave(
    "demo/figures/biasSim1.png", aggPlots$biasPaper,
    width=400, height=280, units = "mm")

(aggPlots$dissDiff <- resultsDF %>%
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

(aggPlots$dissDiffPaper <- resultsDF %>%
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
    guides(color=FALSE) +
        theme(
            strip.text = element_text(size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            title = element_text(size=25),
            axis.text.x = element_text(size=12),
            axis.text.y = element_text(size=15),
            axis.title.x = element_text(size=20)))

ggsave(
    "demo/figures/dissSim1.png", aggPlots$dissDiffPaper,
    width=400, height=280, units = "mm")

pm <- resultsDF %>%
    group_by(covType, covVal, rangeE, M, seed, sampling, polysize) %>%
    mutate(groupC=all(converge == 0)) %>%
    filter(groupC) %>%
    ungroup %>%
    select(covType, covVal, rangeE, M, seed, sampling, polysize, 
           model, runtime) %>%
    spread("model", "runtime") %>%
    mutate(model=1:n()) %>%
    select(Ignore:Utazi) %>%
    ggpairs(aes(alpha=.1), lower=list(continuous="density")) +
    theme_classic()

pm2 <- pm
for(i in 2:pm$nrow) {
    for(j in 1:(i-1)) {
        pm2[i,j] <- pm[i,j] +
            scale_x_continuous(limits = c(0, 10)) +
            scale_y_continuous(limits = c(0, 10)) +
            geom_abline(color="red", linetype=2)
    }
}

pm2

write_rds(aggPlots, "~/Documents/PointPolygon/demo/aggplots.Rds")
