.libPaths(c("~/R3.6/", .libPaths()))
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
library(forcats)

rdsPathList <- list.files("~/Data/spaceTimeTest3/", full.names=TRUE)

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
    mutate(model=gsub("Reimann", "Riemann", model)) %>%
    mutate(model=gsub("IHME Resample", "Resample", model)) %>%
    mutate(model=gsub("Known", "Unmasked", model)) %>%
    mutate(model=gsub("Mixture Model", "Mixture", model)) %>%
    mutate(model=gsub("Utazi", "Ecological", model)) %>%
    mutate(model = fct_relevel(
        model,
        "Ignore", "Resample", "Ecological", "Riemann", "Mixture", "Unmasked"))

aggPlotsDR <- list()

(aggPlotsDR$coverage <- resultsDF %>%
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

(aggPlotsDR$coveragePaper <- resultsDF %>%
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
        theme_bw() +
        facet_grid(rangeE~covType) +
        coord_flip() +
        geom_hline(yintercept=.95, linetype=2) +
        labs(x="Model", y="") +
        ggtitle("95% Coverage of Province Probability(N=32)") +
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


(aggPlotsDR$rmseRelative <- resultsDF %>%
    filter(model=="Utazi") %>%
    select(covType:rmse) %>%
    rename(rmseUtazi=rmse) %>%
    right_join(select(resultsDF, covType:rmse, model, converge)) %>%
    filter(converge == 0 & model != "Ignore") %>%#rmse <.3) %>%
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
    geom_text(aes(y=upr), nudge_y = .04))

(aggPlotsDR$rmseRelativePaper <- resultsDF %>%
        filter(model=="Ecological") %>%
        select(covType:rmse) %>%
        rename(rmseUtazi=rmse) %>%
        right_join(select(resultsDF, covType:rmse, model, converge)) %>%
        filter(converge == 0 & model != "Riemann" & model != "Ignore") %>%
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
        labs(x="Model", y="Relative Improvement") +
        ggtitle("RMSE: Margin of Improvement Over Utazi Model") +
        theme(panel.spacing.y = unit(0, "lines")) +
        guides(color=FALSE) +
        geom_text(aes(y=upr), nudge_y = .04) +
        theme(
            strip.text = element_text(size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            title = element_text(size=25),
            axis.text.x = element_text(size=12),
            axis.text.y = element_text(size=15),
            axis.title.x = element_text(size=20)))

(aggPlotsDR$rmseProvRelative <- resultsDF %>%
        filter(model=="Ecological") %>%
        select(covType:seed, provrmse) %>%
        rename(rmseUtazi=provrmse) %>%
        right_join(select(resultsDF, covType:seed, provrmse, model, converge)) %>%
        filter(converge == 0 & model != "Ignore") %>%
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
        theme_bw() +
        facet_grid(rangeE~covType) +
        coord_flip() +
        geom_hline(yintercept=0, linetype=2) +
        labs(x="Model", y="Relative Improvement") +
        ggtitle("Province RMSE: Margin of Improvement Over Ecological Model") +
        theme(panel.spacing.y = unit(0, "lines")) +
        guides(color=FALSE) +
        geom_text(aes(y=upr), nudge_y = .06) +
        theme(
            strip.text = element_text(size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            title = element_text(size=25),
            axis.text.x = element_text(size=12),
            axis.text.y = element_text(size=15),
            axis.title.x = element_text(size=20)))

(aggPlotsDR$rmseProvRelativePaper <- resultsDF %>%
        filter(model=="Ecological") %>%
        select(covType:seed, provrmse) %>%
        rename(rmseUtazi=provrmse) %>%
        right_join(select(resultsDF, covType:seed, provrmse, model, converge)) %>%
        filter(converge == 0 & model != "Riemann" & model != "Ignore") %>%
        mutate(improveRatio=(rmseUtazi-provrmse)/rmseUtazi) %>%
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
        labs(x="Model", y="Relative Improvement") +
        ggtitle("Province RMSE: Margin of Improvement Over Ecological Model") +
        theme(panel.spacing.y = unit(0, "lines")) +
        guides(color=FALSE) +
        geom_text(aes(y=upr), nudge_y = .06) +
        theme(
            strip.text = element_text(size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            title = element_text(size=25),
            axis.text.x = element_text(size=12),
            axis.text.y = element_text(size=15),
            axis.title.x = element_text(size=20)))

(aggPlotsDR$rmseSingleProvRelativePaper <- resultsDF %>%
        filter(model=="Ecological") %>%
        select(covType:seed, provrmse) %>%
        rename(rmseUtazi=provrmse) %>%
        right_join(select(resultsDF, covType:seed, provrmse, model, converge)) %>%
        filter(converge == 0 & model != "Riemann") %>%
        mutate(improveRatio=(rmseUtazi-provrmse)/rmseUtazi) %>%
        group_by(model) %>%
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
        coord_flip() +
        geom_hline(yintercept=0, linetype=2) +
        labs(x="Model", y="Relative Improvement") +
        ggtitle("Province RMSE: Margin of Improvement Over IHME Resample Model") +
        theme(panel.spacing.y = unit(0, "lines")) +
        guides(color=FALSE) +
        geom_text(aes(y=upr), nudge_y = .06) +
        theme(
            strip.text = element_text(size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            title = element_text(size=25),
            axis.text.x = element_text(size=12),
            axis.text.y = element_text(size=15),
            axis.title.x = element_text(size=20)))

(aggPlotsDR$rmsePaper <- resultsDF %>%
        filter(converge == 0 & model != "Riemann" & rmse < .3) %>%
        group_by(covType, model, rangeE) %>%
        summarize(
            mu = mean(rmse),
            lwr = mean(rmse) - 1.96*(sd(rmse)/sqrt(n())),
            upr = mean(rmse) + 1.96*(sd(rmse)/sqrt(n()))) %>%
        ungroup %>%
        rename(Model=model) %>%
        mutate(txt=round(mu, 4)) %>%
        ggplot(aes(x=Model, ymin=lwr, y=mu, ymax=upr, color=Model, label=txt)) +
        geom_point() +
        geom_errorbar() +
        theme_bw() +
        facet_grid(rangeE~covType) +
        coord_flip() +
        geom_hline(yintercept=0, linetype=2) +
        labs(x="Model", y="RMSE") +
        ggtitle("Province RMSE") +
        theme(panel.spacing.y = unit(0, "lines")) +
        guides(color=FALSE) +
        geom_text(aes(y=upr), nudge_y = .002))

(aggPlotsDR$rmseProvPaper <- resultsDF %>%
        filter(converge == 0 & model != "Riemann" & rmse < .3) %>%
        group_by(covType, model, rangeE) %>%
        summarize(
            mu = mean(provrmse),
            lwr = mean(provrmse) - 1.96*(sd(provrmse)/sqrt(n())),
            upr = mean(provrmse) + 1.96*(sd(provrmse)/sqrt(n()))) %>%
        ungroup %>%
        rename(Model=model) %>%
        mutate(txt=round(mu, 4)) %>%
        ggplot(aes(x=Model, ymin=lwr, y=mu, ymax=upr, color=Model, label=txt)) +
        geom_point() +
        geom_errorbar() +
        theme_classic() +
        facet_grid(rangeE~covType) +
        coord_flip() +
        geom_hline(yintercept=0, linetype=2) +
        labs(x="Model", y="RMSE") +
        ggtitle("Province RMSE") +
        theme(panel.spacing.y = unit(0, "lines")) +
        guides(color=FALSE) +
        geom_text(aes(y=upr), nudge_y = .002))

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
        filter(converge == 0 & model != "Ignore") %>%
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
        mutate(Model=fct_rev(model)) %>%
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

(aggPlotsDR$runtime <- resultsDF %>%
    select(covType:seed, model, runtime) %>%
    filter(!(model %in% c("Unmaksed"))) %>%
    left_join(
        resultsDF %>%
            filter(model %in% c("Unmasked")) %>%
            select(covType:seed, runtime) %>%
            rename(IHME=runtime), 
        by=c("covType", "covVal", "rangeE", "seed")) %>%
    ggplot(aes(x=IHME, y=runtime)) +
    geom_point() +
    geom_abline() +
    theme_bw() +
    facet_wrap(~model, scales="free_y") +
    labs(x="Unmasked Runtime", y="Runtime") + 
    expand_limits(x = 0, y = 0) +
        theme(
            strip.text = element_text(size=15),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            title = element_text(size=25),
            axis.text.x = element_text(size=12),
            axis.text.y = element_text(size=15),
            axis.title.x = element_text(size=20)))

ggsave(
    "demo/figures/dissSim2.png", aggPlotsDR$dissDiffPaper,
    width=400, height=280, units = "mm")

ggsave(
    "demo/figures/biasSim2.png", aggPlotsDR$biasPaper,
    width=400, height=280, units = "mm")

ggsave(
    "demo/figures/covSim2.png", aggPlotsDR$coveragePaper,
    width=400, height=280, units = "mm")

ggsave(
    "demo/figures/provcovSim2.png", aggPlotsDR$provcoveragePaper,
    width=400, height=280, units = "mm")

ggsave(
    "demo/figures/rmseSim2.png", aggPlotsDR$rmseRelativePaper,
    width=400, height=280, units = "mm")

ggsave(
    "demo/figures/provrmseSim2.png", aggPlotsDR$rmseProvRelativePaper,
    width=400, height=280, units = "mm")

write_rds(aggPlotsDR, "~/Documents/PointPolygon/demo/aggplotsDR.Rds")
