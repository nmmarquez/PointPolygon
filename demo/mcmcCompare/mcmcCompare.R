# Compare mcmc estimation of mixture to LA 
.libPaths(c("~/R3.6/", .libPaths()))
rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(PointPolygon)
library(sp)
library(Rcpp)
library(Matrix)
library(TMB)
library(sparseMVN)
library(tidybayes)
library(ggplot2)

rangeE <- .7 # range of spatial prces varies from {.3, .5, .7}
covVal <- -2 # covariate effect in set {.2, .4, -.5, .2, -2}
covType <- "spatial" # either random spatial or cluster 
M <- 50 # number of samples in Polygons chosen from U(50, 300)
seed <- 123 # RNG
nT <- 1 # time periods

modelname <- paste0(
    "range=", rangeE,
    ",cov=", covVal,
    ",covtype=", covType,
    ",M=", M,
    ",seed=", seed, ".Rds"
)

set.seed(seed)

save_file <- "~/Documents/PointPolygon/demo/mcmcCompare/results.RDS"

if(!file.exists(save_file)){
    unitSim <- simField(
        N = 60, # use 60 squares as base field just as in example
        sigmaE = 1, # unit variance for spde process
        rangeE = rangeE,
        shape = NULL, # null shape which by default creates unit square
        beta0 = -2, # intercept
        betaList = list(list(type=covType, value=covVal)), 
        link = arm::invlogit,
        offset = c(0.1, 0.2),
        max.edge = c(0.05,0.2))
    
    rWidthSamples <- c("10"=10)
    pointSamples <- 250
    totalTrials <- pointSamples*M
    
    # sample with mix of nonoverlapping
    mixDFList <- lapply(rWidthSamples, function(x){
        samplePerPolygon <- totalTrials / x^2
        sampleTrailClusterRatio <- pointSamples / M
        clustersPerPolygon <- round(sqrt(samplePerPolygon/sampleTrailClusterRatio))
        trialPerCluster <- round(sqrt(samplePerPolygon*sampleTrailClusterRatio))
        samplePPMix(
            unitSim,
            N = clustersPerPolygon,
            M = trialPerCluster,
            p = .5, # sample all polygons
            rWidth = x)})
    
    
    modelStandard <- runFieldModel(
        unitSim,
        pointDF=mixDFList$`10`$pointDF,
        polyDF=mixDFList$`10`$polyDF,
        verbose=T,
        control=list(),
        moption=0,
        rWidth=10
    )
    
    fes <- modelStandard$opt$par
    
    modelMCMC <- runFieldModel(
        unitSim,
        pointDF=mixDFList$`10`$pointDF,4 9m5h
        polyDF=mixDFList$`10`$polyDF,
        verbose=T,
        mcmc = TRUE,
        chains = 4,
        init = "par",
        iter = 2000,
        warmup = 1000,
        cores = 4,
        start = list(
            z=matrix(unname(modelStandard$sd$par.random), ncol=nT),
            beta=fes[names(fes) == "beta"],
            log_tau = fes[names(fes) == "log_tau"],
            log_kappa = fes[names(fes) == "log_kappa"]))
    
    saveRDS(
        list(
            laModel = modelStandard,
            mcModel = modelMCMC,
            field   = unitSim,
            data    = mixDFList
        ),
        file = save_file)
}

results <- readRDS(save_file)

# fixed effects plot
results$mcModel$opt %>%
    spread_draws(beta[1:2], log_tau, log_kappa) %>%
    spread(`1:2`, beta) %>%
    rename(beta1=`1`, beta2=`2`) %>%
    gather("parameter", "value", log_tau:beta2) %>%
    mutate(.chain = as.factor(.chain)) %>%
    ggplot(aes(x=.iteration, y=value, group=.chain, color=.chain)) +
    geom_line() +
    theme_classic() +
    geom_hline(
        aes(yintercept=value),
        linetype=3,
        data = data.frame(
            value = results$laModel$opt$par,
            parameter = c("beta1", "beta2", "log_tau", "log_kappa"))) +
    facet_wrap(~parameter)

results$mcModel$opt %>%
    spread_draws(beta[1:2], log_tau, log_kappa) %>%
    filter(.chain !=4) %>%
    spread(`1:2`, beta) %>%
    rename(beta1=`1`, beta2=`2`) %>%
    gather("parameter", "value", log_tau:beta2) %>%
    mutate(.chain = as.factor(.chain)) %>%
    ggplot(aes(x=.iteration, y=value, group=.chain, color=.chain)) +
    geom_line() +
    theme_classic() +
    geom_hline(
        aes(yintercept=value),
        linetype=3,
        data = data.frame(
            value = results$laModel$opt$par,
            parameter = c("beta1", "beta2", "log_tau", "log_kappa"))) +
    facet_wrap(~parameter)

maxZ <- 16
results$mcModel$opt %>%
    spread_draws(z[1:1617]) %>%
    rename(zindex=`1:1617`) %>%
    filter(zindex <= maxZ) %>%
    mutate(.chain = as.factor(.chain)) %>%
    ggplot(aes(x=.iteration, y=z, group=.chain, color=.chain)) +
    geom_line(alpha=.5) +
    geom_hline(
        aes(yintercept=z),
        linetype=3,
        data = data.frame(
            z = results$laModel$sd$par.random[1:maxZ],
            zindex = 1:maxZ)) +
    theme_classic() +
    facet_wrap(~zindex)

compareDF <- results$mcModel$opt %>%
    spread_draws(z[1:1617]) %>%
    rename(zindex=`1:1617`) %>%
    filter(.chain==1) %>%
    group_by(zindex) %>%
    summarize(
        mu=mean(z),
        lwr=quantile(z, probs=.025),
        upr=quantile(z, probs=.975)
    ) %>%
    left_join(
        data.frame(zLA = results$laModel$sd$par.random, zindex = 1:1617),
        by = "zindex")

compareDF %>%
    ggplot(aes(x=zLA, y=mu, ymin=lwr, ymax=upr)) +
    geom_point(alpha=.5) + 
    geom_errorbar(alpha=.5) +
    theme_classic() +
    geom_abline(color="red", linetype=3)


laFit <- simulateFieldCI(results$field, results$laModel)
mcFit <- lapply(1:4, function(i){
    simulateFieldCI(results$field, results$mcModel, chain=i)
})
names(mcFit) <- paste0("MCMC ", 1:4)

ggFieldEst(results$field, c(list(TMB=laFit), mcFit)) + 
    facet_wrap(~Type)

ciCompare <- ggFieldEst(results$field, c(list(TMB=laFit), mcFit), sd=TRUE) + 
    facet_wrap(~Type) +
    labs(fill = "Confidence/\nCredible\nInterval", x = "", y = "") +
    theme(
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        axis.text = element_text(size=13),
        axis.title = element_text(size=17),
        title =  element_text(size=20),
        strip.text.x = element_text(size = 17))
    

resCompare <- mcFit$`MCMC 1` %>%
    rename_all(function(x) paste0("MCMC", x)) %>%
    rename(id = MCMCid) %>%
    left_join(laFit) %>%
    ggplot(aes(x=mu, y=MCMCmu, ymin = MCMClwr, ymax = MCMCupr)) +
    geom_point(alpha=.3) +
    geom_errorbar(alpha=.1, size=.2) +
    theme_classic() +
    labs(x="Laplace Approximation", y="MCMC") +
    ggtitle("") +
    theme(
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        axis.text = element_text(size=13),
        axis.title = element_text(size=17),
        title =  element_text(size=20))

ggsave(
    "demo/figures/MCMCvINLAsd.png", ciCompare,
    width=400, height=280, units = "mm")

ggsave(
    "demo/figures/MCMCvINLA.png", resCompare,
    width=400, height=280, units = "mm")
