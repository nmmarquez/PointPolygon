# TODO: 
# 1) Check intercept values across the models
# 2) check the differnence in the gp vs the RW2
# 3) Check possibility of constraints
# 4) check sum of all random effects RW and GP
# 5) Try resmpling more points
# 6) Map out the resampled points to pop surface

.libPaths(c("~/R3.6/", .libPaths()))
rm(list=ls())
library(boot)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(PointPolygon)
library(ggplot2)
library(sp)
library(Rcpp)
library(Matrix)
library(TMB)
library(sparseMVN)
library(SUMMER)
library(Biograph)
library(haven)
library(plotly)
library(forcats)
load("./demo/prepData.rda")

maxYear <- 2014
minYear <- 2000
nYears <- length(minYear:maxYear)
ageVec <- c(NN=0, PNN1=1, PNN2=2, `1yr`=3, `2yr`=4, `3yr`=5, `4yr`=6)

field <- simField(
    N = as.matrix(select(fullDF, long, lat)), rangeE = 1,
    offset = c(0.1, 0.2),
    max.edge = c(0.18,0.2),
    beta0 = -4,
    sigmaE = .3,
    nTimes = nYears,
    rho = .9, shape=spDF)

fullDF <- fullDF %>%
    rename(x=long, y=lat, oldid=id) %>%
    mutate(strat=paste(sprintf("%02d", reg), (2-urban), sep="_")) %>%
    left_join(
        select(filter(as_tibble(field$spdf), tidx==0), x, y, id),
        by = c("x", "y"))

syearWDF <- yearWDF %>%
    rename(oldid=id) %>%
    mutate(tidx=year-maxYear-1+nYears) %>%
    left_join(select(fullDF, oldid, id), by="oldid") %>%
    filter(year <= maxYear) %>%
    select(-year, -oldid) %>%
    arrange(id) %>%
    filter(tidx==max(tidx))

pointDF <- pointDF %>%
    filter((year >= minYear) & (year <= maxYear)) %>%
    mutate(denom=N, obs=died, yid=year-minYear, aid=ageVec[age_group]) %>%
    mutate(sid=as.numeric(as.factor(source))-1) %>%
    mutate(cid=group_indices(., sid, psu)-1) %>%
    rename(oldid=id) %>%
    left_join(select(fullDF, id, oldid), by="oldid") %>%
    select(denom, obs, id, yid, aid, sid, cid)

polyDF <- polyDF %>%
    filter((year >= minYear) & (year <= maxYear)) %>%
    mutate(denom=N, obs=died, yid=year-minYear, aid=ageVec[age_group]) %>%
    mutate(sid=max(pointDF$sid) + 1) %>%
    mutate(cid=group_indices(., psu) + max(pointDF$cid)) %>%
    mutate(polyid=as.numeric(as.factor(strat))-1) %>%
    select(denom, obs, id, yid, aid, sid, cid, polyid, strat, weight)

polyDF %>%
    group_by(yid, strat, aid) %>%
    summarize(
        mu = weighted.mean((obs/denom), weight)) %>%
    summarise(mu=1-prod(1-mu)) %>%
    left_join(
        syearWDF %>%
            group_by(strat) %>%
            summarize(Population = sum(Population))) %>%
    summarise(mu=weighted.mean(mu, Population)) %>%
    mutate(Year=yid+2000) %>%
    ggplot(aes(x=Year, y=mu)) +
    geom_line() +
    theme_classic()

field$spdf <- field$spdf %>%
    left_join(select(fullDF, id, urban), by="id")

oosDF <- bind_rows(lapply(0:14, function(i){
    modRes <- readRDS(sprintf("~/Data/dataRun/model_y%s_r%s_.RDS", i, "NA"))
    print(i)
    
    bind_rows(
        bind_rows(lapply(names(modRes$preds$noRE), function(n){
            pointDF %>%
                filter(yid == i) %>%
                rename(tidx = yid) %>%
                left_join(modRes$preds$noRE[[n]]) %>%
                mutate(nll=-dbinom(obs, denom, mu, log=T)) %>%
                summarize(nll=sum(nll)) %>%
                mutate(pars="noRE", method=n)
        })),
        
        bind_rows(lapply(names(modRes$preds$temporal), function(n){
            pointDF %>%
                filter(yid == i) %>%
                rename(tidx = yid) %>%
                left_join(modRes$preds$temporal[[n]]) %>%
                mutate(nll=-dbinom(obs, denom, mu, log=T)) %>%
                summarize(nll=sum(nll)) %>%
                mutate(pars="temporal", method=n)
        })),
        
        bind_rows(lapply(names(modRes$preds$full), function(n){
            pointDF %>%
                filter(yid == i) %>%
                rename(tidx = yid) %>%
                left_join(modRes$preds$full[[n]]) %>%
                mutate(nll=-dbinom(obs, denom, mu, log=T)) %>%
                summarize(nll=sum(nll)) %>%
                mutate(pars="full", method=n)
        }))) %>%
        mutate(yid=i)
}))

results <- list()

results$bestmodel <- oosDF %>%
    group_by(yid) %>%
    mutate(ismin=min(nll) == nll) %>%
    filter(ismin & yid != 14) %>%
    group_by(method, pars) %>%
    summarize(N=n()) %>%
    arrange(-N)

oosRegDF <- bind_rows(lapply(0:9, function(i){
    modRes <- readRDS(sprintf("~/Data/dataRun/model_y%s_r%s_.RDS", "NA", i))
    print(i)
    
    pointDF_ <- pointDF %>%
        left_join(fullDF) %>%
        mutate(polyid=reg-1)
    
    bind_rows(
        bind_rows(lapply(names(modRes$preds$noRE), function(n){
            pointDF_ %>%
                filter(polyid == i) %>%
                rename(tidx = yid) %>%
                left_join(modRes$preds$noRE[[n]]) %>%
                mutate(nll=-dbinom(obs, denom, mu, log=T)) %>%
                summarize(nll=sum(nll)) %>%
                mutate(pars="noRE", method=n)
        })),
        
        bind_rows(lapply(names(modRes$preds$temporal), function(n){
            pointDF_ %>%
                filter(polyid == i) %>%
                rename(tidx = yid) %>%
                left_join(modRes$preds$temporal[[n]]) %>%
                mutate(nll=-dbinom(obs, denom, mu, log=T)) %>%
                summarize(nll=sum(nll)) %>%
                mutate(pars="temporal", method=n)
        })),
        
        bind_rows(lapply(names(modRes$preds$full), function(n){
            pointDF_ %>%
                filter(polyid == i) %>%
                rename(tidx = yid) %>%
                left_join(modRes$preds$full[[n]]) %>%
                mutate(nll=-dbinom(obs, denom, mu, log=T)) %>%
                summarize(nll=sum(nll)) %>%
                mutate(pars="full", method=n)
        }))) %>%
        mutate(polyid=i)
}))

results$bestmodelreg <- oosRegDF %>%
    group_by(polyid) %>%
    mutate(ismin=min(nll) == nll) %>%
    filter(ismin) %>%
    group_by(method, pars) %>%
    summarize(N=n()) %>%
    arrange(-N)

# Best model is mixture with temporal effects
modelRes <- readRDS("~/Data/dataRun/model_yNA_rNA_.RDS")
aggShapes <- readRDS("./demo/aggShapes.Rds")

results$predDF5q0 <- simulateDataFieldCI(
    modelRes$modelList$temporal$mixture, field, agg5q0 = T)

results$predDF5q0Resample <- simulateDataFieldCI(
    modelRes$modelList$temporal$resample, field, agg5q0 = T)

results$predDF5q0Point <- simulateDataFieldCI(
    modelRes$modelList$temporal$point, field, agg5q0 = T)

results$Map <- results$predDF5q0 %>%
    {left_join(field$spdf, .)} %>%
    ggplot(aes(x, y, fill = mu)) +
    geom_raster() +
    coord_equal() +
    theme_void() +
    scale_fill_distiller(palette = "Spectral") +
    facet_wrap(~tidx)

results$MapResample <- results$predDF5q0Resample %>%
    {left_join(field$spdf, .)} %>%
    ggplot(aes(x, y, fill = mu)) +
    geom_raster() +
    coord_equal() +
    theme_void() +
    scale_fill_distiller(palette = "Spectral") +
    facet_wrap(~tidx)

results$MapPoint <- results$predDF5q0Point %>%
    {left_join(field$spdf, .)} %>%
    ggplot(aes(x, y, fill = mu)) +
    geom_raster() +
    coord_equal() +
    theme_void() +
    scale_fill_distiller(palette = "Spectral") +
    facet_wrap(~tidx)

results$provPredDF5q0 <- arealDataCI(
    modelRes$modelList$temporal$mixture, 
    field,
    agg5q0 = TRUE,
    polygonList = lapply(1:nrow(aggShapes$provShape@data), function(i){
        sf::st_as_sf(aggShapes$provShape[i,])}),
    popDF = mutate(select(syearWDF, -tidx), w=Population),
    draws=1000)

results$provPredDF5q0Resample <- arealDataCI(
    modelRes$modelList$temporal$resample, 
    field,
    agg5q0 = TRUE,
    polygonList = lapply(1:nrow(aggShapes$provShape@data), function(i){
        sf::st_as_sf(aggShapes$provShape[i,])}),
    popDF = mutate(select(syearWDF, -tidx), w=Population),
    draws=1000)

(results$provmap <- st_as_sf(aggShapes$provShape) %>%
    mutate(polyid=as.numeric(Prov)-1) %>%
    right_join(results$provPredDF5q0) %>%
    mutate(Year = tidx + 2000) %>%
    mutate(mu = mu * 1000) %>%
    ggplot() +
    geom_sf(aes(fill=mu)) +
    theme_void() +
    scale_fill_distiller(palette = "Spectral") +
    facet_wrap(~Year) +
    labs(fill="U5MR Per\n1000 Births"))

results$provMapResample <- st_as_sf(aggShapes$provShape) %>%
    mutate(polyid=as.numeric(Prov)-1) %>%
    right_join(results$provPredDF5q0Resample) %>%
    ggplot() +
    geom_sf(aes(fill=mu)) +
    theme_void() +
    scale_fill_distiller(palette = "Spectral") +
    facet_wrap(~tidx)

results$provNatDF5q0 <- arealDataCI(
    modelRes$modelList$temporal$mixture, 
    field,
    agg5q0 = TRUE,
    polygonList = list(sf::st_as_sf(field$bound)),
    popDF = mutate(select(syearWDF, -tidx), w=Population),
    draws=1000)

results$provNatDF5q0Resample <- arealDataCI(
    modelRes$modelList$temporal$resample, 
    field,
    agg5q0 = TRUE,
    polygonList = list(sf::st_as_sf(field$bound)),
    popDF = mutate(select(syearWDF, -tidx), w=Population),
    draws=1000)

results$provNatDF5q0Point <- arealDataCI(
    modelRes$modelList$temporal$point, 
    field,
    agg5q0 = TRUE,
    polygonList = list(sf::st_as_sf(field$bound)),
    popDF = mutate(select(syearWDF, -tidx), w=Population),
    draws=1000)

ihmeDF <- read_csv("demo/5q0Results_estimates.csv") %>%
    filter(age_group_id == 1 & year >= 2000 & year <= 2014) %>%
    rename(mu=mean, upr=upper, lwr=lower, Year = year) %>%
    mutate(method = "IHME GBD Estimates", tidx=Year-2000) %>%
    select(mu, upr, lwr, Year, method, tidx)

u5mHTDF <- "~/Documents/PointPolygon/demo/DRBR61DT/DRBR61FL.DTA" %>%
    getBirths(surveyyear = 2013, year.cut = 2000:2014) %>%
    countrySummary(
        years = unique(.$time), regionVar = "strata", timeVar = "time",
        clusterVar = "~v001+v002", national.only = TRUE) %>%
    filter(region == "All") %>%
    arrange(years) %>%
    mutate(Year=2000:2012, tidx = Year-2000, method = "HT DHS 2013") %>%
    rename(mu=u5m, lwr=lower, upr=upper) %>%
    select(mu, upr, lwr, Year, method, tidx)

micsDF <- paste0(
    "~/Documents/PointPolygon/demo/Dominican Republic_MICS5_Datasets/",
    "Dominican Republic MICS 2014 SPSS Datasets/bh.sav") %>%
    read_sav() %>%
    filter(BH5 != 9) %>%
    rename(PSU = HH1, urban = HH6, region = HH7) %>%
    filter(BH4Y != 9999 & BH4Y != 9997) %>%
    filter(BH4M != 99 & BH4M != 97) %>%
    mutate(dob=sprintf("%i-%02d-01", BH4Y, BH4M)) %>%
    mutate(dob = Date_as_cmc(dob, "%Y-%m-%d")$cmc) %>%
    mutate(alive = ifelse(BH5 == 1, "yes", "no")) %>%
    mutate(`DeathMonth` = floor(c(1/30.5, 1, 12)[BH9U] * BH9N)) %>%
    filter(BH5 == 1 | !is.na(DeathMonth)) %>%
    as.data.frame() %>%
    {getBirths(
        variables = c("PSU", "wmweight"),
        data = ., surveyyear = 2014, date.interview = "WDOI",
        dob = "dob", alive = "alive", age = "DeathMonth",
        strata = c("urban", "region"),
        year.cut = 2000:2014
    )} %>%
    countrySummary(
        years = unique(.$time), clusterVar = "~strata+PSU", 
        weightsVar = "wmweight", national.only=TRUE) %>%
    filter(region == "All") %>%
    arrange(years) %>%
    mutate(Year=2000:2013, tidx = Year-2000, method = "HT MICS 2014") %>%
    rename(mu=u5m, lwr=lower, upr=upper) %>%
    select(mu, upr, lwr, Year, method, tidx)

priors <- simhyper(
    R = 2, nsamp = 1e+05, nsamp.check = 5000, Amat = mat, only.iid = TRUE)

micsSmoothDF <- paste0(
    "~/Documents/PointPolygon/demo/Dominican Republic_MICS5_Datasets/",
    "Dominican Republic MICS 2014 SPSS Datasets/bh.sav") %>%
    read_sav() %>%
    filter(BH5 != 9) %>%
    rename(PSU = HH1, urban = HH6, region = HH7) %>%
    filter(BH4Y != 9999 & BH4Y != 9997) %>%
    filter(BH4M != 99 & BH4M != 97) %>%
    mutate(dob=sprintf("%i-%02d-01", BH4Y, BH4M)) %>%
    mutate(dob = Date_as_cmc(dob, "%Y-%m-%d")$cmc) %>%
    mutate(alive = ifelse(BH5 == 1, "yes", "no")) %>%
    mutate(`DeathMonth` = floor(c(1/30.5, 1, 12)[BH9U] * BH9N)) %>%
    filter(BH5 == 1 | !is.na(DeathMonth)) %>%
    as.data.frame() %>%
    {getBirths(
        variables = c("PSU", "wmweight"),
        data = ., surveyyear = 2014, date.interview = "WDOI",
        dob = "dob", alive = "alive", age = "DeathMonth",
        strata = c("urban", "region"),
        year.cut = 2000:2014
    )} %>%
    countrySummary(
        years = unique(.$time), clusterVar = "~strata+PSU", 
        weightsVar = "wmweight", national.only=TRUE) %>%
    filter(years != "100-100") %>%
    arrange(years) %>%
    fitINLA(
        geo = NULL, Amat = NULL, 
        year_names = .$years, year_range = c(2001, 2013),
        priors = priors, rw = 2,
        is.yearly=TRUE, m = 1, type.st = 0) %>%
    .[["fit"]] %>%
    summary() %>%
    .[["linear.predictor"]] %>%
    as_tibble() %>%
    head(n=13) %>%
    mutate(mu = inv.logit(mean)) %>%
    mutate(lwr = inv.logit(`0.025quant`), upr = inv.logit(`0.975quant`)) %>%
    select(mu, lwr, upr) %>%
    mutate(method = "MICS 2014 Smoothed", Year = 2001:2013)

u5mSmoothHTDF <- "~/Documents/PointPolygon/demo/DRBR61DT/DRBR61FL.DTA" %>%
    getBirths(surveyyear = 2013, year.cut = 2000:2014) %>%
    countrySummary(
        years = unique(.$time), regionVar = "strata", timeVar = "time",
        clusterVar = "~v001+v002", national.only = TRUE) %>%
    filter(region == "All" & years != "100-100") %>%
    arrange(years) %>%
    fitINLA(
        geo = NULL, Amat = NULL, 
        year_names = .$years, year_range = c(2001, 2012),
        priors = priors, rw = 2,
        is.yearly=TRUE, m = 1, type.st = 0) %>%
    .[["fit"]] %>%
    summary() %>%
    .[["linear.predictor"]] %>%
    as_tibble() %>%
    head(n=12) %>%
    mutate(mu = inv.logit(mean)) %>%
    mutate(lwr = inv.logit(`0.025quant`), upr = inv.logit(`0.975quant`)) %>%
    select(mu, lwr, upr) %>%
    mutate(method = "DHS 2013 Smoothed", Year = 2001:2012)

(results$compare5q0 <- bind_rows(
    mutate(results$provNatDF5q0, method="Mixture"),
    mutate(results$provNatDF5q0Resample, method="Resample"),
    u5mHTDF, micsDF) %>%
    mutate(Year=tidx+2000) %>%
    mutate(mu = mu * 1000, lwr = lwr * 1000, upr = upr * 1000) %>%
    ggplot(aes(x=Year, y=mu, ymin=lwr, ymax=upr, group=method)) +
    geom_line(aes(color=method)) +
    geom_ribbon(aes(fill=method), alpha=.3) +
    theme_classic() +
    labs(y="U5MR Per 1000 Births", fill="Method", color="Method") +
    ggtitle("Change In Dominican Republic U5MR: Direct Estimate Compare") +
        theme(
            legend.text = element_text(size=13),
            legend.title = element_text(size=15),
            axis.text = element_text(size=13),
            axis.title = element_text(size=17),
            title =  element_text(size=20)
        ))

results$compareIHME5q0 <- bind_rows(
    mutate(results$provNatDF5q0, method="Mixture"),
    mutate(results$provNatDF5q0Resample, method="Resample"),
    ihmeDF) %>%
    mutate(Year=tidx+2000) %>%
    mutate(mu = mu * 1000, lwr = lwr * 1000, upr = upr * 1000) %>%
    ggplot(aes(x=Year, y=mu, ymin=lwr, ymax=upr, group=method)) +
    geom_line(aes(color=method)) +
    geom_ribbon(aes(fill=method), alpha=.3) +
    theme_classic() +
    labs(y="U5MR Per 1000 Births", fill="Method", color="Method") +
    ggtitle("Change In Dominican Republic U5MR") +
    theme(
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        axis.text = element_text(size=13),
        axis.title = element_text(size=17),
        title =  element_text(size=20))

results$compareSmooth5q0 <- bind_rows(
    mutate(results$provNatDF5q0, method="Mixture"),
    mutate(results$provNatDF5q0Resample, method="Resample"),
    u5mSmoothHTDF %>% mutate(tidx=Year-2000), 
    micsSmoothDF %>% mutate(tidx=Year-2000)) %>%
    mutate(Year=tidx+2000) %>%
    mutate(mu = mu * 1000, lwr = lwr * 1000, upr = upr * 1000) %>%
    ggplot(aes(x=Year, y=mu, ymin=lwr, ymax=upr, group=method)) +
    geom_line(aes(color=method)) +
    geom_ribbon(aes(fill=method), alpha=.3) +
    theme_classic() +
    labs(y="U5MR Per 1000 Births", fill="Method", color="Method") +
    ggtitle("Change In Dominican Republic U5MR: Smooth Direct Compare") +
    theme(
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        axis.text = element_text(size=13),
        axis.title = element_text(size=17),
        title =  element_text(size=20))

ggsave(
    "demo/figures/rmseIHMECompare.png", results$compareIHME5q0,
    width=400, height=280, units = "mm")

ggsave(
    "demo/figures/rmseSmoothCompare.png", results$compareSmooth5q0,
    width=400, height=280, units = "mm")

ggsave(
    "demo/figures/rmseDirectCompare.png", results$compare5q0,
    width=400, height=280, units = "mm")

bind_rows(
    mutate(results$provNatDF5q0, method="Mixture"),
    mutate(results$provNatDF5q0Resample, method="Resample"),
    ihmeDF, u5mHTDF, micsDF) %>%
    mutate(Year=tidx+2000) %>%
    mutate(mu = mu * 1000, lwr = lwr * 1000, upr = upr * 1000) %>%
    filter(
        method == "IHME GBD Estimates" | 
            method == "Resample" |
            method == "Mixture") %>% 
    ggplot(aes(x=Year, y=mu, ymin=lwr, ymax=upr, group=method)) +
    geom_line(aes(color=method)) +
    geom_ribbon(aes(fill=method), alpha=.3) +
    theme_classic() +
    labs(y="U5MR Per 1000 Births", fill="Method", color="Method") +
    ggtitle("Change In Dominican Republic U5MR") +
    theme(
        legend.text = element_text(size=13),
        legend.title = element_text(size=15),
        axis.text = element_text(size=13),
        axis.title = element_text(size=17),
        title =  element_text(size=20)) +
    scale_fill_manual(values=c("#000000", "#b7a57a", "#4b2e83")) +
    scale_color_manual(values=c("#000000", "#b7a57a", "#4b2e83"))

saveRDS(results, "~/Documents/PointPolygon/demo/resultsData.RDS")

# Extra Stuff
extractAge <- function(model){
    ages <- c("0-28 Days", "1-6 Months", "6-12 Months", paste0(1:4, " Years"))
    ages <- factor(ages, levels = ages)
    idx <- c(1, which(names(model$sd$par.fixed) == "beta_age"))
    pars_ <- model$sd$par.fixed
    pars_ <- c(pars_[1], pars_[1] + pars_[2:length(pars_)])
    sds_ <- sqrt(diag(model$sd$cov.fixed))
    for(i in 2:length(sds_)){
        sds_[i] <- sqrt(sds_[i]^2 + sds_[1]^2 + 2*model$sd$cov.fixed[1,i])
    }
    tibble(beta=pars_[idx], std.err=sds_[idx], age=ages)
}

extractWalks <- function(model){
    ages <- c("0-28 Days", "1-6 Months", "6-12 Months", paste0(1:4, " Years"))
    ages <- factor(ages, levels = ages)
    idx <- which(row.names(model$sd$jointPrecision) == "phi")
    subMat <- model$sd$jointPrecision[idx, idx]
    sds <- sqrt(diag(solve(subMat)))
    mu <- model$sd$par.random[names(model$sd$par.random) == "phi"]
    tibble(mu=mu, std.err=sds, age=rep(ages, 15), Year=rep(2000:2014, each=7))
}

bind_rows(
    extractAge(modelRes$modelList$temporal$mixture) %>%
        mutate(Model="Mixture"),
    extractAge(modelRes$modelList$temporal$resample) %>%
        mutate(Model="Resample")) %>%
    mutate(lwr=beta - 1.96*std.err) %>%
    mutate(upr=beta + 1.96*std.err) %>%
    ggplot(aes(x=Model, y=beta, ymin=lwr, ymax=upr, color=Model)) +
    geom_point() +
    geom_errorbar() +
    facet_wrap(~age, ncol = 1) +
    theme_classic() +
    coord_flip()

bind_rows(
    extractWalks(modelRes$modelList$temporal$mixture) %>%
        mutate(Model="Mixture"),
    extractWalks(modelRes$modelList$temporal$resample) %>%
        mutate(Model="Resample")) %>%
    mutate(lwr=mu - 1.96*std.err) %>%
    mutate(upr=mu + 1.96*std.err) %>%
    ggplot(aes(x=Year, y=mu, ymin=lwr, ymax=upr, group=Model)) +
    geom_line(aes(color=Model)) +
    geom_ribbon(aes(fill=Model), alpha=.3) +
    facet_wrap(~age) +
    theme_classic()

bind_rows(
    extractWalks(modelRes$modelList$temporal$mixture) %>%
        mutate(Model="Mixture"),
    extractWalks(modelRes$modelList$temporal$resample) %>%
        mutate(Model="Resample")) %>%
    mutate(lwr=mu - 1.96*std.err) %>%
    mutate(upr=mu + 1.96*std.err) %>%
    filter(age != "4 Years") %>%
    ggplot(aes(x=Year, y=mu, ymin=lwr, ymax=upr, group=Model)) +
    geom_line(aes(color=Model)) +
    geom_ribbon(aes(fill=Model), alpha=.3) +
    facet_wrap(~age) +
    theme_classic()

