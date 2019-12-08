.libPaths(c("~/R3.5/", .libPaths()))
rm(list=ls())
library(arm)
library(rgeos)
library(sp)
require(tidyr)
require(gridExtra)
require(dplyr)
require(ggplot2)
library(PointPolygon)
library(stringr)
library(sf)
library(Rcpp)
library(rgdal)

# range of the underlying spatial process as defined 
rangeE <- .5 # range of spatial prces varies from {.3, .5, .7}
covVal <- 2 # covariate effect in set {.2, .4, -.5, .2, -2}
M <- 100 # number of samples in Polygons chosen from U(50, 300)
seed <- 12345 # RNG

set.seed(seed)
unitSim <- simField(
    N = 60, # how detailed the grid is creates an NxN grid
    offset = c(0.1, 0.2), # all simulations will have this to create pred mesh 
    max.edge = c(0.05,0.2), # all simulations will have this to create pred mesh
    beta0 = -2, # the intercept term
    rangeE = rangeE,
    betaList = list(list(type="random", value=covVal))) # the cov type & value


# just some R plotting 
plotList <- lapply(c("V1", "z", "theta"), function(eff){
    unitSim$spdf %>%
        gather("Effect", "Value", V0:theta) %>%
        filter(Effect==eff) %>%
        ggplot(aes(x, y, fill=Value)) +
        geom_raster() +
        coord_equal() +
        theme_void() +
        scale_fill_distiller(palette = "Spectral") +
        ggtitle(NULL) +
        guides(fill=FALSE)
})


unitSim2 <- simField(
    N = 60, # how detailed the grid is creates an NxN grid
    offset = c(0.1, 0.2), # all simulations will have this to create pred mesh 
    max.edge = c(0.05,0.2), # all simulations will have this to create pred mesh
    beta0 = -2, # the intercept term
    rangeE = rangeE,
    betaList = list(list(type="spatial", value=covVal))) # the cov type & value


# just some R plotting 
plotList2 <- lapply(c("V1", "z", "theta"), function(eff){
    unitSim2$spdf %>%
        gather("Effect", "Value", V0:theta) %>%
        filter(Effect==eff) %>%
        ggplot(aes(x, y, fill=Value)) +
        geom_raster() +
        coord_equal() +
        theme_void() +
        scale_fill_distiller(palette = "Spectral") +
        ggtitle(NULL) +
        guides(fill=FALSE)
})



unitSim3 <- simField(
    N = 60, # how detailed the grid is creates an NxN grid
    offset = c(0.1, 0.2), # all simulations will have this to create pred mesh 
    max.edge = c(0.05,0.2), # all simulations will have this to create pred mesh
    beta0 = -2, # the intercept term
    rangeE = rangeE,
    betaList = list(list(type="cluster", value=covVal))) # the cov type & value


plotList3 <- lapply(c("V1", "z", "theta"), function(eff){
    unitSim3$spdf %>%
        gather("Effect", "Value", V0:theta) %>%
        filter(Effect==eff) %>%
        ggplot(aes(x, y, fill=Value)) +
        geom_raster() +
        coord_equal() +
        theme_void() +
        scale_fill_distiller(palette = "Spectral") +
        ggtitle(NULL) +
        guides(fill=FALSE)
})

unitSim <- simField(
    N = 100, # use 60 squares as base field just as in example
    sigmaE = 1, # unit variance for spde process
    rangeE = .3,
    shape = NULL, # null shape which by default creates unit square
    beta0 = -2, # intercept
    betaList = list(list(type="spatial", value=1.5)), 
    link = arm::invlogit,
    offset = c(0.1, 0.2),
    max.edge = c(0.05,0.2))

mixList <- samplePPMix(unitSim, N = 5, M = 20, rWidth=5)
polyList <- samplePolygons(unitSim, N = 5, M = 20, rWidth=5)

samplePlots <- list()

(samplePlots$fieldPlot <- ggField(unitSim) +
    labs(fill="Probability", title="Continuous Probability of Event Field") + 
    theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
    ))

(samplePlots$mixPlot <- mixList$polyDF %>%
    select(-id) %>%
    rename(id=trueid) %>%
    select(tidx, trials, obs, polyid, id) %>%
    mutate(locKnown=FALSE) %>%
    rbind(mutate(mixList$pointDF, locKnown=TRUE)) %>%
    select(id, locKnown) %>%
    right_join(as_tibble(unitSim$spdf)) %>%
    ggplot(aes(x=x, y=y, fill=locKnown)) +
    geom_raster() +
    theme_classic() +
    scale_fill_discrete(na.value="white") +
    geom_vline(xintercept=seq(0, 1, .2)) +
    geom_hline(yintercept=seq(0, 1, .2)) +
    coord_equal() +
    scale_x_continuous(limits=c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits=c(0, 1), expand = c(0, 0)) +
    theme_void() +
    labs(fill="Geolocated"))

(samplePlots$ovPlot <- polyList %>%
    select(id=trueid) %>%
    mutate(locKnown=as.logical(rbinom(n(), 1, .5))) %>%
    right_join(as_tibble(unitSim$spdf)) %>%
    ggplot(aes(x=x, y=y, fill=locKnown)) +
    geom_raster() +
    theme_classic() +
    scale_fill_discrete(na.value="white") +
    geom_vline(xintercept=seq(0, 1, .2)) +
    geom_hline(yintercept=seq(0, 1, .2)) +
    coord_equal() +
    scale_x_continuous(limits=c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits=c(0, 1), expand = c(0, 0)) +
    theme_void() +
    labs(fill="Geolocated"))

sampUnitResults <- 
    "demo/Results/range=0.3,cov=-0.5,covtype=random,M=300,seed=6.Rds"

rezUnit <- readRDS(sampUnitResults)
rezUnit$sim$spdf <- rezUnit$sim$spdf %>%
    st_as_sf() %>%
    mutate(tidx=0)

predList <- lapply(list(
    `Mixture Model` = rezUnit$pred$`Mixed Model`$overlap$`5`,
    Riemann = rezUnit$pred$Reimann$overlap$`5`,
    Utazi = rezUnit$pred$Utazi$overlap$`5`,
    Resample = rezUnit$pred$Resample$overlap$`5`,
    Known = rezUnit$pred$Known$overlap$`5`), function(df){
        df %>% mutate(tidx=0)
    })

(samplePlots$results <- ggFieldEst(rezUnit$sim, predList) +
    labs(fill="Probability") +
    facet_wrap(~Type))

(samplePlots$resultsSD <- ggFieldEst(rezUnit$sim, predList, sd=T) +
    labs(fill="Std. Error") +
    facet_wrap(~Type))

(samplePlots$resultsPaper <- ggFieldEst(rezUnit$sim, predList[-2]) +
        labs(fill="Probability") +
        facet_wrap(~Type))

(samplePlots$resultsSDPaper <- ggFieldEst(rezUnit$sim, predList[-2], sd=T) +
        labs(fill="Std. Error") +
        facet_wrap(~Type))

#sourceCpp("./demo/dist.cpp")
load("./demo/prepData.rda")

nT <- 6

rangeE <- .5
covVal <- .2
covType <- "random"

set.seed(123)

field2 <- simField(
    N = as.matrix(select(fullDF, long, lat)), rangeE = rangeE,
    offset = c(0.1, 0.2),
    max.edge = c(0.18,0.2),
    beta0 = -2,
    sigmaE = .2,
    betaList = list(list(type=covType, value=covVal)),
    nTimes = nT,
    rho = .75, shape=spDF)

fullDF <- fullDF %>%
    rename(x=long, y=lat, oldid=id) %>%
    mutate(strat=paste(sprintf("%02d", reg), (2-urban), sep="_")) %>%
    left_join(
        select(filter(as_tibble(field2$spdf), tidx==0), x, y, id),
        by = c("x", "y"))

maxYear <- max(c(pointDF$year, polyDF$year))

yearWDF <- yearWDF %>%
    rename(oldid=id) %>%
    mutate(tidx=year-maxYear-1+nT) %>%
    left_join(select(fullDF, oldid, id), by="oldid") %>%
    filter(year <= maxYear) %>%
    select(-year, -oldid) %>%
    arrange(id) %>%
    filter(tidx>=0)

syearWDF <- filter(yearWDF, tidx == max(tidx))

(samplePlots$popDR <- syearWDF %>%
    left_join(select(fullDF, x, y, id), by="id") %>%
    ggplot(aes(x=x, y=y, fill=Population)) +
    geom_raster() +
    scale_fill_distiller(palette = "Spectral") +
    theme_void())

(samplePlots$fieldDR <- ggField(field2) +
        labs(fill="Probability"))

(samplePlots$fieldDR2 <- ggField(field2) +
        labs(fill="Probability") +
        geom_point(
            aes(x=long, y=lat),
            fill=NA,
            size=.1,
            data=pointDF %>%
                filter(grepl("2013", source)) %>%
                select(lat, long) %>%
                unique() %>%
                mutate(key=1) %>%
                left_join(tibble(tidx=0:4, key=1), by="key"))
        )

spDF <- readOGR("~/Documents/DRU5MR/data-extended/SECCenso2010.dbf") %>%
    spTransform(CRS("+proj=longlat"))
spDF$strat <- paste0(spDF$REG, "_", spDF$ZONA)
provShape <- st_as_sf(rgeos::gUnaryUnion(spDF, id=spDF@data$PROV))
regShape <- st_as_sf(rgeos::gUnaryUnion(spDF, id=spDF@data$REG))
regDivideShape <- st_as_sf(rgeos::gUnaryUnion(spDF, id=spDF@data$strat))

(samplePlots$reg <- ggplot(regShape) +
    geom_sf() +
    theme_void() +
    coord_sf(datum=NA))

(samplePlots$regUR <- ggplot(regDivideShape) +
    geom_sf(alpha=0) +
    theme_void() +
    coord_sf(datum=NA))

(samplePlots$popWDR <- syearWDF %>%
        left_join(select(fullDF, x, y, id), by="id") %>%
        ggplot(aes(x=x, y=y, fill=popW)) +
        geom_raster() +
        scale_fill_distiller(palette = "Spectral") +
        geom_sf(aes(x=NULL, y=NULL, fill=NULL), data=regDivideShape, alpha=0) +
        theme_void() +
        labs(fill="Weighted\nPopulation"))

(samplePlots$prov <- ggplot(provShape) +
    geom_sf() +
    theme_void() +
    coord_sf(datum=NA))

(samplePlots$regSamples <- ggplot(regShape) +
    geom_sf() +
    theme_void() +
    coord_sf(datum=NA) +
    geom_point(
        aes(x=long, y=lat),
        fill=NA,
        size=.1,
        data=pointDF %>%
            filter(grepl("2013", source)) %>%
            select(lat, long) %>%
            unique()))

sampDRResults <- 
    "~/Data/spaceTimeTest3/range=0.3,cov=-2,covtype=spatial,seed=1.Rds"
rezDR <- readRDS(sampDRResults)

(samplePlots$drResults <- ggFieldEst(rezDR$sim, rezDR$pred) +
    labs(fill="Probability"))
(samplePlots$drResultsPaper <- ggFieldEst(rezDR$sim, rezDR$pred) +
        labs(fill="Probability"))
(samplePlots$drSD <- ggFieldEst(rezDR$sim, rezDR$pred, sd=T) +
    labs(fill="Std. Error"))
(samplePlots$drSDPaper <- ggFieldEst(rezDR$sim, rezDR$pred[-2], sd=T) +
        labs(fill="Std. Error"))

(samplePlots$drProvError <- 
    do.call(sf:::rbind.sf, lapply(names(rezDR$pred), function(n){
        rezDR$provPred[[n]] %>%
            mutate(Model=n)})) %>%
        mutate(absDiff=abs(mu-trueValue)) %>%
        ggplot() +
        geom_sf(aes(fill = absDiff)) +
        theme_void() +
        coord_sf(datum=NA) +
        facet_grid(Model~tidx) +
        scale_fill_distiller(palette = "Spectral") +
        labs(fill="Absolute\nError"))

(samplePlots$drProvErrorPaper <- 
        do.call(sf:::rbind.sf, lapply(names(rezDR$pred), function(n){
            rezDR$provPred[[n]] %>%
                mutate(Model=n)})) %>%
        mutate(absDiff=abs(mu-trueValue)) %>%
        filter(Model != "Riemann") %>%
        ggplot() +
        geom_sf(aes(fill = absDiff)) +
        theme_void() +
        coord_sf(datum=NA) +
        facet_grid(Model~tidx) +
        scale_fill_distiller(palette = "Spectral") +
        labs(fill="Absolute\nError"))

samplePlots$demoPlotList <- list(
    example1 = plotList,
    example2 = plotList2,
    example3 = plotList3
)

samplePlots$samplesMap <-field2$bound %>%
    st_as_sf %>%
    ggplot() +
    geom_sf(alpha=0) +
    geom_point(
        aes(long, lat, color=Source),
        size=.2,
        alpha=.4,
        data=pointDF %>%
            select(lat, long, Source=source) %>%
            unique()) +
    theme_void() +
    coord_sf(datum=NA) +
    theme(legend.position = c(0.6, 0.2))

saveRDS(samplePlots, file="./demo/plotsForPresent.Rds")

pointDF %>%
    select(age_group, source, year) %>%
    bind_rows(select(polyDF, age_group, source, year)) %>%
    mutate(y=case_when(
        age_group == "NN" ~ 1/12,
        age_group == "PNN1" ~ 5/12,
        age_group == "PNN2" ~ 6/12,
        TRUE ~ 1
    )) %>%
    group_by(source, year) %>%
    summarize(personYears=sum(y)) %>%
    arrange(year, source) %>%
    as.data.frame

pointDF %>%
    select(age_group, source, year) %>%
    bind_rows(select(polyDF, age_group, source, year)) %>%
    mutate(y=case_when(
        age_group == "NN" ~ 1,
        age_group == "PNN1" ~ 5,
        age_group == "PNN2" ~ 6,
        TRUE ~ 12
    )) %>%
    group_by(source) %>%
    summarize(personYears=sum(y)) %>%
    arrange(source) %>%
    as.data.frame %>%
    pull(personYears)
