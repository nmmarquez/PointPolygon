.libPaths(c("~/R3.5/", .libPaths()))
rm(list=ls())
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(PointPolygon)
library(sp)
library(Rcpp)
library(Matrix)
library(ggplot2)
setwd("~/Documents/PointPolygon/")
sourceCpp("./demo/dist.cpp")
load("./demo/prepData.rda")

nT <- 6

# for testing
rangeE <- .3
covVal <- 2
covType <- "random"
seed <- as.integer(1)

modelname <- paste0(
    "range=", rangeE,
    ",cov=", covVal,
    ",covtype=", covType,
    ",seed=", seed, ".Rds"
)

set.seed(seed)

find_closest <- function(x1, x2) {
    
    toRad <- pi / 180
    lat1  <- x1[,2]  * toRad
    long1 <- x1[,1] * toRad
    lat2  <- x2[,2]  * toRad
    long2 <- x2[,1] * toRad
    
    ord1  <- order(lat1)
    rank1 <- match(seq_along(lat1), ord1)
    ord2  <- order(lat2)
    
    ind <- find_closest_point(lat1[ord1], long1[ord1], lat2[ord2], long2[ord2])
    
    fullDF$id[ord2[ind + 1][rank1]]
}

maxYear <- max(c(pointDF$year, polyDF$year))

field <- simField(
    N = as.matrix(select(fullDF, long, lat)), rangeE = rangeE,
    offset = c(0.1, 0.2),
    max.edge = c(0.18,0.2),
    beta0 = -2,
    betaList = list(list(type=covType, value=covVal)),
    sigmaE = .2,
    nTimes = nT,
    rho = .85, shape=spDF)

ggField(field)

fullDF <- fullDF %>%
    rename(x=long, y=lat, oldid=id) %>%
    mutate(strat=paste(sprintf("%02d", reg), (2-urban), sep="_")) %>%
    left_join(
        select(filter(as_tibble(field$spdf), tidx==0), x, y, id),
        by = c("x", "y"))

fullDF %>%
    mutate(urban=urban*.8+.2, reg=as.factor(reg)) %>%
    ggplot(aes(x, y, fill=reg, alpha=urban)) +
    geom_raster() +
    theme_classic()

yearWDF <- yearWDF %>%
    rename(oldid=id) %>%
    mutate(tidx=year-maxYear-1+nT) %>%
    left_join(select(fullDF, oldid, id), by="oldid") %>%
    filter(year <= maxYear) %>%
    select(-year, -oldid) %>%
    arrange(id) %>%
    filter(tidx>=0)

syearWDF <- filter(yearWDF, tidx == max(tidx))

locDF <- syearWDF %>%
    group_by(strat) %>%
    summarize(id=list(id))

polyDF_ <- polyDF %>%
    mutate(tidx=year-maxYear-1+nT) %>%
    group_by(strat, tidx) %>%
    summarize(samples=sum(N)) %>%
    ungroup %>%
    filter(tidx >= 0) %>%
    right_join(
        polyDF %>%
            select(psu, strat) %>%
            unique() %>%
            mutate(trueid=sapply(strat, function(s){
                subDF <- filter(syearWDF, strat==s)
                sample(subDF$id, 1, prob=subDF$popW)
            })) %>%
            left_join(tibble(
                psu=rep(.$psu, nT),
                tidx=rep((1:nT)-1, each=nrow(.))
            ), by="psu"),
        by=c("strat", "tidx")) %>%
    rename(id=trueid) %>%
    # We dont need to weight the trials by the weight of the sample of the 
    # pixel as that has already been done in the PSU selection process
    # left_join(select(syearWDF, id, Population), by="id") %>%
    # group_by(strat, tidx) %>%
    # mutate(w=Population/sum(Population)) %>%
    group_by(strat, tidx) %>%
    mutate(w=1/n()) %>%
    ungroup %>%
    mutate(trials=round(w*samples)) %>%
    left_join(
        select(as_tibble(field$spdf), theta, id, tidx),
        by = c("id", "tidx")) %>%
    mutate(obs=rbinom(n(), trials, prob=theta)) %>%
    rename(trueid=id) %>%
    left_join(locDF, by="strat") %>%
    mutate(polyid=as.numeric(as.factor(strat))-1) %>%
    select(id, tidx, trials, obs, polyid, trueid, strat) %>%
    filter(trials != 0)

polyDF_ %>%
    group_by(tidx, polyid) %>%
    summarize(
        trials=sum(trials),
        obs=sum(obs),
        clusters=n()
    ) %>%
    ungroup %>%
    left_join(
        field$spdf %>%
            as_tibble() %>%
            select(-geometry) %>%
            left_join(select(syearWDF, -tidx)) %>%
            group_by(tidx, strat) %>%
            mutate(w=Population/sum(Population)) %>%
            summarize(theta=sum(w*theta)) %>%
            ungroup %>%
            arrange(tidx, strat) %>%
            mutate(polyid=rep((1:length(unique(.$strat)))-1, nT))
    ) %>%
    mutate(crude=obs/trials, absdiff=abs(crude-theta)) %>%
    arrange(-absdiff) %>%
    left_join(
        fullDF %>%
            group_by(strat) %>%
            summarize(cells=n())
    )
    # mutate(ccratio = clusters/ cells) %>%
    # ggplot(aes(x=ccratio, y=absdiff)) +
    # geom_point() +
    # theme_classic()

# we did really bad on this one example where we had seemingly a large number
# of clusters. I want to see the odds of getting the observed given the true
# field

polyid__ = 19
tidx__ = 1
strat__ = "10_2"
crude__ = .315
true__ = .260
clusters__ = 36

# this doesnt match up with what i tend to think of for central limit stuff
# im going to nee to do more sanity checking

testDF <- select(as_tibble(field$spdf), id, tidx, theta) %>%
    left_join(select(syearWDF, -tidx), by="id") %>%
    filter(tidx == tidx__ & strat == strat__) %>%
    mutate(w=Population/sum(Population))

# these are from the actual observed sample
(sampleP <- polyDF_ %>%
    filter(tidx == tidx__ & strat == strat__) %>%
    select(tidx, id=trueid, obsold=obs, trials) %>%
    left_join(select(syearWDF, -tidx), by="id") %>%
    left_join(
        select(as_tibble(field$spdf), id, tidx, theta),
        by=c("id", "tidx")) %>%
    mutate(w=Population/sum(Population)) %>%
    {mean(.$theta)})

# this replicates true__ above
testDF %>%
    {sum(.$w * .$theta)}

# How wide can these samples get???
sapply(1:1000, function(i){
    testDF %>%
        sample_n(clusters__, replace=T, w) %>%
        mutate(w=Population/sum(Population)) %>%
        {mean(.$theta)}
    }) %>%
    {tibble(x=.)} %>%
    ggplot(aes(x)) +
    geom_density() +
    geom_vline(xintercept=crude__, color="red", linetype=2) +
    geom_vline(xintercept=true__, color="black")

tibble(x=100:200, density=dbinom(100:200, 432, sampleP)) %>%
    ggplot(aes(x=x, y=density)) +
    geom_line() +
    theme_classic() +
    geom_vline(xintercept=136, color="black")

stratOrder <- arrange(
    as.data.frame(unique(select(polyDF_, polyid, strat))), polyid)

pointDF_ <- pointDF %>%
    rename(oldid=id) %>%
    mutate(tidx=year-maxYear-1+nT) %>%
    filter(tidx>=0) %>%
    left_join(select(fullDF, oldid, id), by="oldid") %>%
    left_join(
        select(as_tibble(field$spdf), theta, id, tidx),
        by = c("id", "tidx")) %>%
    rename(trials=N) %>%
    mutate(obs=rbinom(n(), trials, prob=theta)) %>%
    arrange(tidx, id) %>%
    select(id, tidx, trials, obs)

# Turns out IHME only has one point or each large admin level :/
ihmePolyDF <- read_csv("./demo/ihmeResampleDF.csv") %>%
    select(longitude, latitude) %>%
    unique %>%
    as.matrix %>%
    find_closest(
        syearWDF %>%
            left_join(
                select(as_tibble(field$spdf), -geometry), 
                by = c("tidx", "id")) %>%
            select(x, y) %>%
            as.matrix()
    ) %>%
    tibble(id=.) %>%
    mutate(id=id-1) %>%
    left_join(syearWDF, by="id") %>%
    rename(trueid=id) %>%
    mutate(location_code=as.numeric(str_split(strat, "_", simplify=T)[,1])) %>%
    select(trueid, location_code) %>%
    right_join(
        polyDF_ %>%
            mutate(location_code=
                       as.numeric(str_split(strat, "_", simplify=T)[,1])) %>%
            select(-trueid) , 
        by="location_code") %>%
    group_by(trueid, tidx, location_code) %>%
    summarise(
        trials=sum(trials), obs=sum(obs), id=id[1]) %>%
    ungroup() %>%
    mutate(polyid=as.numeric(as.factor(trueid))-1) %>%
    select(id, tidx, trials, obs, polyid, trueid) %>%
    arrange(tidx, polyid)

AprojPoly <- syearWDF %>%
    left_join(stratOrder, by="strat") %>%
    mutate(col=polyid+1, row=id+1) %>%
    {sparseMatrix(
        i = .$row,
        j = .$col,
        x = .$popW,
        dims = c(max(.$row), max(.$col)))}

### The approaches that we want to test are as follows
# 1 Mixture model
# 2 IHME Resample
# 3 Riemman
# 4 Ignore
# 5 Known

modelList <- list(
    `IHME Resample` = runFieldModel(
        field, pointDF_, ihmePolyDF, moption=5, verbose=T
    ),
    `Riemann` = runFieldModel(
        field, pointDF_, select(polyDF_, -strat), moption=3, verbose=T,
        AprojPoly=AprojPoly
    ),
    `Ignore` = runFieldModel(
        field, pointDF_, polyDF_, moption=4, verbose=T, AprojPoly=AprojPoly
    ),
    `Known` = runFieldModel(
        field, pointDF_, polyDF_, moption=5, verbose=T, AprojPoly=AprojPoly
    )
)

modelList[["Mixture Model"]] <- runFieldModel(
    field, pointDF_, polyDF_, moption=0, verbose=T, AprojPoly=AprojPoly,
    control = list(eval.max = 100000, iter.max = 100000),
    start = list(
        z=matrix(unname(modelList$`IHME Resample`$sd$par.random), ncol=nT),
        beta=modelList$`IHME Resample`$opt$par[names(
            modelList$`IHME Resample`$opt$par) == "beta"])
)

unitFitList <- lapply(modelList, function(model){
    simulateFieldCI(field, model)
})

sapply(unitFitList, function(df){
    sapply(0:5, function(i){
        df %>%
            filter(tidx == i) %>%
            {sqrt(mean((.$mu - .$trueValue)^2))}
    })
})

sapply(unitFitList, function(df){
    sqrt(mean((df$mu - df$trueValue)^2))
})

# unitBetaList <- lapply(modelList, function(model){
#     simulateBetas(field, model)
# })

convergeList <- lapply(modelList, function(model){
    model$opt$convergence
})

timeList <- lapply(modelList, function(model){
    as.numeric(model$runtime, units="mins")
})

aggShapes <- readRDS("demo/aggShapes.Rds")

provPredList <- lapply(modelList, function(model){
    arealCI(
        field, 
        model, 
        polygonList = lapply(1:nrow(aggShapes$provShape@data), function(i){
            sf::st_as_sf(aggShapes$provShape[i,])}),
        popDF = mutate(select(syearWDF, -tidx), w=Population))
})

regPredList <- lapply(modelList, function(model){
    arealCI(
        field, 
        model, 
        polygonList = lapply(1:nrow(aggShapes$regShape@data), function(i){
            sf::st_as_sf(aggShapes$regShape[i,])}),
        popDF = mutate(select(syearWDF, -tidx), w=Population))
})

sapply(provPredList, function(df){
    sapply(0:5, function(i){
        df %>%
            filter(tidx == i) %>%
            {sqrt(mean((.$mu - .$trueValue)^2))}
    })
})

cbind(sapply(modelList, function(model){
    model$opt$par
}), True=c(2,covVal, field$fieldPars))


ggFieldEst(field, unitFitList)
ggFieldEst(field, unitFitList, sd=T)

