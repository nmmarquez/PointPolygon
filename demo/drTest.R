rm(list=ls())
library(tidyverse)
library(PointPolygon)
library(sp)
library(Rcpp)
library(Matrix)
set.seed(12345)
setwd("~/Documents/PointPolygon/")
sourceCpp("./demo/dist.cpp")
load("./demo/prepData.rda")

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

field <- simField(
    N = as.matrix(select(fullDF, long, lat)), rangeE = .4,
    offset = c(0.1, 0.2),
    max.edge = c(0.1,0.2),
    beta0 = -2,
    #betaList = list(list(type="spatial", value=1)),
    nTimes = 15,
    rho = .85, shape=spDF)

ggField(field)

fullDF <- fullDF %>%
    rename(x=long, y=lat, oldid=id) %>%
    mutate(strat=paste(sprintf("%02d", reg), (2-urban), sep="_")) %>%
    left_join(
        select(filter(as_tibble(field$spdf), tidx==0), x, y, id),
        by = c("x", "y"))

yearWDF <- yearWDF %>%
    rename(oldid=id) %>%
    mutate(tidx=year-min(year)) %>%
    left_join(select(fullDF, oldid, id), by="oldid") %>%
    select(-year, -oldid) %>%
    arrange(id)

syearWDF <- filter(yearWDF, tidx == 14)

locDF <- syearWDF %>%
    group_by(strat) %>%
    summarize(id=list(id))

# for each psu resample points based off of the population weights of locations
psuDF <- polyDF %>%
    select(psu, strat) %>%
    unique() %>%
    mutate(trueid=sapply(strat, function(s){
        subDF <- filter(syearWDF, strat==s)
        sample(subDF$id, 1, prob=subDF$popW)
    }))

polyDF <- polyDF %>%
    mutate(tidx=year-min(year)) %>%
    group_by(psu, tidx) %>%
    summarize(trials=sum(N)) %>%
    ungroup %>%
    left_join(psuDF, by="psu") %>%
    rename(id=trueid) %>%
    # lets do constant population for now
    left_join(
        select(filter(as_tibble(field$spdf), tidx==14), theta, id),
        by = c("id")) %>%
    mutate(obs=rbinom(n(), trials, prob=theta)) %>%
    rename(trueid=id) %>%
    left_join(locDF, by="strat") %>%
    mutate(polyid=as.numeric(as.factor(strat))-1) %>%
    select(id, tidx, trials, obs, polyid, trueid, strat)

stratOrder <- arrange(
    as.data.frame(unique(select(polyDF, polyid, strat))), polyid)

pointDF <- pointDF %>%
    rename(oldid=id) %>%
    mutate(tidx=year-min(year)) %>%
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
        polyDF %>%
            mutate(location_code=as.numeric(str_split(strat, "_", simplify=T)[,1])) %>%
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

# TODO: We need to pass a custom projection matrixc toaccount for the correct pop weights

modelList <- list(
    `Mixture Model` = runFieldModel(
        field, pointDF, polyDF, moption=0, verbose=T
    ),
    `IHME Resample` = runFieldModel(
        field, pointDF, ihmePolyDF, moption=5, verbose=T
    ),
    `Mixture Model` = runFieldModel(
        field, pointDF, polyDF, moption=3, verbose=T
    ),
    `Mixture Model` = runFieldModel(
        field, pointDF, polyDF, moption=4, verbose=T
    ),
    `Mixture Model` = runFieldModel(
        field, pointDF, polyDF, moption=5, verbose=T
    )
)
