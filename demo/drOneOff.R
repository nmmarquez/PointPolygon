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
setwd("~/Documents/PointPolygon/")
sourceCpp("./demo/dist.cpp")
load("./demo/prepData.rda")

# nT <- 6

# for testing
rangeE <- .3
covVal <- 2
covType <- "random"
seed <- as.integer(123)

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
    #betaList = list(list(type=covType, value=covVal)),
    nTimes = 1,
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
    #mutate(tidx=year-maxYear-1+nT) %>%
    mutate(tidx=ifelse(year==2013, 0, -1)) %>%
    left_join(select(fullDF, oldid, id), by="oldid") %>%
    select(-year, -oldid) %>%
    arrange(id) %>%
    filter(tidx>=0)

syearWDF <- filter(yearWDF, tidx == max(tidx))

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
    #mutate(tidx=year-maxYear-1+nT) %>%
    mutate(tidx=ifelse(year==2013, 0, -1)) %>%
    filter(tidx>=0) %>%
    group_by(psu, tidx) %>%
    summarize(trials=sum(N)) %>%
    ungroup %>%
    left_join(psuDF, by="psu") %>%
    rename(id=trueid) %>%
    # lets do constant population for now
    left_join(
        select(filter(as_tibble(field$spdf), tidx==max(tidx)), theta, id),
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
    #mutate(tidx=year-maxYear-1+nT) %>%
    mutate(tidx=ifelse(year==2013, 0, -1)) %>%
    filter(tidx>=0) %>%
    left_join(select(fullDF, oldid, id), by="oldid") %>%
    left_join(
        select(as_tibble(field$spdf), theta, id, tidx),
        by = c("id", "tidx")) %>%
    rename(trials=N) %>%
    mutate(obs=rbinom(n(), trials, prob=theta)) %>%
    arrange(tidx, id) %>%
    select(id, tidx, trials, obs) %>%
    group_by(id, tidx) %>%
    summarize_all(sum) %>%
    ungroup

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

# lets try the polygons again but make the units the superstates
# syearSSWDF <- syearWDF %>%
#     mutate(location_code=as.numeric(str_split(strat, "_", simplify=T)[,1])) %>%
#     mutate(polyid=location_code-1) %>%
#     select(-location_code) %>%
#     group_by(polyid) %>%
#     mutate(popW=Population/sum(Population))
# 
# polyDF2 <- polyDF %>%
#     mutate(location_code=as.numeric(str_split(strat, "_", simplify=T)[,1])) %>%
#     mutate(polyid=location_code-1)

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

# modelList <- list(
#     `Sub` = runFieldModel(
#         field, dplyr::sample_frac(pointDF, size=.6), moption=0, verbose=T
#     ),
#     `All` = runFieldModel(
#         field, pointDF, moption=0, verbose=T
#     ),
#     `IHME Resample` = runFieldModel(
#         field, pointDF, ihmePolyDF, moption=5, verbose=T, AprojPoly=AprojPoly
#     ),
#     `Known` = runFieldModel(
#         field, pointDF, polyDF, moption=5, verbose=T, AprojPoly=AprojPoly
#     )
# )

modelList <- list(
    # `Mixture Model` = runFieldModel(
    #     field, pointDF, polyDF, moption=0, verbose=T, AprojPoly=AprojPoly,
    #     control = list(eval.max = 100000, iter.max = 100000)
    # ),
    `IHME Resample` = runFieldModel(
        field, pointDF, ihmePolyDF, moption=5, verbose=T, AprojPoly=AprojPoly
    ),
    `Riemann` = runFieldModel(
        field, pointDF, select(polyDF, -strat), moption=3, verbose=T, AprojPoly=AprojPoly
    ),
    `Ignore` = runFieldModel(
        field, pointDF, polyDF, moption=4, verbose=T, AprojPoly=AprojPoly
    ),
    `Known` = runFieldModel(
        field, pointDF, polyDF, moption=5, verbose=T, AprojPoly=AprojPoly
    )
)

modelList[["Mixture Model"]] <- runFieldModel(
        field, pointDF, polyDF, moption=0, verbose=T, AprojPoly=AprojPoly,
        control = list(eval.max = 100000, iter.max = 100000),
        start = list(
            z=matrix(unname(modelList$`IHME Resample`$sd$par.random), ncol=1),
            beta=modelList$`IHME Resample`$opt$par["beta"])
    )

modelList$Riemann <- NULL

unitFitList <- lapply(modelList, function(model){
    simulateFieldCI(field, model)
})

unitBetaList <- lapply(modelList, function(model){
    simulateBetas(field, model)
})

convergeList <- lapply(modelList, function(model){
    model$opt$convergence
})

timeList <- lapply(modelList, function(model){
    as.numeric(model$runtime, units="mins")
})

oneOffResults <- list()

oneOffResults$modelTable <- tibble(
    model = names(unitFitList),
    # rmse of pixels
    rmse = sapply(unitFitList, function(df){
        sqrt(mean((df$mu - df$trueValue)^2))
    }),

    # Confidence intervals coverage
    coverage = sapply(unitFitList, function(x){
        mean(x$trueValue > x$lwr & x$trueValue < x$upr)
    }),

    # Std err of estimates
    `Std. Error` = sapply(unitFitList, function(x){
        mean(x$sd)
    })) %>% arrange(rmse)



oneOffResults$field <- ggField(field) +
    ggplot2::theme(
        strip.background = ggplot2::element_blank(),
        strip.text.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(fill="Probability", title="True Probability Field")

oneOffResults$fieldEst <- ggFieldEst(field, unitFitList) +
    ggplot2::facet_wrap(~Type) +
    ggplot2::labs(fill="Probability", title="Probability Field Estimates")

oneOffResults$fieldSD <- ggFieldEst(field, unitFitList, sd=T) +
    ggplot2::facet_wrap(~Type) +
    ggplot2::labs(fill="Probability\nStd. Err.", title="Probability Field Estimates")

oneOffResults$unitResults <- list(
    sim = field,
    pred = unitFitList,
    # model = modelList, # this takes up a lot of space so ignore for now
    converge = convergeList,
    covType = covType,
    rangeE = rangeE,
    covVal = covVal,
    seed = seed,
    runtime = timeList
)

saveRDS(oneOffResults, file=paste0("./demo/oneOff.Rds"))
