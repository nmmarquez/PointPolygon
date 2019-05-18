rm(list=ls())

paramDF <- expand.grid(
    rangeE = c(.3, .5, .7),
    covVal = c(2, .4, -.5, .2, -2),
    covType = c("random", "spatial", "cluster"),
    M = seq(50, 300, by=50),
    seed = 1:2)#1:10)

for(i in 1:nrow(paramDF)){
    modelname <- paste0(
        "range=", paramDF$rangeE[i],
        ",cov=", paramDF$covVal[i],
        ",covtype=", paramDF$covType[i],
        ",M=", paramDF$M[i],
        ",seed=", paramDF$seed[i])
    
    qsub <- paste(
        "qsub", 
        "-e ~/errors/",
        "-o ~/outputs/",
        "-l mem_free=20G -l m_mem_free=20G -P proj_geo_nodes_u5m",
        "-l fthread=10 -q geospatial.q",
        "-N", modelname,
        "/share/singularity-images/lbd/shells/singR.sh -m 2 -o 4 -e s",
        "~/Documents/PointPolygon/demo/replicateUtazi.R", 
        paramDF$rangeE[i],
        paramDF$covVal[i],
        paramDF$covType[i],
        paramDF$M[i],
        paramDF$seed[i],

        sep=" ")
    
    system(qsub)
}


paramDF2 <- dplyr::filter(unique(dplyr::select(paramDF, -M)), seed %in% 1:10)

for(i in 1:nrow(paramDF2)){
    modelname <- paste0(
        "range=", paramDF2$rangeE[i],
        ",cov=", paramDF2$covVal[i],
        ",covtype=", paramDF2$covType[i],
        ",seed=", paramDF2$seed[i])
    
    qsub <- paste(
        "qsub", 
        "-e ~/errors/",
        "-o ~/outputs/",
        "-l mem_free=100G -l m_mem_free=100G -P proj_geo_nodes_u5m",
        "-l fthread=20 -l h_rt=05:00:00:00 -q geospatial.q",
        "-N", modelname,
        "/share/singularity-images/lbd/shells/singR.sh -m 2 -o 4 -e s",
        "~/Documents/PointPolygon/demo/drTest.R", 
        paramDF$rangeE[i],
        paramDF$covVal[i],
        paramDF$covType[i],
        paramDF$seed[i],
        
        sep=" ")
    
    system(qsub)
}
