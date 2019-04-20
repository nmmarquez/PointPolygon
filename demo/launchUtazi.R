rm(list=ls())

paramDF <- expand.grid(
    rangeE = c(.3, .5, .7),
    covVal = c(2, .4, -.5, .2, -2),
    covType = c("random", "spatial", "cluster"),
    M = seq(50, 300, by=50),
    seed = 1:10)

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
