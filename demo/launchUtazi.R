rm(list=ls())

paramDF <- expand.grid(
    rangeE = c(.3, .5, .7),
    covVal = 2#c(2, .4, -.5, .2, -2),
    covType = c("random", "spatial", "cluster"),
    M = 200#seq(50, 300, by=50),
    seed = 1:200)

for(i in 1:nrow(paramDF)){
    modelname <- paste0(
        "range=", paramDF$rangeE[i],
        ",cov=", paramDF$covVal[i],
        ",covtype=", paramDF$covType[i],
        ",M=", paramDF$M[i],
        ",seed=", paramDF$seed[i])
    
    qsub <- paste(
	# these all worked the last time I ran jobs
        "qsub", 
        "-e ~/errors/",
        "-o ~/outputs/",
        "-l mem_free=20G -l m_mem_free=20G -P proj_geo_nodes_u5m",
        "-l fthread=10 -q geospatial.q",
        "-N", modelname,
        "/share/singularity-images/lbd/shells/singR.sh -m 2 -o 4 -e s",
        # this will need to be changed
	"~/Documents/PointPolygon/demo/replicateUtazi.R", 
        paramDF$rangeE[i],
        paramDF$covVal[i],
        paramDF$covType[i],
        paramDF$M[i],
        paramDF$seed[i],

        sep=" ")

    system(qsub)
}


paramDF2 <- unique(dplyr::select(paramDF, -M))

#for(i in 1:nrow(paramDF2)){
#    modelname <- paste0(
#        "range=", paramDF2$rangeE[i],
#        ",cov=", paramDF2$covVal[i],
#        ",covtype=", paramDF2$covType[i],
#        ",seed=", paramDF2$seed[i])
#    
#    qsub <- paste(
#        "qsub", 
#        "-e ~/errors/",
#        "-o ~/outputs/",
#	# these all worked on the cluster when i last ran these
#        "-l mem_free=100G -l m_mem_free=100G -P proj_geo_nodes_u5m",
#        "-l fthread=20 -l h_rt=05:00:00:00 -q geospatial.q",
#        "-N", modelname,
#	# this may need to be chnaged
#        "/share/singularity-images/lbd/shells/singR.sh -m 2 -o 4 -e s",
#	# This will need to be changed
#        "~/Documents/PointPolygon/demo/drTest.R", 
#        paramDF2$rangeE[i],
#        paramDF2$covVal[i],
#        paramDF2$covType[i],
#        paramDF2$seed[i],
#        
#        sep=" ")
#
#        system(qsub)
#}

#for(i in 0:14){
#    modelname <- paste0("full_data_y", i)
#    
#    qsub <- paste(
#        "qsub", 
#        "-e ~/errors/",
#        "-o ~/outputs/",
#        "-l mem_free=200G -l m_mem_free=200G -P proj_geo_nodes_u5m",
#        "-l fthread=10 -l h_rt=05:00:00:00 -q geospatial.q",
#        "-N", modelname,
#        "/share/singularity-images/lbd/shells/singR.sh -m 10 -o 5 -e s",
#        "~/Documents/PointPolygon/demo/dataRun.R", i, "NA", sep=" ")
#    
#    system(qsub)
#}

#for(i in 0:9){
#    modelname <- paste0("full_data_reg", i)
#    
#    qsub <- paste(
#        "qsub", 
#        "-e ~/errors/",
#        "-o ~/outputs/",
#        "-l mem_free=200G -l m_mem_free=200G -P proj_geo_nodes_u5m",
#        "-l fthread=20 -l h_rt=05:00:00:00 -q geospatial.q",
#        "-N", modelname,
#        "/share/singularity-images/lbd/shells/singR.sh -m 10 -o 5 -e s",
#        "~/Documents/PointPolygon/demo/dataRun.R", "NA", i, sep=" ")
#    
#    system(qsub)
#}

#system(paste(
#    "qsub", 
#    "-e ~/errors/",
#    "-o ~/outputs/",
#    "-l mem_free=200G -l m_mem_free=200G -P proj_geo_nodes_u5m",
#    "-l fthread=20 -l h_rt=05:00:00:00 -q geospatial.q",
#    "-N", "full_data",
#    "/share/singularity-images/lbd/shells/singR.sh -m 10 -o 5 -e s",
#    "~/Documents/PointPolygon/demo/dataRun.R", "NA", "NA", sep=" "
#))

