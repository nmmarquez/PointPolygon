# because of the way that the csde slurm machines are set up we are in a 
# better spot if we request a whole node and run in parallel from there. This 
# script is thus a wrapper around launch Utazi which creates a bunch of those
# jobs in parallel.
library(parallel)

args <- commandArgs(trailingOnly=TRUE)
rangeE <- as.numeric(args[1]) # range of spatial prces varies from {.3, .5, .7}
covVal <- as.numeric(args[2]) # covariate effect in set {.2, .4, -.5, .2, -2}
covType <- args[3] # either random spatial or cluster 
M <- as.integer(args[4]) # number of samples in Polygons chosen from U(50, 300)
s1 <- as.integer(args[5]) # RNG
seeds <- (((s1 - 1) * 20) + 1):(s1*20)

mclapply(seeds, function(seed){
    source("./replicateUtazi.R", local = TRUE)
}, mc.cores = 20)