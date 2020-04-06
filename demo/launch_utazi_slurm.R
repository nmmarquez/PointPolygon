# GOAL: take the original slurm scrip and replace the last value with a new
# simulation numnber for now later it will be replacing more details
# Must launch from the point polygon directory
slurm_script <- readLines("demo/test.slurm", -1)

launch <- paste0(
    "module load singularity gcc_8.2.1-ompi_3.1.4 && singularity ",
    "run --app Rscript singularity-r.simg ./replicateUtazi.R")

paramDF <- expand.grid(
    rangeE = .3,
    covVal = 2,
    covType = "spatial",
    M = as.integer(200),
    seed = 1:500)

for(i in 1:nrow(paramDF)){
    fn <- modelname <- paste0(
        "demo/slurm_files/",
        "range=", paramDF$rangeE[i],
        ",cov=", paramDF$covVal[i],
        ",covtype=", paramDF$covType[i],
        ",M=", paramDF$M[i],
        ",seed=", paramDF$seed[i], ".slurm")
    
    slurm_script_new <- slurm_script
    slurm_script_new[length(slurm_script_new)] <- paste(
        launch, paramDF$rangeE[i], paramDF$covVal[i], paramDF$covType[i],
        paramDF$M[i], paramDF$seed[i], sep = " "
    )
    
    writeLines(slurm_script_new, fn)
    
    system(paste0("sbatch ", fn))
    Sys.sleep(.5)
}
