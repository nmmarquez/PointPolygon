


ggFieldEst <- function(field, predList){
    DF <- field$spdf@data
    allDF <- DF
    allDF$Type <- "True"
    for(i in names(predList)){
        DFpred <- DF
        DFpred$Type <- i
        DFpred$theta <- predList[[i]]$mu
        allDF <- dplyr::bind_rows(allDF, DFpred)
    }
    
    new <- field
    new$spdf@data <- allDF
    ggField(new) +
        facet_wrap(~Type)
}
