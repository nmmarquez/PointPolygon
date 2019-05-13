#' Estimate Confidence Intervals (CI) for Areal Probability Field Via Simulation
#'
#' @description Calculate the mean, 95% confidence intervals, and the standard
#' deviation of areal units from an estimated probability field from a model. 
#' 
#' @param field field object which simulated underlying data
#' @param modelFit fitted model object from the field
#' @param draws int, number of draws for simulation to calculate CI's
#' @param polygonList list, list of polygons to sample from
#' @param rWidth integer, instead of using a polygon list divide the original
#' bounding box into approximately equal sized squares with rWidth number
#' squares along the x axis.
#' @param popDF data.frame, use population weights when aggrgating data
#' @return data.frame, mean, CI, and sd of areal fields.
#' 
#' @examples
#' \dontrun{
#' 
#' unitSim <- simField(
#' N = 200, rangeE = .7,
#' offset = c(0.1, 0.2), 
#' max.edge = c(0.1,0.2),
#' beta0 = -2,
#' betaList = list(
#'     list(type="random", value=2),
#'     list(type="spatial", value=-.5),
#'     list(type="cluster", value=-2)
#' ))
#'
#' pointDF <- samplePoints(unitSim, 500, 100)
#'
#' unitFit <- runFieldModel(unitSim, pointDF)
#'
#' head(arealCI(unitSim, unitFit, rWidth=3))
#'
#' }
#' 
#' @export

arealCI <- function(field, modelFit, polygonList=NULL, rWidth=NULL, draws=1000, popDF=NULL){
    if(is.null(polygonList)){
        if(is.null(rWidth)){
            stop("rWidth and polygonList can not both be NULL")
        }
        sectionedSPDF <- dividePolygon(field$bound, rWidth)
        polygonList <- lapply(1:nrow(sectionedSPDF@data), function(i){
            sf::st_as_sf(sectionedSPDF[i,])
        })
    }
    if(modelFit$moption == 2){
        if(field$nTimes > 1){
            stop("Utazi Model Only Supports Single Year Analysis Currently.")
        }
        post.samples <- inla.posterior.sample(1000, modelFit)
        index.pred <- inla.stack.index(modelFit$stack, "pred")$data
        pSamples <- arm::invlogit(sapply(post.samples, function(z){
            z$latent[index.pred,]}))
        resultsDF <- data.frame(
            mu = apply(pSamples, 1, mean),
            sd = apply(pSamples, 1, sd),
            lwr = apply(pSamples, 1, quantile, probs=.025),
            upr = apply(pSamples, 1, quantile, probs=.975),
            id = field$spdf$id,
            trueValue = field$spdf$theta
        )
        row.names(resultsDF) <- NULL
        return(resultsDF)
    }
    if(is.null(popDF)){
        w <- rep(1, nrow(field$spdf))
    }
    else{
        w <- field$spdf %>%
            as_tibble %>%
            left_join(popDF) %>%
            pull(w)
    }
    modelLRP <- buildModelInputs(field, model=FALSE)
    
    parDraws <- t(ar.matrix::sim.AR(draws, modelFit$sd$jointPrecision)) +
        c(modelFit$opt$par, modelFit$sd$par.random)
    
    betaRows <- row.names(modelFit$sd$jointPrecision) == "beta"
    zRows <- row.names(modelFit$sd$jointPrecision) == "z"
    
    betaDraws <- parDraws[betaRows,]
    zDraws <- parDraws[zRows,]
    nNod <- field$mesh$n
    
    fieldProbs <- arm::invlogit(
        modelLRP$Data$covs %*% betaDraws +
            do.call(rbind, lapply(1:field$nTimes, function(i){as.matrix(
                modelLRP$Data$AprojObs %*% zDraws[((i-1)*nNod + 1):(i*nNod),])
            })))
    
    ciDF <- bind_rows(lapply(1:length(polygonList), function(i){
        subSPDF <- polygonList[[i]]
        subSPDF$isPresent <- TRUE
        pointsDF <- cbind(
            isPresent=sapply(sf::st_intersects(field$spdf,subSPDF), function(z){
                length(z)>0}), 
            field$spdf) %>%
            filter(isPresent) %>%
            dplyr::as_tibble() %>%
            select(tidx, id)
        
        bind_rows(lapply(unique(pointsDF$tidx), function(t){
            ridx <- pointsDF %>%
                filter(tidx == t) %>%
                mutate(present=T) %>%
                right_join(
                    dplyr::as_tibble(field$spdf),
                    by = c("tidx", "id")) %>%
                mutate(present=!is.na(present)) %>%
                pull(present) %>%
                which()
            
            subFieldProbs <- apply(fieldProbs[ridx,], 2, weighted.mean, w=w[ridx])
            
            data.frame(
                mu = mean(subFieldProbs),
                sd = stats::sd(subFieldProbs),
                lwr = stats::quantile(subFieldProbs, probs=.025),
                upr = stats::quantile(subFieldProbs, probs=.975),
                polyid = i-1,
                tidx = t,
                trueValue = weighted.mean(field$spdf$theta[ridx], w=w[ridx])
            )
        }))
    }))
    
    do.call(rbind, polygonList) %>%
        mutate(polyid=(1:n())-1) %>%
        right_join(ciDF, by="polyid")
}
