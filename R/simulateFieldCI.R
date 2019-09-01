#' Estimate Confidence Intervals (CI) for Probability Field Via Simulation
#'
#' @description Calculate the mean, 95% confidence intervals, and the standard
#' deviation of an estimated probability field from a model. 
#' 
#' @param field field object which simulated underlying data
#' @param modelFit fitted model object from the field
#' @param draws int, number of draws for simulation to calculate CI's
#' @param chain int, the chain to use if drawing inference from MCMC model
#' 
#' @return data.frame, mean, CI, and sd of underlying field.
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
#' head(simulateFieldCI(unitSim, unitFit))
#'
#' }
#' 
#' @export

simulateFieldCI <- function(field, modelFit, draws=1000, chain=1){
    if(modelFit$moption == 2){
        post.samples <- inla.posterior.sample(draws, modelFit)
        index.pred <- inla.stack.index(modelFit$stack, "pred")$data
        pSamples <- arm::invlogit(sapply(post.samples, function(z){
            z$latent[index.pred,]}))
        resultsDF <- data.frame(
            mu = apply(pSamples, 1, mean),
            sd = apply(pSamples, 1, sd),
            lwr = apply(pSamples, 1, quantile, probs=.025),
            upr = apply(pSamples, 1, quantile, probs=.975),
            id = field$spdf$id,
            tidx = field$spdf$tidx,
            trueValue = field$spdf$theta
        )
        row.names(resultsDF) <- NULL
        return(resultsDF)
    }
    
    modelLRP <- buildModelInputs(field, model=FALSE)
    
    if(class(modelFit$opt) == "stanfit"){
        sims <- modelFit$opt@sim$samples[[chain]]
        miter <- modelFit$opt@stan_args[[chain]]$iter
        dstart <- miter - draws + 1

        betaRows <- startsWith(names(sims), "beta[")
        zRows <- startsWith(names(sims), "z[")

        betaDraws <- do.call(rbind, sims[betaRows])[,dstart:miter]
        zDraws <- do.call(rbind, sims[zRows])[,dstart:miter]
    }

    else{
        parDraws <- t(ar.matrix::sim.AR(draws, modelFit$sd$jointPrecision)) +
            c(modelFit$opt$par, modelFit$sd$par.random)
        
        betaRows <- row.names(modelFit$sd$jointPrecision) == "beta"
        zRows <- row.names(modelFit$sd$jointPrecision) == "z"
        
        betaDraws <- parDraws[betaRows,]
        zDraws <- parDraws[zRows,]
    }
    nNod <- field$mesh$n

    fieldProbs <- arm::invlogit(modelLRP$Data$covs %*% betaDraws +
        do.call(rbind, lapply(1:field$nTimes, function(i){
            as.matrix(
                modelLRP$Data$AprojObs %*% zDraws[((i-1)*nNod + 1):(i*nNod),])
        })))

    data.frame(
        mu = apply(fieldProbs, 1, mean),
        sd = apply(fieldProbs, 1, stats::sd),
        lwr = apply(fieldProbs, 1, stats::quantile, probs=.025),
        upr = apply(fieldProbs, 1, stats::quantile, probs=.975),
        id = field$spdf$id,
        tidx = field$spdf$tidx,
        trueValue = field$spdf$theta
        )
}
