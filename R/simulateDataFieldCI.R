#' Estimate Confidence Intervals (CI) for Probability Field Via Simulation for
#' data estimated models
#'
#' @description Calculate the mean, 95% confidence intervals, and the standard
#' deviation of an estimated probability field from a model. 
#' 
#' @param field field object which simulated underlying data
#' @param modelFit fitted model object from the field
#' @param draws int, number of draws for simulation to calculate CI's
#' 
#' @return data.frame, mean, CI, and sd of underlying field.
#' 
#' @import dplyr
#' @import tibble
#' @import sf
#' 
#' @export

simulateDataFieldCI <- function(modelFit, field, draws=1000){
  parDraws <- t(ar.matrix::sim.AR(draws, modelFit$sd$jointPrecision)) + 
    c(modelFit$opt$par, modelFit$sd$par.random)
  
  ageGN <- sum(names(modelFit$opt$par) == "beta_age")
  
  predDF <- bind_rows(lapply(0:ageGN, function(x){
    field$spdf %>%
      as_tibble() %>%
      select(id, tidx, urban) %>%
      mutate(aid=x)})) %>%
    mutate(mu=NA, sd=NA, lwr=NA, upr=NA)
  
  betaRows <- row.names(modelFit$sd$jointPrecision) == "beta"
  betaAgeRows <- row.names(modelFit$sd$jointPrecision) == "beta_age"
  zRows <- row.names(modelFit$sd$jointPrecision) == "z"
  phiRows <- row.names(modelFit$sd$jointPrecision) == "phi"
  betaDraws <- parDraws[betaRows, ]
  betaAgeDraws <- parDraws[betaAgeRows, ]
  zDraws <- parDraws[zRows, ]
  phiDraws <- parDraws[phiRows, ]
  nNod <- field$mesh$n
  projDraws <- t(do.call(rbind, lapply(1:field$nTimes, function(i){
    as.matrix(field$AprojField %*% zDraws[((i - 1) * nNod + 1):(i * nNod), ])
  })))
  
  for(i in 0:ageGN){
    t_ <- projDraws + betaDraws[1,] +
      as.matrix(betaDraws[2,])%*%t(filter(predDF, aid==i)$urban)
    if(i != 0){
      t_ <- t_ + betaAgeDraws[i,]
    }
    if(sum(phiRows) > 0){
      rIndex <- seq(i+1, ((ageGN + 1) * field$nTimes) - (ageGN - i), by=ageGN+1)
      phiADraws <- phiDraws[rIndex,]
      nPix <- nrow(filter(predDF, (aid == i) & (tidx==0)))
      for(j in 1:field$nTimes){
        j2 <- j-1
        t_[,(j2*nPix+1):((j2+1)*nPix)] <- t_[,(j2*nPix+1):((j2+1)*nPix)] +
          phiADraws[j,]
      }
    }
    
    predDraws <- arm::invlogit(t_)
    predDF$mu[predDF$aid == i] <- apply(predDraws, 2, mean)
    predDF$sd[predDF$aid == i] <- apply(predDraws, 2, sd)
    predDF$lwr[predDF$aid==i] <- apply(predDraws, 2, quantile, probs=.025)
    predDF$upr[predDF$aid==i] <- apply(predDraws, 2, quantile, probs=.975)
  }
  
  predDF
}
