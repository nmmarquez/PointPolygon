#' Estimate Confidence Intervals (CI) for Probability Field Via Simulation for
#' data estimated models
#'
#' @description Calculate the mean, 95% confidence intervals, and the standard
#' deviation of an estimated probability field from a model. 
#' 
#' @param field field object which simulated underlying data
#' @param modelFit fitted model object from the field
#' @param draws int, number of draws for simulation to calculate CI's
#' @param agg5q0 logical, collpase age groups to 5q0?
#' @param popDF data.frame, use population weights when aggrgating data
#' 
#' @return data.frame, mean, CI, and sd of underlying field.
#' 
#' @import dplyr
#' @import tibble
#' @import sf
#' 
#' @export

arealDataCI <- function(modelFit, field, draws=1000, agg5q0=FALSE, popDF=NULL, polygonList=NULL){
  if(is.null(popDF)){
    w <- rep(1, nrow(field$spdf))
  }
  else{
    w <- field$spdf %>%
      as_tibble %>%
      left_join(popDF) %>%
      pull(w)
  }
  
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
  
  fieldDraws <- 1
  ciDFAll <- tibble()
  
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
    
    if(!agg5q0){
      fieldDraws <- t(arm::invlogit(t_))
      
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
          
          subFieldProbs <- apply(fieldDraws[ridx,], 2, weighted.mean, w=w[ridx])
          
          tibble(
            mu = mean(subFieldProbs),
            sd = stats::sd(subFieldProbs),
            lwr = stats::quantile(subFieldProbs, probs=.025),
            upr = stats::quantile(subFieldProbs, probs=.975),
            polyid = i-1,
            tidx = t
          )}))
        })) %>%
        mutate(aid=i)
      
      ciDFAll <- bind_rows(ciDFAll, ciDF)
      
    }
    else{
      fieldDraws <- fieldDraws * (1-t(arm::invlogit(t_)))
    }
  }
  
  if(agg5q0){
    
    fieldDraws <- 1 - fieldDraws
    
    ciDFAll <- bind_rows(lapply(1:length(polygonList), function(i){
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
        
        subFieldProbs <- apply(fieldDraws[ridx,], 2, weighted.mean, w=w[ridx])
        
        tibble(
          mu = mean(subFieldProbs),
          sd = stats::sd(subFieldProbs),
          lwr = stats::quantile(subFieldProbs, probs=.025),
          upr = stats::quantile(subFieldProbs, probs=.975),
          polyid = i-1,
          tidx = t
        )}))
    })) %>%
      mutate(aid=99)
  }
  
  ciDFAll
}
