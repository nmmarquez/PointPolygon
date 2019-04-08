#' Build model inputs for use with runFieldModel
#'
#' @description Build model inputs for use with runFieldModel
#' 
#' @param field field object which simulated underlying data
#' @param pointDF data simulated from samplePoints
#' @param polyDF data simulated from samplePolygns
#' @param model logical, build model inputs for modeling rather than prediction
#' @param moption int, intger indicating how polygon data should be estimated 0
#' is by Reimann sum approximation, 1 is by redistribution, and 2 is by Utazi
#' approach.
#' 
#' @return List of data and parameter objects to be used in TMB model run.
#' 
#' @export

buildModelInputs <- function(
    field,
    pointDF=NULL,
    polyDF=NULL,
    model=TRUE,
    moption=0){
    idPoly <- vector(mode="integer")
    if((moption == 1) & !is.null(polyDF)){
        polyDF <- dplyr::bind_rows(lapply(1:nrow(polyDF), function(i){
            ids <- polyDF[[i,"id"]][[1]]
            outs <- c(
                rep(1, polyDF[i,"obs"]),
                rep(0, polyDF[i,"trials"] - polyDF[i,"obs"]))
            data.frame(
                id = sample(ids, polyDF[i,"trials"], replace=TRUE),
                trials = 1,
                obs = sample(outs, polyDF[i,"trials"])
            )
        }))
        idPoly <- polyDF$id
    }
    if(is.null(polyDF)){
        empty <- vector("integer")
        polyDF <- data.frame(obs=empty, trials=empty, id=empty, polyid=empty)
    }
    if(is.null(pointDF)){
        empty <- vector("integer")
        pointDF <- data.frame(obs=empty, trials=empty, id=empty)
    }

    idx <- pointDF$id + 1 # have to add one for r index
    covDF <- dplyr::select(field$spdf@data, dplyr::matches("V[0-9]+"))

    if(ncol(covDF) == 1){
        if((nrow(polyDF) == 0) & model){
            covs <- matrix(covDF[idx,1], nrow=length(idx), ncol=1)
        }
        else{
            covs <- matrix(covDF[,1], nrow=nrow(covDF), ncol=1)
        }
    }
    else{
        if((nrow(polyDF) == 0) & model){
            covs <- as.matrix(covDF)[idx,]
        }
        else{
            covs <- as.matrix(covDF)
        }
    }

    if(nrow(polyDF) > 0){
        # hacky way of rearranging index for TMB
        polyIDX <- as.numeric(as.factor(polyDF$polyid)) -1
        uniquePolyid <- polyDF[!duplicated(polyDF$polyid), c("id", "polyid")]
        valID <- dplyr::bind_rows(lapply(1:nrow(uniquePolyid), function(i){
            pID <- uniquePolyid$id[[i]] + 1 # have to add one for r index
            data.frame(row=pID, col=i, val=1/length(pID))
        }))
        AprojPoly <- Matrix::sparseMatrix(
            i = valID$row,
            j = valID$col,
            x = valID$val,
            dims = c(nrow(field$spdf@data), nrow(uniquePolyid)))
        Data <- list(
            yPoint=pointDF$obs, denomPoint=pointDF$trials, idPoint=pointDF$id,
            yPoly=polyDF$obs, denomPoly=polyDF$trials, covs=covs,
            M0=field$spde$param.inla$M0,M1=field$spde$param.inla$M1,
            M2=field$spde$param.inla$M2, AprojObs=field$AprojField,
            AprojPoly=AprojPoly, moption=moption, idPoly=polyIDX)
        if(moption==3){
            reimDF <- polyDF %>% 
                select(-id, -trueid) %>%
                group_by(polyid) %>% 
                summarize_all(sum)
            Data$denomPoly <- reimDF$trials
            Data$yPoly <- reimDF$obs
        }
    }
    else{
        if(model & (moption == 0)){
            AprojObs <- field$AprojField[idx,]
        }
        else{
            AprojObs <- field$AprojField
        }
        AprojPoly <- Matrix::Matrix(data=0, nrow=0, ncol=0, sparse=TRUE)
        Data <- list(
            yPoint=pointDF$obs, denomPoint=pointDF$trials,
            idPoint=0:(nrow(pointDF) - 1),
            yPoly=polyDF$obs, denomPoly=polyDF$trials, covs=covs,
            M0=field$spde$param.inla$M0,M1=field$spde$param.inla$M1,
            M2=field$spde$param.inla$M2, AprojObs=AprojObs,
            AprojPoly=AprojPoly, moption=moption, idPoly=idPoly)
    }

    Params <- list(
        beta=0*(1:ncol(covs)), log_tau=0, log_kappa=0,z=rep(0, field$mesh$n)
    )

    return(list(Data=Data, Params=Params))
}
