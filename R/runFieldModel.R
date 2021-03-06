#' Run a Model and Estimate Underlying Field From Point and Polygon Data
#'
#' @description Run a model and estimate underlying latent and probability
#' field from point and polygon data using TMB. 
#' 
#' @param field field object which simulated underlying data
#' @param pointDF data simulated from samplePoints
#' @param polyDF data simulated from samplePolygns
#' @param moption int, intger indicating how polygon data should be estimated 0
#' is by mixed model approximation, 1 is by redistribution, 2 is by Utazi
#' approach, 3 is by riemman approximation, 4 is ignoring polygon data,
#' 5 is if all data was known.
#' @param verbose logical, print model fitting information
#' @param symbolic logical, use metas reordering in model fitting
#' @param control list, control list passed to nlminb
#' @param rWidth int, only used in moption 2 to build besag prior
#' @param priors logical, default FALSE Should priors be evaluated for top level
#' parameters.
#' @param mcmc logical, default FALSE Should model be fit with MCMC. Not
#' compatible with moption 2.
#' @param AprojPoly sparseMatrix, sparse matrix with population weight information for polygons.
#' @param shape3 shape to build the adjaceny matrix for when moption == 2
#' @param start starting points of parameters.
#' @param ... Further arguments to pass to tmbstan
#' 
#' @return List of fitted model objects.
#' 
#' @import TMB
#' @useDynLib PointPolygon
#' 
#' @examples
#' \dontrun{
#' unitSim <- simField(
#' N = 500, rangeE = .7,
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
#' runFieldModel(unitSim, pointDF)
#'
#' }
#' 
#' @export

runFieldModel <- function(
    field, 
    pointDF = NULL,
    polyDF = NULL,
    moption = 0,
    verbose = FALSE,
    symbolic = TRUE,
    control = list(eval.max=1e4, iter.max=1e4),
    rWidth = NULL,
    priors = FALSE,
    mcmc = FALSE,
    AprojPoly = NULL,
    shape3 = NULL,
    start = list(),
    ...){
    model <- "PointPolygon"
    if(moption == 1 & !is.null(polyDF)){
        pointDF <- polyDF %>% 
            mutate(trueid=sapply(1:nrow(.), function(i){
                sample(polyDF$id[[i]], 1)})) %>%
            select(-id) %>%
            rename(id=trueid) %>%
            bind_rows(pointDF)
        polyDF <- NULL
    }
    if(is.null(polyDF)){
        moption <- 0
    }
    if(moption == 2){
        fit <- runFieldModelUtazi(
            field      = field, 
            pointDF    = pointDF,
            polyDF     = polyDF,
            moption    = moption,
            verbose    = verbose,
            symbolic   = symbolic,
            control    = control,
            rWidth     = rWidth,
            shape3     = shape3)
        
        return(fit)
    }
    if(moption == 5){
        moption <- 0
        pointDF <- polyDF %>%
            select(-id, -polyid) %>%
            rename(id=trueid) %>%
            select(id, tidx, trials, obs) %>%
            rbind(select(pointDF, id, tidx, trials, obs))
        polyDF <- NULL
    }
    moption_ <- ifelse(moption == 4, 0, moption)
    if(moption==4){
        polyDF <- NULL
    }
    inputs <- buildModelInputs(
        field,
        pointDF = pointDF,
        polyDF = polyDF,
        moption = moption_)
    inputs$Data$priors <- as.numeric(priors)
    if(!is.null(AprojPoly)){
        inputs$Data$AprojPoly <- AprojPoly
    }
    if(length(start) > 0){
       for(n in names(start)){
           if(n %in% names(inputs$Params)){
               if(verbose){
                   print(class(inputs$Params[[n]]))
                   print(paste0("replacing start value of ", n))
                   print(paste0("Original Length ", length(inputs$Params[[n]])))
                   print(paste0("New Length ", length(start[[n]])))
               }
               inputs$Params[[n]] <- start[[n]]
           }
       } 
    }
    startTime <- Sys.time()
    Obj <- TMB::MakeADFun(
        data = inputs$Data,
        parameters = inputs$Params,
        map=inputs$Map,
        DLL = model,
        random = "z",
        silent = !verbose)
    Obj$env$tracemgc <- verbose
    Obj$env$inner.control$trace <- verbose
    
    if(!mcmc){
        if(symbolic){
            nah <- utils::capture.output(TMB::runSymbolicAnalysis(Obj))
        }
        Opt <- stats::nlminb(
            start = Obj$par,
            objective = Obj$fn,
            gradient = Obj$gr,
            control = control)
        sdrep <- TMB::sdreport(Obj, getJointPrecision=TRUE)
        runtime <- Sys.time() - startTime
    }
    else{
        inputs$Data$priors <- as.numeric(1) # priors must be used for MCMC
        Opt <- tmbstan::tmbstan(Obj, ...)
        sdrep <- NULL
        runtime <- Sys.time() - startTime
    }

    list(
        obj = Obj,
        opt = Opt,
        runtime = runtime,
        sd = sdrep,
        moption = moption,
        stack=NULL)
}
