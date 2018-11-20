#' Run a Model and Estimate Underlying Field From Point and Polygon Data
#'
#' @description Run a model and estimate underlying latent and probability
#' field from point and polygon data using TMB. 
#' 
#' @param field field object which simulated underlying data
#' @param pointDF data simulated from samplePoints
#' @param polyDF data simulated from samplePolygns
#' @param verbose logical, print model fitting information
#' @param symbolic logical, use metas reordering in model fitting
#' @param modelpath character, folder where model file resides
#' @param control list, control list passed to nlminb
#' @param recompile logical, rebuild model for testing purposes only
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
    verbose = FALSE,
    symbolic = TRUE,
    modelpath = "~/Documents/PointPolygon/isrc/",
    control = list(eval.max=1e4, iter.max=1e4),
    recompile = FALSE){
    model <- "model"
    fullPath <- paste0(modelpath, model)
    if(recompile){
        TMB::compile(paste0(fullPath, ".cpp"))
    }

    inputs <- buildModelInputs(field, pointDF=pointDF, polyDF=polyDF)
    dyn.load(TMB::dynlib(fullPath))
    startTime <- Sys.time()
    Obj <- TMB::MakeADFun(
        data = inputs$Data,
        parameters = inputs$Params,
        DLL = model,
        random = "z",
        silent = !verbose)
    Obj$env$tracemgc <- verbose
    Obj$env$inner.control$trace <- verbose
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
    nah <- utils::capture.output(dyn.unload(TMB::dynlib(fullPath)))

    list(obj=Obj, opt=Opt, runtime=runtime, sd=sdrep)
}
