runFieldModel <- function(
    field, pointDF=NULL, polyDF=NULL, verbose=F, symbolic=T){
    
    inputs <- buildModelInputs(field, pointDF=pointDF, polyDF=polyDF)
    dyn.load(dynlib(model))
    startTime <- Sys.time()
    Obj <- MakeADFun(data=Data, parameters=Params, DLL=model, random="z",
                     silent=!verbose)
    Obj$env$tracemgc <- verbose
    Obj$env$inner.control$trace <- verbose
    symbolic <- T
    if(symbolic){
        nah <- capture.output(runSymbolicAnalysis(Obj))
    }
    Opt <- nlminb(
        start=Obj$par,
        objective=Obj$fn,
        gradient=Obj$gr,
        control=list(eval.max=1e4, iter.max=1e4))
    sdrep <- sdreport(Obj, getJointPrecision=TRUE)
    runtime <- Sys.time() - startTime
    dyn.load(dynlib(model))
    
    list(obj=Obj, opt=Opt, runtime=runtime, sd=sdrep)
}