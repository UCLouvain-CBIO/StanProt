##' @include ProtModel.R
NULL

##' Take an initialised `ProtModel` object, runs the inference engine and 
##' returns an updated `ProtModel` object containing the estimated coefficients.
##' 
##' @title Run a proteomics model
##' 
##' @param x A `ProtModel` object as produced by [makeModel()].
##' 
##' @param modelScript `character(1)` provide the name without extension
##'     to a stan model file. The function will check if either
##'     `stanmodel.rds` or `stanmodel.stan` are available an use the
##'     more recent one. Only used in 'stan' mode
##'     
##' @param mode `character(1)` can be either 'stan' or 'native'
##' 
##' @param algo `character(1)` can be :
##' - 'sampling' for exact inference,
##' - 'mfvb' for approximate inference using mean field variational bayes,
##' - 'MAP' for maximum a posteriori, i.e. mode finding using numerical 
##' optimization.
##'     
##' @param modelSpecPars an optional `ProtModelSpecPars` of parameters tuning 
##' the specification of the proteomics model. 
##' If not provided, will use previously set `ModelSpecPars` within `x`
##' 
##' @param algoArgs an optional list of tuning arguments to be provided to the 
##' engine. If not provided, will use previously set `algoArgs`
##' within `x`
##'
##' @return A updated `ProtModel` object.
##'
##' @export
##' 
##' @author Philippe Hauchamps
runModel <- function(x, 
                     modelScript = NULL,
                     providerPackage = c("stan", "native"),
                     runMode = c("sampling","MFVB","MAP"),
                     modelSpecPars = NULL,
                     algoArgs = NULL,
                     initialValues = NULL,
                     fullStorage = TRUE,
                     forceRecompile = FALSE) {
    stopifnot(inherits(x, "ProtModel"))
    
    if(is.null(x@InputData)){
        stop("InputData is NULL : run should be done on a properly initialized ProtModel object!")
    }
    
    if(is.null(modelSpecPars)){
        modelSpecPars = ModelSpecPars(x)
        if(is.null(modelSpecPars)) stop("No Model Specifications provided!")
    }
    
    if(is.null(algoArgs)){
        algoArgs = algoArgs(x)
        if(is.null(algoArgs)) stop("No algorithm parameters provided!")
    }
    
    # there is a little trick to be done to the model specification parameters :
    # the mean and sd priors of the gamma parameters need to be aggregated
    # and if the priors of the gamma terms are provided as a single number, 
    # they also need to be duplicated. In all cases, the model input is a
    # vector of gamma terms means and a vector of gamma terms std devs
    # (for the set of MAR gamma parameters excluding the intercept)
    # TO DO : do the same for alpha priors.
    MGamma <- max(x@InputData$MGamma)
    gTM <- priorGammaTermsMean(modelSpecPars)
    gTSd <- priorGammaTermsSd(modelSpecPars)
    if(MGamma > 2){
        if(length(gTM) == 1){
            gTM <- rep(gTM, MGamma-1)
            priorGammaTermsMean(modelSpecPars) = gTM
        }
        if(length(gTSd) == 1){
            gTSd <- rep(gTSd, MGamma-1)
            priorGammaTermsSd(modelSpecPars) = gTSd
        }
    } 
    
    
    if(providerPackage == "stan"){
        x@ModelFit <- ProtModelFitStan(.modelScript = modelScript,
                                        .runMode = runMode,
                                        .recompileModelEachTime = forceRecompile,
                                        .ModelSpecPars = modelSpecPars,
                                        .algoArgs = algoArgs,
                                        .fullStorage = fullStorage) 
    } else if (providerPackage == "native"){
        if(runMode != "MFVB"){
            msg <- "[native] provider package currently only supports "
            msg <- paste0(msg,"MFVB run mode")
            stop(msg)
        }
        x@ModelFit  <- ProtModelFitMFVB(.ModelSpecPars = modelSpecPars,
                                        .algoArgs = algoArgs,
                                        .fullStorage = fullStorage)
    } else {
        msg <- "Unrecognized providerPackage, should be [stan] or [native]"
        stop(msg)
    }
    
    x@ModelFit <- runProtModelEngine(ModelFit(x), InputData(x), initialValues)
    
    
    x@run <- x@run + 1
    x@lastModelSummary <- list()
    return(x)
}
