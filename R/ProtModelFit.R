##' @include ProtModelSpecPars.R
NULL

##' @title ProtModelFit class
##'
##' @aliases ProtModelFit-class, ProtModelFit
##'
##' @name ProtModelFit-class
##'
##' @rdname ProtModelFit
##'
##' @description
##'
##' This class encapsulates a generic model fit for proteomics data analysis.
##'
##' Object is created internally when runModel() is called on a ProtModel object
##' and stored as a slot in this object
##'
##' @slot ModelSpecPars A `ProtModelSpecPars` object containing STAN model
##' specification parameters.
##'
##' @slot algoArgs A `list()` containing tuning arguments for the ModelFit engine
##' 
##' @slot fullStorage A 'logical' indicating whether full storage of the 
##' parameters data will be done after engine run
##'
##' @import methods
##
##' @importFrom stats formula
##'
##' @exportClass ProtModelFit
##'
.ProtModelFit <-
  setClass("ProtModelFit",
           slots = c(
             ModelSpecPars = "ProtModelSpecPars",
             algoArgs = "list",
             fullStorage = "logical"))

##' 
##' @param .ModelSpecPars A `ProtModelSpecPars` object containing the model
##' specification parameters.
##' 
##' @param .algoArgs A 'list` containing tuning arguments for the ModelFit engine
##' 
##' @param .fullStorage A 'logical' indicating whether full storage of the 
##' parameters data will be done after engine run
##' 
##' @return An object of class `ProtModelFit`
##'
##' @rdname ProtModelFit
##'
ProtModelFit <- function(.ModelSpecPars = ProtModelSpecPars(),
                         .algoArgs = list(),
                         .fullStorage = FALSE) {
  
  ## Any data preparation needed before instantiating the new object
  .ProtModelFit(ModelSpecPars = .ModelSpecPars,
                algoArgs = .algoArgs,
                fullStorage = .fullStorage)
}

setMethod("show", "ProtModelFit",
          function(object) {
            cat("Message of class ProtModelFit:\n")
            cat(" ModelSpecPars:\n")
            show(object@ModelSpecPars)
            cat(" algoArgs:\n")
            show(object@algoArgs)
            cat(" fullStorage:", fullStorage, "\n")
          })

setGeneric("runProtModelEngine", 
           function(x, InputData, initialValues) 
             standardGeneric("runProtModelEngine"))

##' @export
setGeneric("getModelFitParamSample", 
           function(x, paramName) standardGeneric("getModelFitParamSample"))

##' @export
setGeneric("getModelFitFinalValues", 
           function(x, whichValue) standardGeneric("getModelFitFinalValues"))

# getters and setters

##' @export
ModelSpecPars <- function(object){
  stopifnot(inherits(object, "ProtModelFit"))
  return(object@ModelSpecPars)
}

##' @export
"ModelSpecPars<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelFit"))
  stopifnot(inherits(value, "ProtModelSpecPars"))
  object@ModelSpecPars <- value
  return(object)
}

##' @export
algoArgs <- function(object){
  stopifnot(inherits(object, "ProtModelFit"))
  return(object@algoArgs)
}

##' @export
"algoArgs<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelFit"))
  stopifnot(is.list(value))
  object@algoArgs <- value
  return(object)
}

##' @export
fullStorage <- function(object){
  stopifnot(inherits(object, "ProtModelFitStan"))
  return(object@fullStorage)
}



