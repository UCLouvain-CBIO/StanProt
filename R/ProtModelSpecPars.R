##' @title ProtModelSpecPars class
##'
##' @aliases ProtModelSpecPars-class, ProtModelSpecPars
##'
##' @name ProtModelSpecPars-class
##'
##' @rdname ProtModelSpecPars
##'
##' @description
##'
##' This class encapsulates a set of model specifications for a ProtModelFit.
##'
##' Object can be created manually with `ProtModelSpecPars()`.
##'
##' @slot priorIntensityMean The mean of the (base) intensity prior distribution.
##'
##' @slot priorFoldChangeMean The mean of the log2Fold changes prior distribution.
##' 
##' @slot priorIntensitySd The std dev of the (base) intensity prior distribution.
##' 
##' @slot priorFoldChangeSd The std dev of the log2Fold changes prior distribution.
##' 
##' @slot priorGammaInterceptMean The prior mean of the intercept coefficient
##' in the missingness model
##' 
##' @slot priorGammaTermsMean The prior mean of the other MAR coefficients
##' in the missingness model. Either there is only one mean, and then it is 
##' interpreted as a common prior mean for all MAR gamma terms (except intercept), 
##' or it is a vector. 
##' 
##' @slot priorGammaInterceptSd The prior std dev of the intercept coefficient
##' in the missingness model
##' 
##' @slot priorGammaTermsSd The prior std dev of the other MAR coefficients
##' in the missingness model. Either there is only one mean, and then it is 
##' interpreted as a common prior mean for all MAR gamma terms (except intercept), 
##' or it is a vector. 
##' 
##' @slot priorGammaYMean The prior mean of the NMAR gamma coefficients in the
##' missingness model
##' 
##' @slot priorGammaYSd the prior std dev of the NMAR gamma coefficient in the
##' missingness model
##' 
##' @slot priorNuEpsilon the Nu coefficient of the prior inv chi2 distribution 
##' for the intensity model errors variances (used in non hierarchical set-up 
##' only)
##' 
##' @slot priorScaleEpsilon the Scale coefficient of the prior inv chi2 
##' distribution for the intensity model errors variances (used in non 
##' hierarchical set-up only)
##' 
##' @slot priorNuRdm the Nu coefficient of the prior inv chi2 distribution for
##' the intensity model random effects variances (used in non hierarchical 
##' set-up only)
##' 
##' @slot priorScaleRdm the Scale coefficient of the prior inv chi2 distribution 
##' for the intensity model random effects variances (used in non hierarchical 
##' set-up only)
##' 
##' @slot t_likelihood A `logical` indicating whether we use a student-t 
##' distribution for the peptide intensities (`TRUE`) or normal (`FALSE`)
##' 
##' @slot t_fixedDF A `numeric` designing the degrees of freedom of the student
##' intensity data distribution (if used)
##' 
##' @slot modelMiss A `logical` indicating whether we model missingness or not.
##' 
##' @slot NMAR A `logical` indicating whether we model a NMAR missingness
##' component or not.
##' 
##' @slot sigmaHierarch A `logical` indicating whether we use a hierarchical 
##' set up for the error residual variances, or not.
##' 
##' @slot sigmaRdmHierarch A `logical` indicating whether, for the random 
##' effects (if any), we use a hierarchical set-up for the random effect
##' variances, or not.
##' 
##' @slot gammaHierarch A `logical` indicating whether we use a hierarchical 
##' set up for the missingness model coefficients (gamma & gammY), or not.
##' 
##'
##' @importFrom methods new
##' 
##' @importFrom stats qnorm
##'
##' @exportClass ProtModelSpecPars
##'
##' @examples
##' ProtModelSpecPars()
##' 

.ProtModelSpecPars <-
  setClass("ProtModelSpecPars",
           slots = c(
             priorIntensityMean = "numeric",
             priorFoldChangeMean = "numeric",
             priorIntensitySd = "numeric",
             priorFoldChangeSd = "numeric",
             priorGammaInterceptMean = "numeric",
             priorGammaTermsMean = "numeric",
             priorGammaInterceptSd = "numeric",
             priorGammaTermsSd = "numeric",
             priorGammaYMean = "numeric",
             priorGammaYSd = "numeric",
             priorNuEpsilon = "numeric",
             priorScaleEpsilon = "numeric",
             priorNuRdm = "numeric",
             priorScaleRdm = "numeric",
             t_likelihood = "logical",
             t_fixedDF = "numeric",
             modelMiss = "logical",
             NMAR = "logical",
             sigmaHierarch = "logical",
             sigmaRdmHierarch = "logical",
             gammaHierarch = "logical"))

##' @export
ProtModelSpecPars <- function( 
  .priorIntensityLow = 0,
  .priorIntensityHigh = 30,
  .priorFoldChangeLow = -10,
  .priorFoldChangeHigh = 10,
  .priorIntensityMean = 0.5*(.priorIntensityLow+.priorIntensityHigh),
  .priorFoldChangeMean = 0.5*(.priorFoldChangeLow+.priorFoldChangeHigh),
  .priorIntensitySd = 0.5*(.priorIntensityHigh-.priorIntensityLow)/qnorm(0.995),
  .priorFoldChangeSd = 0.5*(.priorFoldChangeHigh-.priorFoldChangeLow)/qnorm(0.995),
  .priorGammaInterceptLow = -30,
  .priorGammaInterceptHigh = 30,
  .priorGammaTermsLow = -30,
  .priorGammaTermsHigh = 30,
  .priorGammaYLow = -10,
  .priorGammaYHigh = 10,
  .priorGammaInterceptMean = 0.5*(.priorGammaInterceptLow+.priorGammaInterceptHigh),
  .priorGammaTermsMean = 0.5*(.priorGammaTermsLow+.priorGammaTermsHigh),
  .priorGammaInterceptSd = 0.5*(.priorGammaInterceptHigh-.priorGammaInterceptLow)/qnorm(0.995),
  .priorGammaTermsSd = 0.5*(.priorGammaTermsHigh-.priorGammaTermsLow)/qnorm(0.995),
  .priorGammaYMean = 0.5*(.priorGammaYLow+.priorGammaYHigh),
  .priorGammaYSd = 0.5*(.priorGammaYHigh-.priorGammaYLow)/qnorm(0.995),
  .priorNuEpsilon = 1,
  .priorScaleEpsilon = 1,
  .priorNuRdm = 1,
  .priorScaleRdm = 1,
  .t_likelihood = TRUE,
  .t_fixedDF = 2,
  .modelMiss = TRUE,
  .NMAR = TRUE,
  .sigmaHierarch = FALSE,
  .sigmaRdmHierarch = FALSE,
  .gammaHierarch = FALSE){
  
  
  .ProtModelSpecPars(
    priorIntensityMean = .priorIntensityMean,
    priorFoldChangeMean = .priorFoldChangeMean,
    priorIntensitySd = .priorIntensitySd,
    priorFoldChangeSd = .priorFoldChangeSd,
    priorGammaInterceptMean = .priorGammaInterceptMean,
    priorGammaTermsMean = .priorGammaTermsMean,
    priorGammaInterceptSd = .priorGammaInterceptSd,
    priorGammaTermsSd = .priorGammaTermsSd,
    priorGammaYMean = .priorGammaYMean,
    priorGammaYSd = .priorGammaYSd,
    priorNuEpsilon = .priorNuEpsilon,
    priorScaleEpsilon = .priorScaleEpsilon,
    priorNuRdm = .priorNuRdm,
    priorScaleRdm = .priorScaleRdm,
    t_likelihood = .t_likelihood,
    t_fixedDF = .t_fixedDF,
    modelMiss = .modelMiss,
    NMAR = .NMAR,
    sigmaHierarch = .sigmaHierarch,
    sigmaRdmHierarch = .sigmaRdmHierarch,
    gammaHierarch = .gammaHierarch
  )
}

setMethod("show", "ProtModelSpecPars",
          function(object) {
            cat(" Object of class ProtModelSpecPars:\n")
            cat("  priorIntensityMean:", round(object@priorIntensityMean, 4), "\n")
            cat("  priorFoldChangeMean:", round(object@priorFoldChangeMean, 4), "\n")
            cat("  priorIntensitySd:", round(object@priorIntensitySd, 4), "\n")
            cat("  priorFoldChangeSd:", round(object@priorFoldChangeSd, 4), "\n")
            cat("  priorGammaInterceptMean:", round(object@priorGammaInterceptMean, 4), "\n")
            cat("  priorGammaTermsMean:", round(object@priorGammaTermsMean, 4), "\n")
            cat("  priorGammaInterceptSd:", round(object@priorGammaInterceptSd, 4), "\n")
            cat("  priorGammaTermsSd:", round(object@priorGammaTermsSd, 4), "\n")
            cat("  priorGammaYMean:", round(object@priorGammaYMean, 4), "\n")
            cat("  priorGammaYSd:", round(object@priorGammaYSd, 4), "\n")
            cat("  priorNuEpsilon:", round(object@priorNuEpsilon, 4), "\n")
            cat("  priorScaleEpsilon:", round(object@priorScaleEpsilon, 4), "\n")
            cat("  priorNuRdm:", round(object@priorNuRdm, 4), "\n")
            cat("  priorScaleRdm:", round(object@priorScaleRdm, 4), "\n")
            cat("  t_likelihood: ", object@t_likelihood, "\n")
            cat("  t_fixedDF: ", object@t_fixedDF, "\n")
            cat("  modelMiss: ", object@modelMiss, "\n")
            cat("  NMAR: ", object@NMAR, "\n")
            cat("  sigmaHierarch: ", object@sigmaHierarch, "\n")
            cat("  sigmaRdmHierarch: ", object@sigmaRdmHierarch, "\n")
            cat("  gammaHierarch: ", object@gammaHierarch, "\n")
          }
)

setMethod("as.list", "ProtModelSpecPars",
          function(x) {
            pGTM <- x@priorGammaTermsMean
            if(length(pGTM) == 1){
              pGTM <- array(pGTM, dim=1)
            }
            pGTSd <- x@priorGammaTermsSd
            if(length(pGTSd) == 1){
              pGTSd <- array(pGTSd, dim=1)
            }
            return(list( priorIntensityMean = x@priorIntensityMean,
                         priorFoldChangeMean = x@priorFoldChangeMean,
                         priorIntensitySd = x@priorIntensitySd,
                         priorFoldChangeSd = x@priorFoldChangeSd,
                         priorGammaInterceptMean = x@priorGammaInterceptMean,
                         priorGammaTermsMean = pGTM,
                         priorGammaInterceptSd = x@priorGammaInterceptSd,
                         priorGammaTermsSd = pGTSd,
                         priorGammaYMean = x@priorGammaYMean,
                         priorGammaYSd = x@priorGammaYSd,
                         priorNuEpsilon = x@priorNuEpsilon,
                         priorScaleEpsilon = x@priorScaleEpsilon,
                         priorNuRdm = x@priorNuRdm,
                         priorScaleRdm = x@priorScaleRdm,
                         t_likelihood = x@t_likelihood,
                         t_fixedDF = x@t_fixedDF,
                         modelMiss = x@modelMiss,
                         NMAR = x@NMAR,
                         sigmaHierarch = x@sigmaHierarch,
                         sigmaRdmHierarch = x@sigmaRdmHierarch,
                         gammaHierarch = x@gammaHierarch))
          }
)

##' @export
as.list.ProtModelSpecPars <- function(x) {
  stopifnot(inherits(x, "ProtModelSpecPars"))
  pGTM <- x@priorGammaTermsMean
  if(length(pGTM) == 1){
    pGTM <- array(pGTM, dim=1)
  }
  pGTSd <- x@priorGammaTermsSd
  if(length(pGTSd) == 1){
    pGTSd <- array(pGTSd, dim=1)
  }
  return(list( priorIntensityMean = x@priorIntensityMean,
               priorFoldChangeMean = x@priorFoldChangeMean,
               priorIntensitySd = x@priorIntensitySd,
               priorFoldChangeSd = x@priorFoldChangeSd,
               priorGammaInterceptMean = x@priorGammaInterceptMean,
               priorGammaTermsMean = pGTM,
               priorGammaInterceptSd = x@priorGammaInterceptSd,
               priorGammaTermsSd = pGTSd,
               priorGammaYMean = x@priorGammaYMean,
               priorGammaYSd = x@priorGammaYSd,
               priorNuEpsilon = x@priorNuEpsilon,
               priorScaleEpsilon = x@priorScaleEpsilon,
               priorNuRdm = x@priorNuRdm,
               priorScaleRdm = x@priorScaleRdm,
               t_likelihood = x@t_likelihood,
               t_fixedDF = x@t_fixedDF,
               modelMiss = x@modelMiss,
               NMAR = x@NMAR,
               sigmaHierarch = x@sigmaHierarch,
               sigmaRdmHierarch = x@sigmaRdmHierarch,
               gammaHierarch = x@gammaHierarch))
}

# getters and setters

##' @export
priorIntensityMean <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@priorIntensityMean)
}

##' @export
"priorIntensityMean<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@priorIntensityMean <- value
  return(object)
}

##' @export
priorFoldChangeMean <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@priorFoldChangeMean)
}

##' @export
"priorFoldChangeMean<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@priorFoldChangeMean <- value
  return(object)
}

##' @export
priorIntensitySd <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@priorIntensitySd)
}

##' @export
"priorIntensitySd<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@priorIntensitySd <- value
  return(object)
}

##' @export
priorFoldChangeSd <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@priorFoldChangeSd)
}

##' @export
"priorFoldChangeSd<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@priorFoldChangeSd <- value
  return(object)
}

##' @export
priorGammaInterceptMean <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@priorGammaInterceptMean)
}

##' @export
"priorGammaInterceptMean<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@priorGammaInterceptMean <- value
  return(object)
}

##' @export
priorGammaTermsMean <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@priorGammaTermsMean)
}

##' @export
"priorGammaTermsMean<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@priorGammaTermsMean <- value
  return(object)
}

##' @export
priorGammaInterceptSd <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@priorGammaInterceptSd)
}

##' @export
"priorGammaInterceptSd<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@priorGammaInterceptSd <- value
  return(object)
}

##' @export
priorGammaTermsSd <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@priorGammaTermsSd)
}

##' @export
"priorGammaTermsSd<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@priorGammaTermsSd <- value
  return(object)
}

##' @export
priorGammaYMean <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@priorGammaYMean)
}

##' @export
"priorGammaYMean<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@priorGammaYMean <- value
  return(object)
}

##' @export
priorGammaYSd <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@priorGammaYSd)
}

##' @export
"priorGammaYSd<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@priorGammaYSd <- value
  return(object)
}

##' @export
priorNuEpsilon <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@priorNuEpsilon)
}

##' @export
"priorNuEpsilon<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@priorNuEpsilon <- value
  return(object)
}

##' @export
priorScaleEpsilon <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@priorScaleEpsilon)
}

##' @export
"priorScaleEpsilon<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@priorScaleEpsilon <- value
  return(object)
}

##' @export
priorNuRdm <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@priorNuRdm)
}

##' @export
"priorNuRdm<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@priorNuRdm <- value
  return(object)
}

##' @export
priorScaleRdm <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@priorScaleRdm)
}

##' @export
"priorScaleRdm<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@priorScaleRdm <- value
  return(object)
}

##' @export
t_likelihood <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@t_likelihood)
}

##' @export
"t_likelihood<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@t_likelihood <- value
  return(object)
}

##' @export
t_fixedDF <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@t_fixedDF)
}

##' @export
"t_fixedDF<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@t_fixedDF <- value
  return(object)
}

##' @export
modelMiss <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@modelMiss)
}

##' @export
"modelMiss<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@modelMiss <- value
  return(object)
}

##' @export
NMAR <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@NMAR)
}

##' @export
"NMAR<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@NMAR <- value
  return(object)
}

##' @export
sigmaHierarch <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@sigmaHierarch)
}

##' @export
"sigmaHierarch<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@sigmaHierarch <- value
  return(object)
}

##' @export
sigmaRdmHierarch <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@sigmaRdmHierarch)
}

##' @export
"sigmaRdmHierarch<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@sigmaRdmHierarch <- value
  return(object)
}

##' @export
gammaHierarch <- function(object){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  return(object@gammaHierarch)
}

##' @export
"gammaHierarch<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModelSpecPars"))
  object@gammaHierarch <- value
  return(object)
}