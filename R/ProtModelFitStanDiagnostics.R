##' @include ProtModelFit.R
NULL

##' @include stanfitDiagnostics.R
NULL

##' @title Check transitions that ended with a divergence, by model
##' 
##' @param x a `ProtModelFitStan` object
##' 
##' @export
checkDiv <- function(x, displayMessages = TRUE) {
  stopifnot(inherits(x, "ProtModelFitStan"))
  nModels <- length(x@stanResList)
  isProb <- rep(FALSE, nModels)
  for(k in seq_along(x@stanResList)){
    if(x@fullStorage){
      res <- check_div(x@stanResList[[k]]$Stanfit)
    } else {
      res <- x@stanResList[[k]]$StanDiagnostics$Div
    }
    isProb[k] <- res$isProb
    if(displayMessages){
      if(nModels > 1){
        cat("For model ", k, ":\n")
        cat("******************\n")
      }
      cat(res$msg, "\n")
    }
  } 
  isProb
}

##' @title Check transitions that ended prematurely due to maximum tree depth 
##' limit, by model
##' 
##' @param x a `ProtModelFitStan` object
##' 
##' @export
checkTreeDepth <- function(x, displayMessages = TRUE) {
  stopifnot(inherits(x, "ProtModelFitStan"))
  nModels <- length(x@stanResList)
  isProb <- rep(FALSE, nModels)
  for(k in seq_along(x@stanResList)){
    if(x@fullStorage){
      res <- check_treedepth(x@stanResList[[k]]$Stanfit)
    } else {
      res <- x@stanResList[[k]]$StanDiagnostics$TreeDepth
    }
    isProb[k] <- res$isProb
    if(displayMessages){
      if(nModels > 1){
        cat("For model ", k, ":\n")
        cat("******************\n")
      }
      cat(res$msg, "\n")
    }
  } 
  isProb
}

##' @title Check energy Bayesian fraction of missing information (E-BFMI), by
##' model
##' 
##' @param x a `ProtModelFitStan` object
##' 
##' @importFrom stats var
##' 
##' @export
checkEnergy <- function(x, displayMessages = TRUE) {
  stopifnot(inherits(x, "ProtModelFitStan"))
  nModels <- length(x@stanResList)
  isProb <- rep(FALSE, nModels)
  for(k in seq_along(x@stanResList)){
    if(x@fullStorage){
      res <- check_energy(x@stanResList[[k]]$Stanfit)
    } else {
      res <- x@stanResList[[k]]$StanDiagnostics$Energy
    }
    isProb[k] <- res$isProb
    if(displayMessages){
      if(nModels > 1){
        cat("For model ", k, ":\n")
        cat("******************\n")
      }
      cat(res$msg, "\n")
    }
  } 
  isProb
}


##' @title Check effective sample size per iteration, by model
##' 
##' @param x a `ProtModelFitStan` object
##' 
##' @export
checkNEff <- function(x, displayMessages = TRUE) {
  stopifnot(inherits(x, "ProtModelFitStan"))
  nModels <- length(x@stanResList)
  isProb <- rep(FALSE, nModels)
  for(k in seq_along(x@stanResList)){
    if(x@fullStorage){
      res <- check_n_eff(x@stanResList[[k]]$Stanfit)
    } else {
      res <- x@stanResList[[k]]$StanDiagnostics$NEff
    }
    isProb[k] <- res$isProb
    if(displayMessages){
      if(nModels > 1){
        cat("For model ", k, ":\n")
        cat("******************\n")
      }
      cat(res$msg, "\n")
    }
  } 
  isProb
}

##' @title Checks potential scale reduction factors, by model
##' 
##' @param x a `ProtModelFitStan` object
##' 
##' @export
checkRhat <- function(x, displayMessages = TRUE) {
  stopifnot(inherits(x, "ProtModelFitStan"))
  nModels <- length(x@stanResList)
  isProb <- rep(FALSE, nModels)
  for(k in seq_along(x@stanResList)){
    if(x@fullStorage){
      res <- check_rhat(x@stanResList[[k]]$Stanfit)
    } else {
      res <- x@stanResList[[k]]$StanDiagnostics$Rhat
    }
    isProb[k] <- res$isProb
    if(displayMessages){
      if(nModels > 1){
        cat("For model ", k, ":\n")
        cat("******************\n")
      }
      cat(res$msg, "\n")
    }
  } 
  isProb
}


##' @title Run all diagnostics one at a time, by model
##' 
##' @param x a `ProtModelFitStan` object
##' 
##' @export
checkAllDiagnostics <- function(x, displayMessages = TRUE) {
  stopifnot(inherits(x, "ProtModelFitStan"))
  nModels <- length(x@stanResList)
  isProb <- rep(FALSE, nModels)
  for(k in seq_along(x@stanResList)){
    if(x@fullStorage){
      res <- check_all_diagnostics(x@stanResList[[k]]$Stanfit)
      currentIsProb <- res$isProb
      msg <- res$msg
    } else {
      resNEff <- x@stanResList[[k]]$StanDiagnostics$NEff
      resRhat <- x@stanResList[[k]]$StanDiagnostics$Rhat
      resDiv <- x@stanResList[[k]]$StanDiagnostics$Div
      resTreeDepth <- x@stanResList[[k]]$StanDiagnostics$TreeDepth
      resEnergy <- x@stanResList[[k]]$StanDiagnostics$Energy
      
      currentIsProb <- resNEff$isProb ||
        resRhat$isProb ||
        resDiv$isProb ||
        resTreeDepth$isProb ||
        resEnergy$isProb
      
      msg <- paste0(resNEff$msg, "\n",
                    resRhat$msg, "\n",
                    resDiv$msg, "\n",
                    resTreeDepth$msg, "\n",
                    resEnergy$msg)
      
    }
    isProb[k] <- currentIsProb
    if(displayMessages){
      if(nModels > 1){
        cat("For model ", k, ":\n")
        cat("******************\n")
      }
      cat(msg,"\n")
    }
  } 
  isProb
}



##' @title Returns total elapsed time that was needed to run a ProtModelFitStan engine
##'
##' @param x a `ProtModelFitStan` object
##'
##' @export
getElapsedTime <- function(x) {
  
  stopifnot(inherits(x, "ProtModelFitStan"))
  nModels <- length(x@stanResList)
  elapsed <- rep(0., nModels)
  for(k in seq_along(x@stanResList)){
    if(x@fullStorage){
      elapsed[k] <- get_total_elapsed_time(x@stanResList[[k]]$Stanfit)
    } else {
      elapsed[k] <- x@stanResList[[k]]$StanDiagnostics$Elapsed
    }
  } 
  elapsed
}




