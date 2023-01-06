##' @include ProtModelFit.R
NULL



##' @title ProtModel class
##'
##' @aliases ProtModel-class, ProtModel
##'
##' @name ProtModel-class
##'
##' @rdname ProtModel
##'
##' @description
##'
##' This class encapsulates a model fit for proteomics data analysis.
##'
##' Object can be created manually with `ProtModel()` or using the
##' [makeModel()] function.
##'
##' @slot modelName A `character(1)` containing the name of the model.
##'
##' @slot data The proteomics data as a `SummarizedExperiment`.
##'
##' @slot formula The formula defining the main intensity linear model.
##'
##' @slot formulaM The formula defining the missingness linear model.
##'
##' @slot run `numeric(1)` where 0 denotes an initialised model that
##'     hasn't been run. Any value > 1 indicates that the parameters
##'     have been estimated.
##'
##' @slot InputData A `list()` of input data that will be passed to the model
##' engine upon a run.
##'
##' @slot proteinInfo A `data.frame` of per protein information (nb peptides,
##' nb observations, %missing data, ...)
##'
##' @slot paramNames A `list` of model parameter names by category, generated from
##' the parsing of the formulas
##'
##' @slot ModelFit A `ProtModelFit` object which is created at run time
##' 
##' @slot lastModelSummary A `list` containing the args and the summary results
##' of the last computed summary 
##' 
##' @import SummarizedExperiment
##'
##' @importFrom methods new
##'
##' @importFrom stats formula
##'
##' @exportClass ProtModel
##'
##' @examples
##' ProtModel()
.ProtModel <-
    setClass("ProtModel",
         slots = c(
             modelName = "character",
             data = "SummarizedExperiment",
             formula = "formula",
             formulaM = "formula",
             run = "numeric",
             InputData = "list",
             ModelFit = "ProtModelFit",
             proteinInfo = "data.frame",
             paramNames = "list",
             lastModelSummary = "list"))

##' @param .modelName A `character(1)` containing the name of the model.
##'
##' @param .data A `SummarizedExperiment` object.
##'
##' @param .formula A `formula` defining the main intensity linear model.
##'
##' @param .formulaM A `formula` defining the missingness linear model.
##'
##' @param .InputData A `list` of input data that will be passed to the model
##' engine upon a run.
##'
##' @param .proteinInfo A `data.frame` of per protein information (nb peptides,
##' nb observations, %missing data, ...).
##'
##' @param .paramNames A `list` of model parameter names by category.
##'
##' @return An object of class `ProtModel`
##'
##' @rdname ProtModel
##'
##' @export
ProtModel <- function(.modelName = "ProtModelName",
                      .data = SummarizedExperiment(),
                      .formula = formula(),
                      .formulaM = formula(),
                      .InputData = list(),
                      .proteinInfo = data.frame(),
                      .paramNames = list()) {

  ## Any data preparation needed before instantiating the new object
  .ProtModel(modelName = .modelName,
             data = .data,
             formula = .formula,
             formulaM = .formulaM,
             InputData = .InputData,
             proteinInfo = .proteinInfo,
             paramNames = .paramNames,
             run = 0,
             lastModelSummary = list()
             )
}

setMethod("show", "ProtModel",
          function(object) {
              cat("Message of class ProtModel:\n")
              cat(" modelName:", paste(object@modelName), "\n")
              cat(" Data:", paste(dim(object@data), collapse = " x "), "\n")
              cat(" Formula:", deparse(object@formula), "\n")
              cat(" FormulaM:", deparse(object@formulaM), "\n")
              cat(" InputData:", length(object@InputData), "\n")
              cat(" ProteinInfo:, ", paste(dim(object@data), collapse = " x "),
                  "\n")
              cat(" paramNames:", length(object@paramNames), "\n")
              cat(ifelse(object@run == 0, " Not run", " Run"), "\n")
              if(object@run > 0){
                cat("  ModelFit:\n")
                show(object@ModelFit)
              }
          })

# getters and setters

##' @export
modelName <- function(object){
  stopifnot(inherits(object, "ProtModel"))
  return(object@modelName)
}

##' @export
"modelName<-" <- function(object, value){
  stopifnot(inherits(object, "ProtModel"))
  object@modelName <- value
  return(object)
}

##' @export
InputData <- function(object){
  stopifnot(inherits(object, "ProtModel"))
  return(object@InputData)
}

##' @export
SummarizedExperimentData <- function(object){
  stopifnot(inherits(object, "ProtModel"))
  return(object@data)
}

##' @export
dataAssay <- function(object, proteins = NULL){
  stopifnot(inherits(object, "ProtModel"))
  if(is.null(proteins)){
    data.selected <- object@data
  } else {
    data.selected <- object@data[rowData(object@data)$protein %in% proteins,]
  }
  return(assay(data.selected))
}

##' @title Get data aggregation of per protein information
##'
##' @param object a `ProtModel` object
##'
##' @export
proteinInfo <- function(object, proteins = NULL){
  stopifnot(inherits(object, "ProtModel"))
  if(is.null(proteins)){
    return(object@proteinInfo)
  } else {
    return(object@proteinInfo[object@proteinInfo$protein %in% proteins,])
  }
  
}
##' @export
ModelFit <- function(object){
  stopifnot(inherits(object, "ProtModel"))
  return(object@ModelFit)
}

##' @importFrom stringr str_detect str_locate str_locate_all
##'
##' @export
paramDisplayNames <- function(object, paramInternalNames){
  stopifnot(inherits(object, "ProtModel"))
  newNames <- vector(mode = "character", length = length(paramInternalNames))
  i <- 0
  nProt <- object@InputData$K
  for(p in paramInternalNames){
    i <- i+1
    nameVector <- vector()
    indexes <- vector()
    hyperParam <- FALSE
    contrast <- FALSE
    if(str_detect(p, "alpha")){
      prefix <- "alpha"
      nameVector <- object@paramNames$paramNamesAlpha
      indexes <- object@InputData$MAlpha
    } else if(str_detect(p, "beta")){
      prefix <- "beta"
      nameVector <- object@paramNames$paramNamesBeta
      indexes <- object@InputData$MBeta
    } else if(str_detect(p, "gamma")){
      prefix <- "gamma"
      nameVector <- object@paramNames$paramNamesGamma
      indexes <- object@InputData$MGamma
    } else if(str_detect(p, "mu0gm")){
      prefix <- "mu0gm"
      nameVector <- c(object@paramNames$paramNamesGamma[[1]], "Y")
      indexes <- object@InputData$MGamma[1] + 1
      hyperParam <- TRUE
    } else if(str_detect(p, "tau0gm")){
      prefix <- "tau0gm"
      nameVector <- c(object@paramNames$paramNamesGamma[[1]], "Y")
      indexes <- object@InputData$MGamma[1] + 1
      hyperParam <- TRUE
    } else if(str_detect(p, "sigma20gm")){
      prefix <- "sigma20gm"
      nameVector <- c(object@paramNames$paramNamesGamma[[1]], "Y")
      indexes <- object@InputData$MGamma[1] + 1
      hyperParam <- TRUE
    } else if(str_detect(p, "tauRdm")){
      prefix <- "tauRdm"
      nameVector <- object@paramNames$paramNamesSigma2Rdm
      indexes <- rep(object@InputData$nRdm, nProt)
    } else if(str_detect(p, "sigma2Rdm")){
      prefix <- "sigma2Rdm"
      nameVector <- object@paramNames$paramNamesSigma2Rdm
      indexes <- rep(object@InputData$nRdm, nProt)
    } else if(str_detect(p, "sigmaRdm")){
      prefix <- "sigmaRdm"
      nameVector <- object@paramNames$paramNamesSigma2Rdm
      indexes <- rep(object@InputData$nRdm, nProt)
    } else if(str_detect(p, "s0Eps2")){
      prefix <- "s0Eps2"
      nameVector <- c("resid")
      indexes <- c(1)
      hyperParam <- TRUE
    } else if(str_detect(p, "nu0Eps")){
      prefix <- "nu0Eps"
      nameVector <- c("resid")
      indexes <- c(1)
      hyperParam <- TRUE
    } else if(str_detect(p, "s0Rdm2")){
      prefix <- "s0Rdm2"
      nameVector <- object@paramNames$paramNamesSigma2Rdm
      indexes <- object@InputData$nRdm
      hyperParam <- TRUE
    } else if(str_detect(p, "nu0Rdm")){
      prefix <- "nu0Rdm"
      nameVector <- object@paramNames$paramNamesSigma2Rdm
      indexes <- object@InputData$nRdm
      hyperParam <- TRUE
    } else if(str_detect(p, "yMis")){
      prefix <- "yMis"
      nameVector <- object@paramNames$paramNamesYMis
      indexes <- object@InputData$NMis
    } else if(str_detect(p, "invBObs")){
      prefix <- "invBObs"
      nameVector <- object@paramNames$paramNamesYObs
      indexes <- object@InputData$NObs
    } else if(str_detect(p, "invBMis")){
      prefix <- "invBMis"
      nameVector <- object@paramNames$paramNamesYMis
      indexes <- object@InputData$NMis
    } else if(str_detect(p, "aObs")){
      prefix <- "aObs"
      nameVector <- object@paramNames$paramNamesYObs
      indexes <- object@InputData$NObs
    } else if(str_detect(p, "aMis")){
      prefix <- "aMis"
      nameVector <- object@paramNames$paramNamesYMis
      indexes <- object@InputData$NMis
    } else if(str_detect(p, "contrast")){
      contrast <- TRUE
      prefix <- ""
      nameVector <- object@paramNames$contrastNames
      indexes <- object@InputData$nContrast
    }else {
      # no conversion needed
      newNames[i] <- p
      next
    }
    breaksPos <- str_locate_all(p, "\\[")[[1]]
    if(nrow(breaksPos) == 0){
      stop("error : trying to convert internal parameter name which does not end with [xx]")
    }
    lastBreakPos <- breaksPos[nrow(breaksPos),1]

    searchIndexStr <- str_sub(p, start = lastBreakPos+1, end = -2)
    if(searchIndexStr == ""){
      stop("error : trying to convert internal parameter name which does not end with [xx]")
    }
    searchIndex<- as.integer(searchIndexStr)
    if(searchIndex > sum(indexes) || searchIndex < 1){
      stop("error : trying to convert internal parameter name with wrong index")
    }
    
    if(hyperParam){
      newNames[i] <- paste0(prefix, "_", nameVector[searchIndex])
    } else {
      csindexes <- cumsum(indexes)
      prot <- which(csindexes>=searchIndex)[1]
      newNameIndex <- searchIndex - ifelse(prot>1, csindexes[prot-1], 0)
      
      if(contrast){
        newNames[i] <- paste0(nameVector[[prot]][newNameIndex], "[",
                              prot, "]")
      } else {
        newNames[i] <- paste0(prefix, "_",
                              nameVector[[prot]][newNameIndex], "[",
                              prot, "]")
      }
    }
    
  }
  newNames
}

##'
##' @export
##'
##' @importFrom stringr str_detect str_locate str_locate_all
paramInternalNames <- function(object, paramDisplayNames){
  stopifnot(inherits(object, "ProtModel"))
  newNames <- vector(mode = "character", length = length(paramInternalNames))
  i <- 0
  nProt <- object@InputData$K
  for(p in paramDisplayNames){
    i <- i+1
    nameVector <- vector()
    indexes <- vector()
    hyperParam <- FALSE
    contrast <- FALSE
    if(str_detect(p, "alpha")){
      prefix <- "alpha"
      nameVector <- object@paramNames$paramNamesAlpha
      indexes <- object@InputData$MAlpha
    } else if(str_detect(p, "beta")){
      prefix <- "beta"
      nameVector <- object@paramNames$paramNamesBeta
      indexes <- object@InputData$MBeta
    } else if(str_detect(p, "gamma")){
      prefix <- "gamma"
      nameVector <- object@paramNames$paramNamesGamma
      indexes <- object@InputData$MGamma
    } else if(str_detect(p, "mu0gm")){
      prefix <- "mu0gm"
      nameVector <- c(object@paramNames$paramNamesGamma[[1]], "Y")
      indexes <- object@InputData$MGamma[1] + 1
      hyperParam <- TRUE
    } else if(str_detect(p, "tau0gm")){
      prefix <- "tau0gm"
      nameVector <- c(object@paramNames$paramNamesGamma[[1]], "Y")
      indexes <- object@InputData$MGamma[1] + 1
      hyperParam <- TRUE
    } else if(str_detect(p, "sigma20gm")){
      prefix <- "sigma20gm"
      nameVector <- c(object@paramNames$paramNamesGamma[[1]], "Y")
      indexes <- object@InputData$MGamma[1] + 1
      hyperParam <- TRUE
    } else if(str_detect(p, "tauRdm")){
      prefix <- "tauRdm"
      nameVector <- object@paramNames$paramNamesSigma2Rdm
      indexes <- rep(object@InputData$nRdm, nProt)
    } else if(str_detect(p, "sigma2Rdm")){
      prefix <- "sigma2Rdm"
      nameVector <- object@paramNames$paramNamesSigma2Rdm
      indexes <- rep(object@InputData$nRdm, nProt)
    } else if(str_detect(p, "sigmaRdm")){
      prefix <- "sigmaRdm"
      nameVector <- object@paramNames$paramNamesSigma2Rdm
      indexes <- rep(object@InputData$nRdm, nProt)
    } else if(str_detect(p, "s0Eps2")){
      prefix <- "s0Eps2"
      nameVector <- c("resid")
      indexes <- c(1)
      hyperParam <- TRUE
    } else if(str_detect(p, "nu0Eps")){
      prefix <- "nu0Eps"
      nameVector <- c("resid")
      indexes <- c(1)
      hyperParam <- TRUE
    } else if(str_detect(p, "s0Rdm2")){
      prefix <- "s0Rdm2"
      nameVector <- object@paramNames$paramNamesSigma2Rdm
      indexes <- object@InputData$nRdm
      hyperParam <- TRUE
    } else if(str_detect(p, "nu0Rdm")){
      prefix <- "nu0Rdm"
      nameVector <- object@paramNames$paramNamesSigma2Rdm
      indexes <- object@InputData$nRdm
      hyperParam <- TRUE
    } else if(str_detect(p, "yMis")){
      prefix <- "yMis"
      nameVector <- object@paramNames$paramNamesYMis
      indexes <- object@InputData$NMis
    } else if(str_detect(p, "invBObs")){
      prefix <- "invBObs"
      nameVector <- object@paramNames$paramNamesYObs
      indexes <- object@InputData$NObs
    } else if(str_detect(p, "invBMis")){
      prefix <- "invBMis"
      nameVector <- object@paramNames$paramNamesYMis
      indexes <- object@InputData$NMis
    } else if(str_detect(p, "aObs")){
      prefix <- "aObs"
      nameVector <- object@paramNames$paramNamesYObs
      indexes <- object@InputData$NObs
    } else if(str_detect(p, "aMis")){
      prefix <- "aMis"
      nameVector <- object@paramNames$paramNamesYMis
      indexes <- object@InputData$NMis
    } else { # check taylor made contrast names
      for(conName in unique(unlist(object@paramNames$contrastNames))){
        if(str_detect(p, conName)){
          contrast = TRUE
          break
        }
      }
      if(contrast){
        prefix <- "contrast"
        nameVector <- object@paramNames$contrastNames
        indexes <- object@InputData$nContrast
      } else {
        # no conversion needed
        newNames[i] <- p
        next
      }
    }
    
    if(contrast){
      firstIndexAfterUnderscore <- 1 # no underscore
    } else {
      firstUnderscore <- str_locate(p, "_")
      if(is.na(firstUnderscore[[1]])){
        # don't know how to make the conversion
        newNames[i] <- p
        next
      }
      firstIndexAfterUnderscore <- firstUnderscore[1,1]+1
    }
    
    if(hyperParam){
      internalName <- str_sub(p, start=firstIndexAfterUnderscore, 
                              end = str_length(p))
      internalNameIndex <- which(nameVector == internalName)[1]
      if(length(internalNameIndex) == 0) next
      newNames[i] <- paste0(prefix, "[",
                            internalNameIndex, "]")
    } else {
      breaksPos <- str_locate_all(p, "\\[")[[1]]
      if(nrow(breaksPos) == 0){
        # don't know how to make the conversion
        newNames[i] <- p
        next
      }
      lastBreakPos <- breaksPos[nrow(breaksPos),1]
      
      protNbStr <- str_sub(p, start = lastBreakPos+1, end = -2)
      if(protNbStr == ""){
        # don't know how to make the conversion
        newNames[i] <- p
        next
      }
      protNb<- as.integer(protNbStr)
      
      internalName <- str_sub(p, start=firstIndexAfterUnderscore, end = lastBreakPos-1)
      internalNameIndex <- which(nameVector[[protNb]] == internalName)[1]
      if(length(internalNameIndex) == 0) next
      csum <- ifelse(protNb>1, sum(indexes[1:(protNb-1)]), 0)
      internalIndex <- csum + internalNameIndex
      
      newNames[i] <- paste0(prefix, "[",
                            internalIndex, "]")
    }

    
  }
  newNames
}

##' @export
computeProtModelSummary<- function(object, ...){
  stopifnot(inherits(object, "ProtModel"))
  if(object@run == 0){
    stop("ProtModel object has not run yet")
  }
  
  args <- list(...)
  if(is.null(args$probs)){
    probs <-  c(0.025,0.25,0.5,0.75,0.975)
  } else {
    probs <- args$probs
  }
  
  args["probs"] <- NULL
  
  summ <- do.call(what = summary, 
                  args = c(list(object = object@ModelFit, probs = probs),
                           args))
  #summ <- summary(object@ModelFit, probs = probs )
  rownames(summ) <-
    paramDisplayNames(object,
                      rownames(summ))
  # update summary storage
  object@lastModelSummary$probs <- probs
  object@lastModelSummary$summ <- summ
  
  return(object)
}

##' @importFrom stringr str_locate_all str_sub str_detect
##' @export
getProtModelSummary <- function(object, ...){
  stopifnot(inherits(object, "ProtModel"))
  if(object@run == 0){
    stop("ProtModel object has not run yet")
  }
  
  if(is.null(object@lastModelSummary)){
    msg <- "model summary has not been extracted yet."
    msg <- c(msg, "\nPls call computeModelSummary() first!")
    stop(msg)
  }
  args <- list(...)
  if(is.null(args$probs)){
    probs <-  object@lastModelSummary$probs
  } else {
    probs <- args$probs
  }
  
  # if(identical(object@lastModelSummary$probs, probs)){
  #   # summary has already been called with the same set of probs => use it
  #   summ <-  object@lastModelSummary$summ
  if(sum(round(probs,12) %in% round(object@lastModelSummary$probs, 12)) == length(probs)){
    # all asked probs have been computed in last summary 
    # => return (subset of) last summary
    summ <- object@lastModelSummary$summ
    summNCol <- ncol(summ)
    summColNames <- colnames(summ)
    keepCol <- rep(TRUE, summNCol)
    j <- 0
    for(k in seq_along(summColNames)){
      if(str_detect(summColNames[k], "%")){
        j <- j+1
        if (!round(object@lastModelSummary$probs[j], 12) %in% round(probs,12)){
          keepCol[k] <- FALSE
        }
      }
    }
    summ <- summ[, keepCol]
  } else {
    msg <- "model summary has not been extracted with all the provided "
    msg <- c(msg, "quantiles argument.\nPls call computeModelSummary() first!")
    stop(msg)
  }
  
  # additional filtering if instructions are provided
  
  pattern <- args$pattern
  if(!is.null(pattern)){
    nR <- nrow(summ)
    keepRow <- rep(FALSE, nR)
    for(pat in pattern){
      patternIndexes <- which(str_detect(rownames(summ), pat))
      keepRow[patternIndexes] <- TRUE
    }
    summ <- summ[keepRow, , drop = FALSE]
  }
  whichProtein <- args$whichProtein
  if(!is.null(whichProtein)){
    nR <- nrow(summ)
    keepRow <- rep(TRUE, nR)
    for(i in 1:nR){
      p <- rownames(summ)[i]
      breaksPos <- str_locate_all(p, "\\[")[[1]]
      if(nrow(breaksPos) == 0){ # protein indicator not found => keep
        next
      }
      lastBreakPos <- breaksPos[nrow(breaksPos),1]
      
      protIndexStr <- str_sub(p, start = lastBreakPos+1, end = -2)
      if(protIndexStr == ""){ # protein indicator not found => keep
        next
      }
      protIndex <- as.integer(protIndexStr)
      if(!protIndex %in% whichProtein){
        keepRow[i] <- FALSE
      }
    }
    summ <- summ[keepRow, , drop = FALSE]
  }
  
  summ
}

#setMethod("summary", "ProtModel", sumProtModel)

##' @export
getFinalValues <- function(x, whichValue = "mean"){
  stopifnot(inherits(x, "ProtModel"))
  if(x@run == 0){
    stop("ProtModel object has not run yet")
  }
  
  res <- getModelFitFinalValues(ModelFit(x), whichValue)
  res
}

##' @export
getParamSample <- function(x, paramName){
  stopifnot(inherits(x, "ProtModel"))
  if(x@run == 0){
    stop("ProtModel object has not run yet")
  }
  
  # internalParamName <- paramInternalNames(x, paramName)
  # sample <- getModelFitParamSample(ModelFit(x), internalParamName)
  
  tryCatch(
    error = function(cnd){
      tryCatch(
        error = function(cnd){
          msg <- paste0(paramName, " not found in object;\n")
          msg <- paste0(msg, "Error : ", conditionMessage(cnd))
          stop(msg)
        },
        {
          internalParamName <- paramInternalNames(x, paramName)
          mySample <<- getModelFitParamSample(ModelFit(x), internalParamName)
        }
      )
    },
    mySample <<- getModelFitParamSample(ModelFit(x), paramName)
  )
  
  mySample
}

##' @export
computeRMSE <- function(x, 
                        variableName, 
                        perProtTargetValues,
                        estimateType = "median"){
  stopifnot(inherits(x, "ProtModel"))
  if(x@run == 0){
    stop("ProtModel object has not run yet")
  }
  mySum <- getProtModelSummary(x, 
                               probs = x@lastModelSummary$probs,
                               pattern = variableName)
  
  if(nrow(mySum) != length(perProtTargetValues)){
    msg <- paste0("summary of protModel, using [",variableName,"] as pattern, ")
    msg <- paste0(msg, "does not return same number of lines than the input ")
    msg <- paste0(msg, "targets")
    stop(msg)
  }
  
  if(estimateType == "mean"){
    obtainedValues <- mySum[,"mean"]
  } else if(estimateType == "median"){
    obtainedValues <- mySum[,"50%"]
  } else {
    stop("estimateType not recognized, allowed values are [mean] or [median]")
  }
  
  errors <- obtainedValues - perProtTargetValues
  
  details <- 
    data.frame(name = proteinInfo(x)$protein,
               target = perProtTargetValues,
               obtained = obtainedValues,
               error = errors,
               sqrError = errors * errors)
  RMSE <- sqrt(mean(details$sqrError, na.rm = TRUE))
  RMSE.sd <- RMSE * sqrt(0.5/sum(!is.na(details$sqrError)))
  RMSE.lb <- RMSE - 1.96 * RMSE.sd
  RMSE.ub <- RMSE + 1.96 * RMSE.sd
  
  res <- list(RMSE = RMSE,
              RMSE.sd = RMSE.sd,
              RMSE.lb = RMSE.lb,
              RMSE.ub = RMSE.ub,
              details = details)
  
  res
}

##' @importFrom stringr str_detect str_locate str_sub
##' 
##' @importFrom Hmisc approxExtrap
##' 
##' @export
computeDAProbabilities <- function(x,
                                   variableName,
                                   lFCThreshold,
                                   testType = c("bilateral", "right", "left")){
                                 
  stopifnot(inherits(x, "ProtModel"))
  if(x@run == 0){
    stop("ProtModel object has not run yet")
  }
  if(lFCThreshold < 0){
    lFCThreshold <- - lFCThreshold
  }
  mySum <- getProtModelSummary(x, 
                               probs = x@lastModelSummary$probs,
                               pattern = variableName)
  
  pctFlags <- str_detect(colnames(mySum),"%")
  pctColInSum <- which(pctFlags)
  nExistingProbs <- sum(pctFlags)
  existingProbs <- rep(0., nExistingProbs)
  k <- 0
  for(col in pctColInSum){
    k <- k+1
    theColName <- colnames(mySum)[col]
    pctCharLoc <- str_locate(theColName, "%")
    if(is.na(pctCharLoc[[1]])){
      stop("Summary col name without [%] = unexpected inconsistency")
    }
    probStr <- str_sub(theColName, 1, pctCharLoc[[1]]-1)
    existingProbs[k] <- as.numeric(probStr)/100
  }
  
  nProts <- nrow(mySum)
  DAProbs <- vector(mode = "numeric", length = nProts)
  for(j in 1:nProts){
    theQuantiles <- 
      approxExtrap( y = existingProbs, 
                    x = unlist(mySum[j, pctColInSum]),
                    xout = c(-lFCThreshold, lFCThreshold),
                    method = "linear")$y
    theQuantiles <- sapply(theQuantiles, 
                           FUN = function(x){
                             max(c(0., min(c(1., x))))
                           }) 
    if(testType[1] == "bilateral"){
      DAProbs[j] <- max(c(
        theQuantiles[1], 1-theQuantiles[2]))
    } else if(testType[1] == "right"){
      DAProbs[j] <- 
        1-theQuantiles[2]
    } else if(testType[1] == "left"){
      DAProbs[j] <- 
        theQuantiles[1]
    } else {
      msg <- "unrecognized testType, possible values are [bilateral], "
      msg <- paste0(msg, "[right] or [left]")
      stop(msg)
    }
  }
  
  DAProbs
}

##' @export
computeROCCurve <- function(x, 
                            variableName, 
                            perProtTargetDAs,
                            lFCThreshold,
                            probThresholds = c(seq(0.5,0.95,0.05),c(0.975,0.98,0.99,0.995,0.999)),
                            testType = c("bilateral","right","left")){
  
  DAProbs <- computeDAProbabilities(x, 
                                    variableName = variableName,
                                    lFCThreshold = lFCThreshold,
                                    testType = testType
                                    )
  
  if(length(DAProbs) != length(perProtTargetDAs)){
    msg <- paste0("summary of protModel, using [",variableName,"] as pattern, ")
    msg <- paste0(msg, "does not return same number of lines than the input targets")
    stop(msg)
  }
  
  nQ <- length(probThresholds)
  nProts <- length(DAProbs)
  details <- data.frame(matrix(nrow = nProts, ncol = nQ, data = NA))
  colnames(details) <- paste0(probThresholds*100,"%")
  
  for(j in 1:nProts){
    if(!is.na(perProtTargetDAs[j])){
      details[j, 1:nQ] <-
        DAProbs[j]>probThresholds
    }
  }
  
  ROC <- data.frame(matrix(nrow = nQ, ncol = 10))
  #colnames(ROC) <- c("TP", "FP", "TN", "FN", "NAs",)
  rownames(ROC) <- colnames(details)
  ROC <- data.frame(t(apply(details,
                            MARGIN = 2,
                            FUN = function(x, target){
                              TP = sum(x & target, na.rm = TRUE)
                              FP = sum(x & !target, na.rm = TRUE)
                              TN = sum(!x & !target, na.rm = TRUE)
                              FN = sum(!x & target, na.rm = TRUE)
                              NAs = sum(is.na(x))
                              sensi = TP/(TP+FN)
                              speci = TN/(TN+FP)
                              FDP = FP/(TP+FP)
                              NAprop = NAs/length(x)
                              c(TP = TP,
                                FP = FP,
                                TN = TN,
                                FN = FN,
                                NAs = NAs,
                                sensi = sensi,
                                speci = speci,
                                FPR = 1-speci,
                                FDP = FDP,
                                NAprop = NAprop)
                            },
                            target = perProtTargetDAs)))
               
  
  res <- list(ROC = ROC,
              details = details)
}

##' @export
computeVolcanoPlotData <- function(x, 
                                   variableName, 
                                   lFCThreshold,
                                   testType = c("bilateral","right","left"),
                                   estimateType = "median"){
  DAProbs <- computeDAProbabilities(x, 
                                    variableName = variableName,
                                    lFCThreshold = lFCThreshold,
                                    testType = testType)
  mySum <- getProtModelSummary(x, 
                               probs = c(0.5),
                               pattern = variableName)
  
  if(estimateType == "mean"){
    obtainedValues <- mySum[,"mean"]
  } else if(estimateType == "median"){
    obtainedValues <- mySum[,"50%"]
  } else {
    stop("estimateType not recognized, allowed values are [mean] or [median]")
  }
  
  details <- data.frame(lFC = obtainedValues,
                        DAProb = DAProbs)
  details
  
}





