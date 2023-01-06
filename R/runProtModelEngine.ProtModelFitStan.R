##' @include ProtModelFitStan.R
NULL


##' @importFrom stringr str_detect str_locate str_locate_all
##' 
##' @export
convertInitValuesIntoSTANList <- function(DF){
  STANParamValueList <- list()
  
  nInitValues <- nrow(DF)
  
  for(i in 1:nInitValues){
    currentInitValueName <- DF[i,1]
    currentInitValue <- DF[i,2]
    # try to find [] in the name end, if found => remove it
    breaksPos <- str_locate_all(currentInitValueName, "\\[")[[1]]
    if(nrow(breaksPos) == 0){
      currentName <- currentInitValueName
    } else {
      lastBreakPos <- breaksPos[nrow(breaksPos),1]
      currentName <- str_sub(currentInitValueName, 
                             start = 1, 
                             end = lastBreakPos - 1)
    }
    if(is.null(STANParamValueList[[currentName]])){
      STANParamValueList[[currentName]] <- c(currentInitValue)
    } else {
      STANParamValueList[[currentName]] <- c(
        STANParamValueList[[currentName]], currentInitValue)
    }
  }
  # now convert all single values into arrays of dim 1
  for(i in seq_along(STANParamValueList)){
    if(length(STANParamValueList[[i]]) == 1){
      STANParamValueList[[i]] <- array(data = STANParamValueList[[i]],
                                       dim = 1)
    }
  }
  STANParamValueList
}

##' @export
runSTANAllProtsAtOnce <- function(x, InputData, compiledSTANModel, 
                                  initialValues){
  # prepare stan input data : 
  # merge stan input data list with additional model specification params
  InputData <- c(InputData, as.list.ProtModelSpecPars(ModelSpecPars(x)))
  
  # prepare stan sampling args :
  # remove algoArgs slots that are not stan specific
  stanSamplingArgs <- algoArgs(x)
  stanSamplingArgs$byProt <- NULL
  
  # prepare initial values
  if(is.null(initialValues)){
    initArg <- "random"
  } else {
    STANInitValueList <- convertInitValuesIntoSTANList(initialValues)
    if(runMode(x) == "sampling"){
      nChains <- algoArgs(x)$chains
      if(is.null(nChains)) nChains <- 4
      initArg <- list()
      
      for(c in 1:nChains){
        initArg[[c]] <- STANInitValueList
      }
    } else {
      initArg <- STANInitValueList
    }
  }
  
  # run stan
  genericProbs <- c(0.001, 0.0025, 0.005, 0.01, 0.02, 0.025, 0.03, 0.04,
                    seq(0.05, 0.95, by = 0.05),
                    0.96, 0.97, 0.975, 0.98, 0.99, 0.995, 0.9975, 0.999)
  if(runMode(x) == "sampling"){
    sampling_args = c(list(object = compiledSTANModel,
                           data = InputData,
                           init = initArg),
                      stanSamplingArgs)
    
    cat("**********************************************************************************\n")
    cat("* CALLING Stan simulation engine... (maybe a good timing for a coffee break :-)) *\n")
    cat("**********************************************************************************\n")
    
    Stanfit <- do.call("sampling", args = sampling_args)
    x@stanResList <- list()
    if(x@fullStorage){
      x@stanResList[[1]] <- list(Stanfit = Stanfit)
    } else {
      x@stanResList[[1]] <- 
        list(StanSummary = summary(Stanfit, probs = genericProbs)$summary)
    }
  } else if (runMode(x) == "MFVB"){
    vb_args = c(list(object = compiledSTANModel,
                     data = InputData,
                     init = initArg))
    
    cat("***************************************\n")
    cat("* CALLING Stan engine in mfvb mode... *\n")
    cat("***************************************\n")
    
    Stanfit <- do.call("vb", args = vb_args)
    x@stanfitList <- list()
    if(x@fullStorage){
      x@stanResList[[1]] <- list(Stanfit = Stanfit)
    } else {
      x@stanResList[[1]] <- 
        list(StanSummary = summary(Stanfit, probs = genericProbs)$summary)
    }
  } else { 
    optimizing_args = c(list(object = compiledSTANModel,
                             data = InputData,
                             init = initArg),
                        stanSamplingArgs)
    
    cat("*********************************************\n")
    cat("* CALLING Stan engine in optimizing mode... *\n")
    cat("*********************************************\n")
    
    StanMAPRes <- do.call("optimizing", args = optimizing_args)
    x@stanResList <- list()
    x@stanResList[[1]] <- list()
    #if(x@fullStorage){
      x@stanResList[[1]] <- list(StanMAPRes = StanMAPRes)
    # } else {
    #   x@stanResList[[1]] <- 
    #     list(StanSummary = summary(StanMAPRes, probs = genericProbs)$summary)
    # }
  }
  
  # in case full storage is FALSE, and mode is either "sampling", or "MFVB"
  # => store dignostics for further use
  if(!x@fullStorage && (runMode(x) == "sampling" || runMode(x) == "MFVB")){
    diag <- list()
    diag$NEff <- check_n_eff(Stanfit)
    diag$Rhat <- check_rhat(Stanfit)
    diag$Div <- check_div(Stanfit)
    diag$TreeDepth <- check_treedepth(Stanfit)
    diag$Energy <- check_energy(Stanfit)
    diag$Elapsed <- get_total_elapsed_time(Stanfit)
    x@stanResList[[1]]$StanDiagnostics <- diag
  }
  x
}

##' @export
runSTANOneProt <- function(protIndex,
                           InputData, 
                           compiledSTANModel,
                           runMode,
                           modelSpecPars,
                           algoArgs, 
                           fullStorage,
                           initArg){
  
  # select the InputData that are part of the selected protein
  OneProtInputData <- list()
  k <- protIndex
  
  if(k == 1){
    yobspos <- 1
    xobspos <- 1
    zobspos <- 1
    tobspos <- 1
    ymispos <- 1
    xmispos <- 1
    zmispos <- 1
    tmispos <- 1
    apos <- 1
    bpos <- 1
    cpos <- 1
    spos <- 1
    conpos <- 1
  } else {
    yobspos <- sum(InputData$NObs[1:(k-1)]) + 1
    xobspos <- sum(InputData$NObs[1:(k-1)] * InputData$MAlpha[1:(k-1)]) + 1
    zobspos <- sum(InputData$NObs[1:(k-1)] * InputData$MBeta[1:(k-1)]) + 1
    tobspos <- sum(InputData$NObs[1:(k-1)] * InputData$MGamma[1:(k-1)]) + 1
    ymispos <- sum(InputData$NMis[1:(k-1)]) + 1
    xmispos <- sum(InputData$NMis[1:(k-1)] * InputData$MAlpha[1:(k-1)]) + 1
    zmispos <- sum(InputData$NMis[1:(k-1)] * InputData$MBeta[1:(k-1)]) + 1
    tmispos <- sum(InputData$NMis[1:(k-1)] * InputData$MGamma[1:(k-1)]) + 1
    
    apos <- sum(InputData$MAlpha[1:(k-1)]) + 1
    bpos <- sum(InputData$MBeta[1:(k-1)]) + 1
    cpos <- sum(InputData$MGamma[1:(k-1)]) + 1
    spos <- (k-1) * InputData$nRdm + 1
    conpos <- sum((InputData$nContrast*InputData$MAlpha)[1:(k-1)]) + 1
  }
  
  OneProtInputData$K <- 1
  OneProtInputData$nRdm <- InputData$nRdm
  OneProtInputData$NObs <- array(InputData$NObs[k], dim=1)
  OneProtInputData$NMis <- array(InputData$NMis[k], dim=1)
  OneProtInputData$MAlpha <- array(InputData$MAlpha[k], dim=1)
  OneProtInputData$MBeta <- array(InputData$MBeta[k], dim=1)
  if(InputData$nRdm > 1){
    OneProtInputData$MBetaPerRdmEffect <- 
      InputData$MBetaPerRdmEffect[spos:(spos+InputData$nRdm-1)]
  } else if(InputData$nRdm == 1){
    OneProtInputData$MBetaPerRdmEffect <- 
      array(InputData$MBetaPerRdmEffect[spos], dim=1)
  } else {
    OneProtInputData$MBetaPerRdmEffect <- InputData$MBetaPerRdmEffect #numeric(0)
  }
  OneProtInputData$MGamma <- array(InputData$MGamma[k], dim=1)
  OneProtInputData$NObsSum <- InputData$NObs[k]
  OneProtInputData$NMisSum <- InputData$NMis[k]
  OneProtInputData$XObsLen <- InputData$NObs[k] * InputData$MAlpha[k]
  OneProtInputData$XMisLen <- InputData$NMis[k] * InputData$MAlpha[k]
  OneProtInputData$ZObsLen <- InputData$NObs[k] * InputData$MBeta[k]
  OneProtInputData$ZMisLen <- InputData$NMis[k] * InputData$MBeta[k]
  OneProtInputData$TObsLen <- InputData$NObs[k] * InputData$MGamma[k]
  OneProtInputData$TMisLen <- InputData$NMis[k] * InputData$MGamma[k]
  OneProtInputData$XObs <- 
    InputData$XObs[xobspos:(xobspos + InputData$NObs[k] * InputData$MAlpha[k]-1)]
  OneProtInputData$XMis <- 
    InputData$XMis[xmispos:(xmispos + InputData$NMis[k] * InputData$MAlpha[k]-1)]
  OneProtInputData$yObs <- 
    InputData$yObs[yobspos:(yobspos + InputData$NObs[k]-1)]
  if(InputData$nRdm >= 1){
    OneProtInputData$ZObs <- 
      InputData$ZObs[zobspos:(zobspos + InputData$NObs[k] * InputData$MBeta[k]-1)]
    OneProtInputData$ZMis <- 
      InputData$ZMis[zmispos:(zmispos + InputData$NMis[k] * InputData$MBeta[k]-1)]
  } else {
    OneProtInputData$ZObs <- InputData$ZObs #numeric(0)
    OneProtInputData$ZMis <- InputData$ZMis #numeric(0)
  }
  
  OneProtInputData$TObs <- 
    InputData$TObs[tobspos:(tobspos + InputData$NObs[k] * InputData$MGamma[k]-1)]
  OneProtInputData$TMis <- 
    InputData$TMis[tmispos:(tmispos + InputData$NMis[k] * InputData$MGamma[k]-1)]
  OneProtInputData$refIntensity <- array(InputData$refIntensity[k], 1)
  
  # contrasts data
  OneProtInputData$nContrast <- array(InputData$nContrast[k], dim=1)
  if(InputData$nContrast[k] >= 1){
    OneProtInputData$contrastCoefs <- 
      InputData$contrastCoefs[conpos:(conpos + InputData$nContrast[k] * InputData$MAlpha[k]-1)]
  } else {
    OneProtInputData$contrastCoefs <- numeric(0)
  }
  OneProtInputData$contrastCoefLen <- InputData$nContrast[k] * InputData$MAlpha[k]
  
  # prepare stan input data : 
  # merge stan input data list with additional model specification params
  OneProtInputData <- c(OneProtInputData, as.list.ProtModelSpecPars(modelSpecPars))
  
  # prepare stan sampling args :
  # remove algoArgs slots that are not stan specific
  stanSamplingArgs <- algoArgs
  stanSamplingArgs$byProt <- NULL
  # as we run by prot, we replace cores in algoArgs(x) by one
  # 1 will then be the number of cores used by STAN engine itself
  stanSamplingArgs$cores <- 1
  
  # extract what is needed from init values if they are provided in initArg
  if(is.list(initArg)){
    # init values are provided => extract for the current protein
    if(runMode == "sampling"){
      initialValues <- initArg[[1]]
    } else {
      initialValues <- initArg
    }
    
    # now extract initial values that belong to the current protein
    protInitialValues <- list()
    
    alphas <- initialValues[["alpha"]]
    if(!is.null(alphas)){
      if(InputData$MAlpha[k] > 1){
        protInitialValues[["alpha"]] <- 
          alphas[apos:(apos+InputData$MAlpha[k]-1)]
      } else {
        protInitialValues[["alpha"]] <- 
          array(alphas[apos], dim=1)
      }
    }
    betas <- initialValues[["beta"]]
    if(!is.null(betas)){
      if(InputData$MBeta[k] > 1){
        protInitialValues[["beta"]] <- 
          betas[bpos:(bpos+InputData$MBeta[k]-1)]
      } else {
        protInitialValues[["beta"]] <- 
          array(betas[bpos], dim=1)
      }
    }
    sigma2Eps <- initialValues[["sigma2Eps"]]
    if(!is.null(sigma2Eps)){
      protInitialValues[["sigma2Eps"]] <- array(sigma2Eps[k], dim=1)
    }
    sigma2Rdm <- initialValues[["sigma2Rdm"]]
    if(!is.null(sigma2Rdm)){
      if(InputData$nRdm > 1)
      {
        protInitialValues[["sigma2Rdm"]] <- 
          sigma2Rdm[spos:(spos+InputData$nRdm-1)]
      } else {
        protInitialValues[["sigma2Rdm"]] <- 
          array(sigma2Rdm[spos], dim=1)
      }
    }
    gammas <- initialValues[["gamma"]]
    if(!is.null(gammas)){
      if(InputData$MGamma[k]>1){
        protInitialValues[["gamma"]] <- 
          gammas[cpos:(cpos+InputData$MGamma[k]-1)]
      } else {
        protInitialValues[["gamma"]] <- 
          array(gammas[cpos], dim=1)
      }
    }
    gammYs <- initialValues[["gammY"]]
    if(!is.null(gammYs)){
      protInitialValues[["gammY"]] <- 
        array(gammYs[k], dim=1)
    }
    yMisses <- initialValues[["yMis"]]
    if(!is.null(yMisses)){
      if(InputData$NMis[k] > 1){
        protInitialValues[["yMis"]] <-
          yMisses[ymispos:(ymispos+InputData$NMis[k]-1)]
      } else {
        protInitialValues[["yMis"]] <-
          array(yMisses[ymispos], dim=1)
      }
      
    }
    
    # reinject new initial values in initArg
    if(runMode == "sampling"){
      for(l in seq_along(initArg)){
        initArg[[l]] <- protInitialValues
      }
    } else {
      initArg <- protInitialValues
    }
  }
  
  # run stan
  genericProbs <- c(0.001, 0.0025, 0.005, 0.01, 0.02, 0.025, 0.03, 0.04,
                    seq(0.05, 0.95, by = 0.05),
                    0.96, 0.97, 0.975, 0.98, 0.99, 0.995, 0.9975, 0.999)
  if(runMode == "sampling"){
    sampling_args = c(list(object = compiledSTANModel,
                           data = OneProtInputData,
                           init = initArg),
                      stanSamplingArgs)
    
    cat("******************************************************\n")
    cat(paste0("* CALLING Stan simulation engine for protein ", k, " ...   *\n"))
    cat("******************************************************\n")
    
    Stanfit <- do.call("sampling", args = sampling_args)
    if(fullStorage){
      res <- list(Stanfit = Stanfit)
    } else {
      res <- list(StanSummary = summary(Stanfit, probs = genericProbs)$summary)
    }
    
  } else if (runMode == "MFVB"){
    vb_args = c(list(object = compiledSTANModel,
                     data = OneProtInputData,
                     init = initArg))
    
    cat("******************************************************\n")
    cat(paste0("* CALLING Stan engine in mfvb mode for protein ", k, " ... *\n"))
    cat("******************************************************\n")
    
    Stanfit <- do.call("vb", args = vb_args)
    if(fullStorage){
      res <- list(Stanfit = Stanfit)
    } else {
      res <- list(StanSummary = summary(Stanfit, probs = genericProbs)$summary)
    }
  } else { 
    optimizing_args = c(list(object = compiledSTANModel,
                             data = OneProtInputData,
                             init = initArg),
                        stanSamplingArgs)
    
    cat("************************************************************\n")
    cat(paste0("* CALLING Stan engine in optimizing mode for protein ", k, " ... *\n"))
    cat("************************************************************\n")
    
    stanMAPRes <- do.call("optimizing", args = optimizing_args)
    #if(fullStorage){
      res <- list(StanMAPRes = stanMAPRes)
    # } else {
    #   res <- list(StanSummary = summary(Stanfit, probs = genericProbs)$summary)
    # }
  }
  
  # in case full storage is FALSE, and mode is either "sampling", or "MFVB"
  # => store dignostics for further use
  if(!fullStorage && (runMode == "sampling" || runMode == "MFVB")){
    diag <- list()
    diag$NEff <- check_n_eff(Stanfit)
    diag$Rhat <- check_rhat(Stanfit)
    diag$Div <- check_div(Stanfit)
    diag$TreeDepth <- check_treedepth(Stanfit)
    diag$Energy <- check_energy(Stanfit)
    diag$Elapsed <- get_total_elapsed_time(Stanfit)
    res$StanDiagnostics <- diag
  }
  res
}



##' @export
runSTANByProt <- function(x, 
                          InputData, 
                          compiledSTANModel, 
                          initialValues){
  cores <- algoArgs(x)$cores
  if(is.null(cores)){
    cores = 1
  }
  
  nProt <- InputData$K
  protIndexes <- c(1:nProt)
  #protIndexes <- c(269:280)
  
  # prepare initial values for all proteins
  if(is.null(initialValues)){
    initArg <- "random"
  } else {
    STANInitValueList <- convertInitValuesIntoSTANList(initialValues)
    if(runMode(x) == "sampling"){
      nChains <- algoArgs(x)$chains
      if(is.null(nChains)) nChains <- 4
      initArg <- list()
      
      for(c in 1:nChains){
        initArg[[c]] <- STANInitValueList
      }
    } else {
      initArg <- STANInitValueList
    }
  }
  
  if(cores > 1 && nProt > 1){
    # take registered bpparam if suitable for parallel computation
    # if not, create SnowParam by default as suitable for any platform
    used_bp <- bpparam() 
    if(is.null(used_bp)){# || inherits(used_bp, "SerialParam"))
      used_bp <- SnowParam(workers = cores,
                           tasks = nProt,
                           progressbar = TRUE,
                           type = "SOCK",
                           #RNGseed = my.seed,
                           timeout = 30L * 24L * 60L * 60L,
                           exportglobals = TRUE,
                           log = FALSE,
                           logdir = NA_character_,
                           jobname = "STAN")
    }
    
    bpstopOnError(used_bp) <- FALSE

    cat("Launching all protein STAN models in parallel...\n")
    res <- bptry({bpmapply(FUN = runSTANOneProt, 
                    protIndex = protIndexes,
                    MoreArgs = list(InputData = InputData,
                                    compiledSTANModel = compiledSTANModel,
                                    runMode = runMode(x),
                                    modelSpecPars = ModelSpecPars(x),
                                    algoArgs = algoArgs(x),
                                    fullStorage = fullStorage(x),
                                    initArg = initArg),
                    SIMPLIFY = FALSE,
                    BPPARAM = used_bp)})
    resOK <- bpok(res)
    failedJobs <- which(!resOK)
    if(length(failedJobs) > 0){
      cat("Found failed tasks : ", failedJobs, "\n")
      cat("Relaunching failed tasks sequentially...\n")
      bpmapply(FUN = runSTANOneProt, BPREDO = res, BPPARAM = SerialParam(),
               protIndex = protIndexes,
               MoreArgs = list(InputData = InputData,
                               compiledSTANModel = compiledSTANModel,
                               runMode = runMode(x),
                               modelSpecPars = ModelSpecPars(x),
                               algoArgs = algoArgs(x),
                               fullStorage = fullStorage(x),
                               initArg = initArg),
               SIMPLIFY = FALSE)
      
      cat("Relaunching failed tasks sequentially...DONE!\n")
    }
    
  } else {
    res <- mapply(FUN = runSTANOneProt, 
                  protIndex = protIndexes,
                  MoreArgs = list(InputData = InputData,
                                  compiledSTANModel = compiledSTANModel,
                                  runMode = runMode(x),
                                  modelSpecPars = ModelSpecPars(x),
                                  algoArgs = algoArgs(x),
                                  fullStorage = fullStorage(x),
                                  initArg = initArg),
                  SIMPLIFY = FALSE)
  }
  
  # store everything back into object x
  
  # if(runMode(x) == "sampling" || runMode(x) == "MFVB"){
  #   x@stanfitList <- res
  # } else {
  #   x@stanMAPResList <- res
  # }
  x@stanResList <- res
  
  x@protIndexList <- list()
  x@protIndexList$MAlpha <- InputData$MAlpha
  x@protIndexList$MBeta <- InputData$MBeta
  x@protIndexList$MGamma <- InputData$MGamma
  x@protIndexList$nRdm <- InputData$nRdm
  x@protIndexList$NObs <- InputData$NObs
  x@protIndexList$NMis <- InputData$NMis
  x@protIndexList$nContrast <- InputData$nContrast
  
  x
  
}


##' @export
runProtModelEngine.ProtModelFitStan <- function(x, 
                                                InputData, 
                                                initialValues){
  stopifnot(inherits(x, "ProtModelFitStan"))
  
  if(is.null(modelScript(x))){
    stop("stan needs a modelScript")
  } 
  
  if(is.null(ModelSpecPars(x))){
    stop("Model Specs not initialized!")
  }
  
  if(is.null(algoArgs(x))){
    stop("Algorithm arguments not initialized!")
  }
  
  ## 0. check that `modelScript.stan` exists
  stanScriptFile <- paste0(x@modelScript, ".stan")
  if(!file.exists(stanScriptFile))
    stop(paste0(stanScriptFile, " does not exist!"))
  ## 1. check if `modelScript.rds` exists.
  ## 2. if not, compile it. Then save it as rds.
  ## 3. if `modelScript.rds` exists, make sure it is more recent 
  #     than `modelScript.stan`.
  ## 4  if more recent, load it, otherwise execute step 2
  stanModelFile <- paste0(x@modelScript, ".rds")
  compile <- TRUE
  if (file.exists(stanModelFile) && !recompileModelEachTime(x)){
    fileTimes <- file.mtime(c(stanScriptFile, stanModelFile))
    if(!(fileTimes[2] > fileTimes[1])){
      compile <- FALSE
    }
  }
  
  if(compile)
  {
    cat(paste0("Compiling Stan script : ", stanScriptFile, "\n"))
    stanc_ret <- stanc(file = stanScriptFile, verbose = TRUE)
    
    stan_mod <- stan_model(stanc_ret = stanc_ret,
                           verbose = TRUE,
                           auto_write = TRUE)
    cat("Model compilation successful! Wrighting model on disk...\n")
    saveRDS(object = stan_mod, file = stanModelFile)
    cat("Done!\n")
  } else {
    cat(paste0("Found an updated Stan model : ", stanModelFile, "\n"))
    cat("Uploading...")
    stan_mod <- readRDS(file = stanModelFile)
    cat("Done!\n")
  }
  
  byProt <- algoArgs(x)$byProt
  if(is.null(byProt)){
    byProt = FALSE
  }
  
  hierarchicalModel <- sigmaHierarch(ModelSpecPars(x)) || 
    sigmaRdmHierarch(ModelSpecPars(x)) ||
    gammaHierarch(ModelSpecPars(x))
  
  if(byProt && hierarchicalModel){
    msg <- "detected inconsistency in AlgoArgs : byProt cannot be [TRUE] when "
    msg <- paste0(msg, "model is hierarchical")
    stop(msg)
  }
  
  if(byProt && InputData$K > 1){
    x <- runSTANByProt(x = x, 
                       InputData = InputData,
                       initialValues = initialValues,
                       compiledSTANModel = stan_mod)
  } else {
    x <- runSTANAllProtsAtOnce(x = x, 
                               InputData = InputData,
                               initialValues = initialValues,
                               compiledSTANModel = stan_mod)
  }
  
  x
}

setMethod("runProtModelEngine", "ProtModelFitStan",
          runProtModelEngine.ProtModelFitStan)