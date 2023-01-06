##' @include ProtModel.R
NULL

##' This function Updates the model parameter names of a ProtModel object, 
##' organized by category, as a result of the formulas parsing.
##'
##' @title Update parameter names and mappings.
##' 
##' @param x ProtModel object
##' 
##' @return `x` with an updated model parameter names mappings
##'
##' @import dplyr
##' 
##' @importFrom magrittr %>%
##' 
##' @export
##' 
##' @author Philippe Hauchamps
updateParamNames <- function(x){
    
    # sanity checks of user input
    stopifnot(inherits(x, "ProtModel"))
    
    if( length(x@InputData) == 0) 
        stop("Call to updateParamNames needs a proprely initialized ProtModel class")
    
    if( length(x@paramNames) == 0) 
        stop("Call to updateParamNames needs a proprely initialized ProtModel class")
    
    # build vector of beta coefficient names
    alphaNames <- vector(mode = "character")
    betaNames <- vector(mode = "character")
    gammaNames <- vector(mode = "character")
    sigma2RdmNames <- vector(mode = "character")
    tauRdmNames <- vector(mode = "character")
    
    #browser()
    for(i in seq_along(x@paramNames$paramNamesAlpha)){
        alphaNames <- c(alphaNames, 
                        paste0("alpha_", x@paramNames$paramNamesAlpha[[i]], 
                               "[", i, "]"))
    }
    for(i in seq_along(x@paramNames$paramNamesBeta)){
        betaNames <- c(betaNames, 
                       paste0("beta_", x@paramNames$paramNamesBeta[[i]], 
                              "[", i, "]"))
    }
    for(i in seq_along(x@paramNames$paramNamesGamma)){
        gammaNames <- c(gammaNames, 
                        paste0("gamma_", x@paramNames$paramNamesGamma[[i]], 
                               "[", i, "]"))
    }
    for(i in seq_along(x@paramNames$paramNamesSigma2Rdm)){
        sigma2RdmNames <- c(sigma2RdmNames, 
                            paste0("sigma2_", x@paramNames$paramNamesSigma2Rdm[[i]], 
                                   "[", i, "]"))
        tauRdmNames <- c(tauRdmNames,
                         paste0("tau_", x@paramNames$paramNamesSigma2Rdm[[i]], 
                                "[", i, "]"))
    }
    x@paramNames$alphaNames <- alphaNames
    x@paramNames$betaNames <- betaNames
    x@paramNames$gammaNames <- gammaNames
    x@paramNames$sigma2RdmNames <- sigma2RdmNames
    x@paramNames$tauRdmNames <- tauRdmNames
    
    
    # build coefficients names for missing observations
    #browser()
    nProt <- x@InputData$nProt
    yMisNames <- vector(mode = "character")
    for(p in 1:nProt){
        nMis <- x@InputData$NMis[p]
        currentYMisNames <- paste0("yMis[",
                                   p,
                                   "][",
                                   x@paramNames$misRuns[[p]],
                                   "][",
                                   x@paramNames$misPepts[[p]],
                                   "]")
        yMisNames <- c(yMisNames, currentYMisNames)
    }
    x@paramNames$yMisNames <- yMisNames
    

    # here do it only once as s0 and nu0 coefficients are shared between proteins
    nu0RdmNames <- paste0("nu0_", x@paramNames$paramNamesSigma2Rdm[[1]])
    x@paramNames$nu0RdmNames <- nu0RdmNames
    s0RdmNames <- paste0("s0_", x@paramNames$paramNamesSigma2Rdm[[1]])
    x@paramNames$s0RdmNames <- s0RdmNames
    
    x
}



##' This function takes raw data input and builds ProtModel input data.
##'
##' @title Builds ProtModel input data
##' 
##' @param data the Summarized Experiment
##' 
##' @param form the formula used to define the main intensity linear model
##' 
##' @param formM the formula used to define the missingness model
##' 
##' @param protColName A character(1) which designates the name of the column
##' in `data@rowData` that contains the protein identifiers
##' 
##' @param peptColName A character(1) which designates the name of the column
##' in `data@rowData` that contains the peptide identifiers
##' 
##' @param cleanNonIdentifiable TRUE if we clean multicollinearity explicitly
##' 
##' @return A `list()` containing various STAN model input data
##'
##' @import dplyr
##' 
##' @importFrom magrittr %>%
##' 
##' @importFrom stats model.matrix formula as.formula reformulate
##' 
##' @importFrom lme4 findbars nobars
##' 
##' @importFrom stringr str_replace_all
##' 
##' @importFrom utils head
##' 
##' @export
##' 
##' @author Philippe Hauchamps
buildModelInputs <- function( data, 
                              form,
                              formM,
                              protColName, 
                              peptColName,
                              cleanNonIdentifiable = FALSE,
                              contrasts){
    
    #browser()
    protColName <- sym(protColName)
    peptColName <- sym(peptColName)
    
    .rowData <- data %>% rowData %>% as.data.frame
    .colData <- data %>% colData %>% as.data.frame
    .assay <- data %>% assay %>% as.data.frame
    
    dummy <- nPept <- endPeptideIndex <- NULL # CMD CHECK NOTE hack
    nObsPeptide <- nTot <- nObs <- NULL # CMD CHECK NOTE hack
    
    # runs
    nRun <- nrow(.colData)
    if(!("run" %in% colnames(.colData))){
        # add a run column explicitely
        .colData <- cbind(.colData,
                          "run" = factor(paste0("R", 1:nRun)))
        colData(data) <- .colData
    }
    
    
    proteinSet <- .rowData  %>% select(!!protColName) %>% distinct %>% 
        pull(!!protColName)
    
    
    
    nProt <- length(proteinSet)
    
    cat("Building STAN model from SummarizedExperiment with :\n")
    cat(paste0("nPeptides = ", dim(data)[[1]], 
               "; nProteins = ", nProt,
               "; nCols = ", dim(data)[[2]],"\n"))
    cat("Head(Proteins Set) : ")
    print(head(proteinSet))
    cat("..................................................\n")
    
    #browser()
    proteinInfo <- cbind(.rowData, 
                         .assay, 
                         nObsPeptide = rowSums(!is.na(.assay)),
                         sumIntensityPeptide = rowSums(.assay, na.rm = TRUE)) %>% 
        group_by(!!protColName) %>%
        summarize( nPept = n(), 
                   nObs = sum(nObsPeptide), 
                   nTot = nPept * ncol(.assay),
                   pctNA = round((nTot-nObs)*100/nTot,1),
                   avgObsIntensity = sum(sumIntensityPeptide) / sum(nObsPeptide),
                   .groups = 'drop') %>% 
        ungroup %>%
        mutate( endPeptideIndex = cumsum(nPept),
                startPeptideIndex = endPeptideIndex - nPept + 1,
                dummy = 1,
                index = cumsum(dummy)) %>%
        select(-dummy) %>% as.data.frame
    
    
    buildOneModelInputs <- function( data,
                                     form,
                                     formM,
                                     startIndex,
                                     nRows,
                                     cleanNonIdentifiable,
                                     contrasts){

      
        #cat("Building One Model Inputs with startIndex = ", startIndex, "\n")      
        y <- value <- NULL # CMD CHECK NOTE hack

        .assay <- assay(data)[startIndex:(startIndex+nRows-1),] %>% 
            as.data.frame
        .colData <- data %>% colData %>% as.data.frame

        df <- .assay %>% tibble::rownames_to_column(var = "peptide") %>%
            reshape2::melt(id.vars = c("peptide"), variable.name = "run", na.rm = FALSE) %>%
            inner_join(.colData, by = "run") %>% rename(y = value)
        
        formFixed <- lme4::nobars(form)
        
        options(contrasts=c('contr.treatment','contr.poly'))
        
        X <- as.matrix(model.matrix(object = formFixed, data = df))
        
        keepAlpha <- rep(FALSE, ncol(X))
        # apply tolerance slightly lower than default (1e-7) to keep the 
        # slightly-non-multicollinear column
        qr.X <- qr(X, tol=1e-9, LAPACK = FALSE)
        rnkX <- qr.X$rank
        
        #browser()
        if (cleanNonIdentifiable){
            keepAlpha[qr.X$pivot[seq_len(rnkX)]] <- TRUE
            X <- X[, keepAlpha, drop=FALSE]
            MAlpha <- rnkX
        } else {
            if(rnkX < ncol(X)){
                msg <- "A protein intensity model is not identifiable! "
                msg <- paste0(msg, "RowData first line:\n")
                cat(msg)
                print(rowData(data)[startIndex,])
                stop()
            } else {
              keepAlpha <- rep(TRUE, ncol(X))
                MAlpha <- ncol(X) 
            }
        }
        
        #browser()
        TT <- as.matrix(model.matrix(object = formM, data = df))
        
        keepGamma <- rep(FALSE, ncol(TT))
        # apply tolerance slightly lower than default (1e-7) to keep the 
        # slightly-non-multicollinear column
        qr.TT <- qr(TT, tol=1e-9, LAPACK = FALSE)
        rnkTT <- qr.TT$rank
        
        # if (cleanNonIdentifiable){
        #     keep[qr.TT$pivot[seq_len(rnkTT)]] <- TRUE
        #     TT <- TT[, keep, drop = FALSE]
        #     MGamma <- rnkTT
        # } else {
            if(rnkTT < ncol(TT)){
                msg <- "A protein missingness model is not identifiable! "
                msg <- paste0(msg, "RowData first line:\n")
                cat(msg)
                print(rowData(data)[startIndex,])
                stop()
            } else {
                keepGamma <- rep(TRUE, ncol(TT))
                MGamma <- ncol(TT) 
            }
        #}
        
        #browser()
        rdmEffectsStruct <- lme4::findbars(term = form)
        nRdm <- length(rdmEffectsStruct)
        
        Z <- matrix(nrow = 0, ncol = 0)
        MBetaPerRdmEffect <- numeric(0)
        paramNamesSigma2Rdm <- character(0)
        
        #specific contrast options for random effects
        options(contrasts=c('contr.sum','contr.poly'))
        
        for(r in seq_along(rdmEffectsStruct)){
            if(deparse(rdmEffectsStruct[[r]][1]) != "`|`()" || 
               deparse(rdmEffectsStruct[[r]][2]) != "1()"){
                msg <- "Error in model formula : currently only random effects "
                msg <- paste0(msg, "of the form (1|variable) are supported.")
                stop(msg)
            }
            currentTerm <- stringr::str_replace_all(
                deparse(rdmEffectsStruct[[r]][3]), "[()]", "")
                
            currentForm <- reformulate(termlabels = currentTerm, intercept = FALSE)
            currentZ <- model.matrix(currentForm, df)
            if(nrow(Z) == 0){
                Z <- currentZ
            } else {
                Z <- cbind(Z, currentZ)
            }
            
            MBetaPerRdmEffect <- c(MBetaPerRdmEffect, ncol(currentZ))
            paramNamesSigma2Rdm <- c(paramNamesSigma2Rdm, currentTerm)
        }
        
        MBeta <- ncol(Z)
        
        paramNamesAlpha <- stringr::str_replace_all(colnames(X), "[()]", "")
        paramNamesBeta <- stringr::str_replace_all(colnames(Z), "[()]", "")
        paramNamesGamma <- stringr::str_replace_all(colnames(TT), "[()]", "")
        
        #browser()
        obsIndex <- !is.na(df$y)
        misIndex <- is.na(df$y)
        obsRuns <- df$run[obsIndex]
        misRuns <- df$run[misIndex]
        obsPepts <- df$peptide[obsIndex]
        misPepts <- df$peptide[misIndex]
        paramNamesYObs <- paste0(obsRuns, "_", obsPepts)
        paramNamesYMis <- paste0(misRuns, "_", misPepts)
        
        NObs <- sum(obsIndex)
        NMis <- sum(misIndex)
        
        yObs <- df$y[obsIndex]
        
        XObs <- X[obsIndex,]
        XMis <- X[misIndex,]
        XObs <- as.vector(XObs) # note as.vector(matrix) uses column major order
        XMis <- as.vector(XMis)
        
        ZObs <- numeric(0)
        ZMis <- numeric(0)
        if( nRdm > 0){
            ZObs <- Z[obsIndex,]
            ZMis <- Z[misIndex,]
            ZObs <- as.vector(ZObs)
            ZMis <- as.vector(ZMis)
        }
        
        
        TObs <- TT[obsIndex,]
        TMis <- TT[misIndex,]
        TObs <- as.vector(TObs)
        TMis <- as.vector(TMis)
        
        # handle contrasts
        
        options(contrasts=c('contr.treatment','contr.poly'))
        
        nTentativeContrast <- length(contrasts)
        nActiveConstrast = 0
        contrastNames = c()
        contrastCoefs = numeric(0)
        if(nTentativeContrast > 0){
          
          for(c in 1:nTentativeContrast){
            contrastName <- contrasts[[c]]$contrastName
            if(is.null(contrastName)){
              msg <- paste0("specified contrast without name detected, skipped!")
              warning(msg)
            }
            factorName <- contrasts[[c]]$factor
            factorCol <- as.factor(df[,factorName])
            if(is.null(factorCol)){
              msg <- paste0("issue with definition of constrast [", contrastName,
                            "] : variable [", factorName, 
                            "] does not exist or can't be converted to a factor!")
              warning(msg)
              next
            }
            factorLevels <- levels(factorCol)
            df.contrast <- rbind(df[1,,drop=FALSE], df[1,,drop=FALSE])
            levels2Compare <- contrasts[[c]]$levels
            if(length(levels2Compare) != 2){
              msg <- paste0("issue with definition of constrast [", contrastName,
                            "] : levels should be on length 2!")
              warning(msg)
              next
            }
            if(sum(levels2Compare %in% factorLevels) != 2){
              msg <- paste0("issue with definition of constrast [", contrastName,
                            "] : levels to compare not found in dataframe!")
              warning(msg)
              next
            }
            
            df.contrast[, factorName] <- c(contrasts[[c]]$levels[1], 
                                           contrasts[[c]]$levels[2])
            
            # now convert all cols to factors if possible, using the initial
            # levels
            varsInFormula <- all.vars(formFixed)
            
            for(v in varsInFormula){
              varCol <- as.factor(df[,v])
              if(is.null(varCol)){
                msg <- paste0("Definition of constrast [", contrastName,
                              "] : variable [", v, 
                              "] does not exist or can't be converted to a factor!")
                warning(msg)
                next
              }
              varFactorLevels <- levels(varCol)
              df.contrast[, v] <- factor(df.contrast[, v],
                                         levels = varFactorLevels)
            }
            
            # df.contrast[,factorName] <- factor(c(contrasts[[c]]$levels[1], 
            #                                      contrasts[[c]]$levels[2]), 
            #                                    levels = factorLevels)
            #ff <- reformulate(termlabels = factorName)
            modelMatrix.contrast <- model.matrix(
              formFixed,
              df.contrast)
            contrast.coef <- modelMatrix.contrast[1,] - modelMatrix.contrast[2,]
            
            nonNullCoefs <- which(contrast.coef != 0)
            if(sum(keepAlpha[nonNullCoefs]) != length(nonNullCoefs)){
              msg <- paste0("contrast [", contrastName,
                            "] was referring to some removed alpha coefficients",
                            " (collinearity) => contrast skipped!")
              warning(msg)
              next
            }
            
            
            # update outputs
            contrastNames <- c(contrastNames, contrastName)
            nActiveConstrast <- nActiveConstrast + 1
            contrastCoefs <- c(contrastCoefs, contrast.coef[keepAlpha])
          }
        }
        

        res <- list(NObs = NObs,
                    NMis = NMis,
                    MAlpha = MAlpha,
                    MBeta = MBeta,
                    MBetaPerRdmEffect = MBetaPerRdmEffect,
                    MGamma = MGamma,
                    nRdm = nRdm,
                    yObs = yObs,
                    XObs = XObs,
                    XMis = XMis,
                    ZObs = ZObs,
                    ZMis = ZMis,
                    TObs = TObs,
                    TMis = TMis,
                    paramNamesAlpha = paramNamesAlpha,
                    paramNamesBeta = paramNamesBeta,
                    paramNamesGamma = paramNamesGamma,
                    paramNamesSigma2Rdm = paramNamesSigma2Rdm,
                    paramNamesYObs = paramNamesYObs,
                    paramNamesYMis = paramNamesYMis,
                    obsRuns = obsRuns,
                    misRuns = misRuns,
                    obsPepts = obsPepts,
                    misPepts = misPepts,
                    contrastNames = contrastNames,
                    nContrast = nActiveConstrast,
                    contrastCoefs = contrastCoefs)
        res
    } # end function definition

    #browser()
    if(nProt == 1){
        res <- buildOneModelInputs(data = data,
                                   form = form,
                                   formM = formM,
                                   startIndex = 1,
                                   nRows = nrow(data),
                                   cleanNonIdentifiable = cleanNonIdentifiable,
                                   contrasts = contrasts)
        if(res$nRdm == 1){
            MBetaPerRdmEffect = array(res$MBetaPerRdmEffect, dim=1)
        } else {
            MBetaPerRdmEffect = res$MBetaPerRdmEffect
        }
        inputs <- list(nProt = nProt,
                    nRun = nRun,
                    nRdm = res$nRdm,
                    nPept = array(proteinInfo$nPept, dim=1),
                    MAlpha = array(res$MAlpha, dim=1),
                    MBeta = array(res$MBeta, dim=1),
                    MBetaPerRdmEffect = MBetaPerRdmEffect,
                    MGamma = array(res$MGamma, dim=1),
                    NObs = array(res$NObs, dim=1),
                    NMis = array(res$NMis, dim=1),
                    yObs = res$yObs,
                    XObs = res$XObs,
                    XMis = res$XMis,
                    ZObs = res$ZObs,
                    ZMis = res$ZMis,
                    TObs = res$TObs,
                    TMis = res$TMis,
                    nContrast = array(res$nContrast, dim=1),
                    contrastCoefs = res$contrastCoefs)
        #browser()
        paramNames <- list(paramNamesAlpha = list(res$paramNamesAlpha),
                           paramNamesBeta = list(res$paramNamesBeta),
                           paramNamesGamma = list(res$paramNamesGamma),
                           paramNamesSigma2Rdm = list(res$paramNamesSigma2Rdm),
                           paramNamesYObs = list(res$paramNamesYObs),
                           paramNamesYMis = list(res$paramNamesYMis),
                           obsRuns = list(res$obsRuns),
                           misRuns = list(res$misRuns),
                           obsPepts = list(res$obsPepts),
                           misPepts = list(res$misPepts),
                           contrastNames = list(res$contrastNames))
                           
    } else {
        res <- mapply( buildOneModelInputs,
                       startIndex = proteinInfo$startPeptideIndex,
                       nRows = proteinInfo$nPept,
                       MoreArgs = list( data = data,
                                        form = as.formula(form),
                                        formM = as.formula(formM),
                                        cleanNonIdentifiable = cleanNonIdentifiable,
                                        contrasts = contrasts),
                       SIMPLIFY = TRUE)

        # sum_if_not_zero <- function(x,y){
        #     ifelse(x == 0, 0, sum(x,y,na.rm = TRUE))
        # }
        
        
        inputs <- list(nProt = nProt,
                    nRun = nRun,
                    nRdm = unlist(res["nRdm",]),
                    nPept = proteinInfo$nPept,
                    MAlpha = unlist(res["MAlpha",]),
                    MBeta = unlist(res["MBeta",]),
                    MBetaPerRdmEffect = unlist(res["MBetaPerRdmEffect",]),
                    MGamma = unlist(res["MGamma",]),
                    NObs = unlist(res["NObs",]),
                    NMis = unlist(res["NMis",]),
                    yObs = unlist(res["yObs",]),
                    XObs = unlist(res["XObs",]),
                    XMis = unlist(res["XMis",]),
                    TObs = unlist(res["TObs",]),
                    TMis = unlist(res["TMis",]),
                    ZObs = unlist(res["ZObs",]),
                    ZMis = unlist(res["ZMis",]),
                    nContrast = unlist(res["nContrast",]),
                    contrastCoefs = unlist(res["contrastCoefs",]))
        paramNames <- list(paramNamesAlpha = res["paramNamesAlpha",],
                           paramNamesBeta = res["paramNamesBeta",],
                           paramNamesGamma = res["paramNamesGamma",],
                           paramNamesSigma2Rdm = res["paramNamesSigma2Rdm",],
                           paramNamesYObs = res["paramNamesYObs",],
                           paramNamesYMis = res["paramNamesYMis",],
                           obsRuns = res["obsRuns",],
                           misRuns = res["misRuns",],
                           obsPepts = res["obsPepts",],
                           misPepts = res["misPepts",],
                           contrastNames = res["contrastNames",]
        )
                             
    }
    
    res <- list(inputs = inputs, 
                paramNames = paramNames, 
                proteinInfo = proteinInfo)
}


##' This function takes a data set and a formula and initiates a
##' `ProtModel`.
##'
##' @title Create a ProtModel object
##' 
##' @param modelName A `character(1)` containing the name of the model.
##' 
##' @param data A `SummarizedExperiment` object.
##' 
##' @param formula A `formula` used to define the main intensity linear model
##' 
##' @param formulaM A `formula` used to define the missingness model
##' 
##' @param protColName A character(1) which designates the name of the column
##' in `data@rowData` that contains the protein identifiers
##' 
##' @param peptColName A character(1) which designates the name of the column
##' in `data@rowData` that contains the peptide identifiers
##' 
##' @param cleanNonIdentifiable A boolean indicating whether the engine should
##' try to remove regressors for those protein models that are not identifiable
##' 
##' @return An initialized `StanModel` object.
##'
##' @import dplyr
##' 
##' @importFrom magrittr %>%
##' 
##' @importFrom stats qnorm terms
##' 
##' @export
##' 
##' @author Laurent Gatto
makeModel <- function(modelName, data, formula, formulaM, protColName, peptColName,
                      cleanNonIdentifiable = FALSE, refIntensity = 0.,
                      contrasts = list()) {
    
    # sanity checks of user input
    stopifnot(inherits(data, "SummarizedExperiment"))
    if(!(protColName %in% colnames(rowData(data))))
        stop(paste0(protColName, " not found in rowData(data)"))
    if(!(peptColName %in% colnames(rowData(data))))
        stop(paste0(peptColName, " not found in rowData(data)"))
    
    form = as.formula(formula)
    formM = as.formula(formulaM)
    # nor formula can contain term 'y'
    if("y" %in% colnames(attr(terms(form),"factors")))
    {
        stop("intensity model formula can't contain [y] term");
    }
    if("y" %in% colnames(attr(terms(formM),"factors")))
    {
        msg <- "missing model formula can't contain explicit [y] term, "
        msg <- paste0(msg, "use [NMAR] argument in runModel() instead")
        stop(msg)
    }
    
    # store contrasts options
    old.contrasts <- options('contrasts')$contrasts
    # options(contrasts=c('contr.sum','contr.poly'))
    # 
    # order data by protein, peptide
    .dataOrdered <- data[order(rowData(data)[,protColName],
                               rowData(data)[,peptColName]),]


    #cleanNonIdentifiable <- FALSE
    
    res <- buildModelInputs( .dataOrdered, 
                             form = form,
                             formM = formM,
                             protColName = protColName,
                             peptColName = peptColName,
                             cleanNonIdentifiable = cleanNonIdentifiable,
                             contrasts = contrasts)
    
    InputData <- res$inputs
    if(InputData$nProt == 1){
      InputData$refIntensity <- array(refIntensity, dim = 1) 
    } else {
      InputData$refIntensity <- rep(refIntensity, InputData$nProt)
    }
    
    InputData$NObsSum <- sum(InputData$NObs)
    InputData$NMisSum <- sum(InputData$NMis)
    InputData$XObsLen <- length(InputData$XObs)
    InputData$XMisLen <- length(InputData$XMis)
    InputData$ZObsLen <- length(InputData$ZObs)
    InputData$ZMisLen <- length(InputData$ZMis)
    InputData$TObsLen <- length(InputData$TObs)
    InputData$TMisLen <- length(InputData$TMis)
    
    # the following vector is in fact a scalar, not dependent on protein
    InputData$nRdm <- InputData$nRdm[1]
    
    InputData$contrastCoefLen <- length(InputData$contrastCoefs)
    
    theNames <- names(InputData)
    theNames[which(theNames=="nProt")] <- "K"
    theNames[which(theNames=="nRun")] <- "R"
    theNames[which(theNames=="nPept")] <- "P"
    names(InputData) <- theNames
    
    #restore contrasts options
    options(contrasts=old.contrasts)
    
    x <- ProtModel(.modelName = modelName,
                   .data = .dataOrdered,
                   .formula = formula,
                   .formulaM = formulaM,
                   .paramNames = res$paramNames,
                   .InputData = InputData,
                   .proteinInfo = res$proteinInfo)
    
    #x <- updateCoefficients(x)
    
    x
              
}
