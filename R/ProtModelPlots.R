##' @include ProtModel.R
NULL

##' Plot either one or two distributions for a given parameter.
##' 
##' @title Plots and compare the distributions for a given parameter
##' 
##' @param x1 A `ProtModel` object which should have run already
##' 
##' @param x2 An optional `ProtModel` object which should have run already 
##' 
##' @param paramName A `character` containing the name of the parameter to plot
##'
##' @return a ggplot graph showing both distributions for the given parameter
##' 
##' @import ggplot2
##' @import patchwork
##' 
##' @importFrom KernSmooth dpih
##'
##' @export
##' 
##' @author Philippe Hauchamps
plotParamDistr <- function(x1, x2 = NULL, x3 = NULL, paramName, 
                           xlim = NULL, ylim = NULL, keepBinWidth = TRUE,
                           vline_xintercept = NULL)
{
  stopifnot(inherits(x1, "ProtModel"))
  if(x1@run == 0) 
    stop("Parameter distribution comparison needs input ProtModel objects to have run")
  if(!is.null(x2)){
    stopifnot(inherits(x2, "ProtModel"))
    if(x2@run == 0) 
      stop("Parameter distribution comparison needs input ProtModel objects to have run")
    if(!is.null(x3)){
      stopifnot(inherits(x3, "ProtModel"))
      if(x3@run == 0) 
        stop("Parameter distribution comparison needs input ProtModel objects to have run")
      nHistograms <- 3
    } else {
      nHistograms <- 2
    }
  } else {
    nHistograms <- 1
  }
  
  samples <- list()
  protModelNames <- rep("no_name",  nHistograms)
  samples[[1]] <- getParamSample(x1, paramName)
  protModelNames[1] <- modelName(x1)
  if(!is.null(x2)){
    samples[[2]] <- getParamSample(x2, paramName)
    protModelNames[2] <- modelName(x2)
  }
  if(!is.null(x3)){
    samples[[3]] <- getParamSample(x3, paramName)
    protModelNames[3] <- modelName(x3)
  }
  
  if(nHistograms > 1){
    g <- list()
    # binwidth ruled (arbitrarily) by histogram 1
    refBinWidth <- dpih(samples[[1]])
    for(i in 1:nHistograms){
      STAN_df <- data.frame(paramValue = samples[[i]])
      if(keepBinWidth == FALSE && i>1){
        binWidth <- dpih(samples[[i]])
      } else {
        binWidth <- refBinWidth
      }
      g[[i]] <- ggplot(data=STAN_df, aes(x=paramValue)) +
        geom_histogram(mapping = aes(y=after_stat(density)),
                       binwidth = binWidth,
                       alpha = 0.5,
                       colour = "white",
                       fill = "blue") +
        labs(x=paramName, y="Density",
             title = paste0("Distribution of ", paramName),
             subtitle = paste0(" for object ", protModelNames[i]))
      if(!is.null(xlim) || !is.null(ylim)){
        g[[i]] <- g[[i]] +
          coord_cartesian(xlim = xlim, ylim = ylim)
      }
      if(!is.null(vline_xintercept)){
        g[[i]] <- g[[i]] +
          geom_vline(xintercept = vline_xintercept, 
                     linetype = "dashed", 
                     color = "red")
      }
    }
    if(nHistograms == 3){
      g[[1]] + g[[2]] + g[[3]]
    } else {
      g[[1]] + g[[2]]
    }
    
  } else {
    STAN_df <- data.frame(paramValue = samples[[1]])
    binWidth <- dpih(samples[[1]])
    g1 <- ggplot(data=STAN_df, aes(x=paramValue)) +
      geom_histogram(mapping = aes(y=after_stat(density)),
                     binwidth = binWidth,
                     alpha = 0.5,
                     colour = "white",
                     fill = "blue") +
      labs(x=paramName, y="Density",
           title = paste0("Distribution of ", paramName))
    if(!is.null(xlim) || !is.null(ylim)){
      g1 <- g1 +
        coord_cartesian(xlim = xlim, ylim = ylim)
    }
    if(!is.null(vline_xintercept)){
      g1 <- g1 +
        geom_vline(xintercept = vline_xintercept, 
                   linetype = "dashed", 
                   color = "red")
    }
    g1
  }
  
}

##' @import ggplot2
##' 
##' @export
##' 
##' @author Philippe Hauchamps
ProtModelVolcanoPlot <- function(x,
                          variableName,
                          lFCThreshold,
                          testType,
                          estimateType = "median",
                          targetsDA = NULL,
                          xlim = NULL){
  DF <- computeVolcanoPlotData(x,
                               variableName = variableName,
                               lFCThreshold = lFCThreshold,
                               testType = testType,
                               estimateType = estimateType
  )
  
  if(is.null(targetsDA)){
    graph <- ggplot(data = DF,
                    aes(x = lFC, y = DAProb)) +
      geom_point(alpha = 1.0) +
      labs(x = "logFoldChange",
           y = "prob(DA)",
           title = paste0("Volcano plot of ",modelName(x)),
           subtitle = paste0("lFC threshold = ", lFCThreshold,
                             "; testType = ", testType))
  } else {
    DF2 <- data.frame(lFC = DF$lFC,
                      DAProb = DF$DAProb,
                      bench = ifelse(targetsDA,"DA","not DA"),
                      alpha = ifelse(targetsDA, 1., 0.5))
    
    DF2 <- DF2[!is.na(targetsDA),]
    
    graph <- ggplot(data = DF2,
                    aes(x = lFC, y = DAProb, 
                        colour = bench, 
                        alpha = alpha)) +
      geom_point() +
      labs(x = "logFoldChange",
           y = "prob(DA)",
           title = paste0("Volcano plot of ",modelName(x)),
           subtitle = paste0("lFC threshold = ", lFCThreshold,
                             "; testType = ", testType)) + 
      guides(alpha=FALSE) 
  }
  if(is.null(xlim)){
    graph <- graph + coord_cartesian(ylim = c(0,1))
  } else {
    graph <- graph + coord_cartesian(xlim = xlim,
                                     ylim = c(0,1))
  }
  
  graph
}

##' @import ggplot2
##' 
##' @export
##' 
##' @author Philippe Hauchamps
ProtModelROCPlot <- function(x,
                             variableName,
                             lFCThreshold,
                             testType,
                             targetsDA,
                             xabs = c("FPR","FDP"),
                             xlim = NULL,
                             ylim = NULL){
  
  if(is.list(x)){
    Models <- x
  } else {
    Models <- list(x)
  }
  
  chosenX <- xabs[1]
  
  if(!(chosenX == "FPR" | chosenX == "FDP")){
    msg <- paste0("Unrecognized xasb parameter passed : [", xabs[1], "].")
    msg <- paste0(msg,"Possible values are : [FPR] and [FDP]")
    stop(msg)
  }
  
  chosenY <- "sensi"
  
  
  
  resROC <- lapply(Models,
                   FUN = computeROCCurve,
                   variableName = variableName,
                   lFCThreshold = lFCThreshold,
                   testType = testType,
                   perProtTargetDAs = targetsDA)
  
  for(i in seq_along(Models)){
    currentDF <- resROC[[i]]$ROC
    # force the right handling of constant abscissa values
    currentDF[, chosenX] <- currentDF[, chosenX] + seq(nrow(currentDF),1,-1)*1e-12
    currentDF <- cbind(currentDF, model = modelName(Models[[i]]))
    if(i>1){
      DF <- rbind(DF, currentDF)
    } else {
      DF <- currentDF
    }
  }
                            
  graph <- ggplot(data = DF, 
                  aes_string(x = chosenX, y = chosenY)) +
    aes(colour = model) + 
    geom_line(size = 1.5) + 
    coord_cartesian(xlim = xlim, ylim = ylim) + 
    labs(title = "ROC Curves",
         subtitle = paste0("lFC threshold = ", lFCThreshold,
                           "; testType = ", testType))
  
  
  graph
    
}

##' @import ggplot2
##' 
##' @importFrom KernSmooth dpih
##' 
##' @export
##' 
##' @author Philippe Hauchamps
ProtModelVarHistogram <- function(x, 
                              variableName, 
                              statistic = "mean",
                              xlim = NULL,
                              ylim = NULL){
  mySum <- getProtModelSummary(x, 
                               pattern = variableName)
  
  if(nrow(mySum) == 0){
    msg <- paste0("Variable name [", variableName, "] not found")
    stop(msg)
  }
  
  if(statistic == "median"){
    chosenStatistic <- "50%"
  } else {
    chosenStatistic <- statistic
  }
  
  if(!chosenStatistic %in% colnames(mySum)){
    msg <- paste0("variable [", statistic, "] not in available columns")
    stop(msg)
  }
  
  bw <- KernSmooth::dpih(mySum[,chosenStatistic])
  
  if(str_detect(chosenStatistic, "%")){
    chosenStatistic <- paste0("`", chosenStatistic, "`")
  }
  
  graph <- ggplot(data = mySum, aes_string(x=chosenStatistic)) +
    geom_histogram(mapping = aes(y = after_stat(density)),
                   binwidth = bw,
                   alpha = 0.5,
                   colour = "white",
                   fill = "blue") +
    labs(title = paste0(variableName, " per protein for ProtModel ", 
                        modelName(x)),
         subtitle = paste0("statistic : ", statistic),
         x = variableName,
         y = "density") +
    coord_cartesian(xlim = xlim, ylim = ylim)
  
  graph
  
}
                 
  