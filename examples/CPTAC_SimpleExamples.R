require(StanProt)
require(SummarizedExperiment)
require(tools) # needed for file_path_sans_ext


summary(CPTAC1347prot)
medianIntensity <- median(assay(CPTAC1347prot), na.rm = TRUE)
medianIntensity
mySeed <- 24052021


############################################################################
# Example 1 : run on one single protein, using Single Protein Model        #
#             compare STAN NUTS with MFVB and Laplace Approximation        #
############################################################################

selectedProts <- c("P08311ups")

data.selected <- 
  CPTAC1347prot[rowData(CPTAC1347prot)$protein %in% selectedProts ,]


# STEP 1 : make model
# intensity part
form <- ~ condition + lab + (1|peptide) 
# missingness part
formM <- ~ 1 + lab

# define contrasts of interest
myContrasts <- list()
myContrasts[[1]] <- list(contrastName = "contrast_6B-6A",
                         factor = "condition",
                         levels = c("6B", "6A"))

myContrasts[[2]] <- list(contrastName = "contrast_6C-6B",
                         factor = "condition",
                         levels = c("6C", "6B"))

myProtModel <- makeModel( modelName = "Stan NUTS",
                          data = data.selected,
                          formula = form,
                          formulaM = formM,
                          protColName = "protein",
                          peptColName = "sequence",
                          refIntensity = medianIntensity,
                          contrasts = myContrasts)

show(myProtModel)
proteinInfo(myProtModel)

# STEP 2 : run model

# prior distributions of alpha coefficients(in log2) => N(mu,sd)
intensityLow <- 0
intensityHigh <- 30
foldChangeLow <- -10
foldChangeHigh <- 10

# prior distribution of variances => ScaledInvChi2(nu,scale)
NuEpsilon = 0.002
ScaleEpsilon = 1
NuRdm = 0.002
ScaleRdm = 1

gammaInterceptMean = 0
gammaInterceptSd = 10
gammaTermsMean = 0
gammaTermsSd = 10
gammaYMean = 0
gammaYSd = 10

# model ingredients
# student-t likelihood(1) or normal(0)
t_likelihood <- TRUE
# fixed degrees of freedom for student-t
t_fixedDF <- 4
# model missingness or not
modelMiss <- TRUE
# NMAR component or not
NMAR <- TRUE
# hierarch. set-up for sigma
sigmaHierarch <- FALSE
# hierarch. set-up for sigma of random effects(if any)
sigmaRdmHierarch <- FALSE
# hierarch. set-up for gamma (if modelMiss)
gammaHierarch <- FALSE

modelSpecPars = ProtModelSpecPars( 
  .priorIntensityLow = intensityLow,
  .priorIntensityHigh = intensityHigh,
  .priorFoldChangeLow = foldChangeLow,
  .priorFoldChangeHigh = foldChangeHigh,
  .priorGammaInterceptMean = gammaInterceptMean,
  .priorGammaTermsMean = gammaTermsMean,
  .priorGammaInterceptSd = gammaInterceptSd,
  .priorGammaTermsSd = gammaTermsSd,
  .priorGammaYMean = gammaYMean,
  .priorGammaYSd = gammaYSd,
  .priorNuEpsilon = NuEpsilon,
  .priorScaleEpsilon = ScaleEpsilon,
  .priorNuRdm = NuRdm,
  .priorScaleRdm = ScaleRdm,
  .t_likelihood = t_likelihood,           
  .t_fixedDF = t_fixedDF,                 
  .modelMiss = modelMiss,                 
  .NMAR = NMAR,                           
  .sigmaHierarch = sigmaHierarch,         
  .sigmaRdmHierarch = sigmaRdmHierarch,   
  .gammaHierarch = gammaHierarch          
)

# STAN NUTS arguments
stanSamplingArgs <- list(seed = mySeed,
                         cores = 1,
                         iter = 2000,
                         chains = 4,
                         refresh = 100)#,
                         # control = list(
                         #   max_treedepth = 10,
                         #   adapt_delta = 0.85)
                         #byProt = FALSE

# Stan model script file => provided with StanProt package
StanModelScript <- file_path_sans_ext(
  system.file("Stan", 
              "FullHierarchical.stan", 
              package = "StanProt"))
                             
# MODEL 1 : exact inference using STAN NUTS 

myProtModel <- runModel(myProtModel,
                        providerPackage = "stan",
                        modelScript = StanModelScript,
                        runMode = "sampling",
                        modelSpecPars = modelSpecPars,
                        algoArgs = stanSamplingArgs)

getElapsedTime(ModelFit(myProtModel))
checkAllDiagnostics(ModelFit(myProtModel))
show(myProtModel)


myProtModel <- 
  computeProtModelSummary(myProtModel)
                                       
(sm1 <- round(getProtModelSummary(myProtModel,
                                  pattern = c("alpha", "gamm", "beta","nu0",
                                              "s0", "mu", "sigma2", "contr"),
                                  whichProtein = 1), 2))


# MODEL 2 : approximate inference using MFVB
MFVBAlgoArgs <- list(nIterMax = 50,
                     relTol = 1e-6,
                     ELBO_each = 1,
                     refresh = 10,
                     cores = 1)

myProtModel2 <- myProtModel
modelName(myProtModel2) <- "MFVB"

myProtModel2 <- runModel(myProtModel2,
                        providerPackage = "native",
                        runMode = "MFVB",
                        modelSpecPars = modelSpecPars,
                        algoArgs = MFVBAlgoArgs)

showFitResults(ModelFit(myProtModel2), whichProtein = 1)

myProtModel2 <- computeProtModelSummary(myProtModel2)
(sm2 <- round(getProtModelSummary(myProtModel2), 3))


# MODEL 3 : Approximation of posterior distribution with Laplace approximation, 
# based on model parameters split (psi, phi)

# Stan model script file => provided with StanProt package
StanModelScript2 <- file_path_sans_ext(
  system.file("Stan", 
              "FullHierarchicalLaplaceApproxi2.stan", 
              package = "StanProt"))

stanMAPArgs <- list(seed = mySeed,
                    draws = 4000,
                    hessian = FALSE,
                    algorithm = "LBFGS",
                    verbose = TRUE)
                     
myProtModel3 <- myProtModel
modelName(myProtModel3) <- "Laplace Approxi (cond/marg)"

start_time <- Sys.time()
myProtModel3 <- runModel(myProtModel3,
                         providerPackage = "stan",
                         modelScript = StanModelScript2,
                         runMode = "MAP",
                         modelSpecPars = modelSpecPars,
                         algoArgs = stanMAPArgs)
end_time <- Sys.time()
elapsed = end_time - start_time
cat("Elapsed time : ", round(elapsed,1), attributes(elapsed)$units)

myProtModel3 <- computeProtModelSummary(myProtModel3)
(sm3 <- round(getProtModelSummary(myProtModel3, 
                      pattern = c("alpha", "gamm", "beta","nu0","s0", "mu", 
                                  "sigma2", "contr"),
                      whichProtein = 1),
              3))


# compare posterior distributions

mod1 <- myProtModel
mod2 <- myProtModel2
mod3 <- myProtModel3


plotParamDistr(mod1, mod2, mod3,
               paramName = "contrast_6B-6A[1]",
               xlim = c(1.2, 2.0), ylim = c(0.0, 7.0))

plotParamDistr(mod1, mod2, mod3,
               paramName = "contrast_6C-6B[1]",
               xlim = c(1.2, 2.1), ylim = c(0.0, 7.0))

plotParamDistr(mod1, mod2, mod3,
               paramName = "alpha_Intercept[1]",
              xlim = c(13, 19), ylim = c(0.0, 1.0))

plotParamDistr(mod1, mod2, mod3,
               paramName = "alpha_lab2[1]",
               xlim = c(-0.4,0.8), ylim = c(0,4))

plotParamDistr(mod1, mod2, mod3,
               paramName = "sigma2Eps[1]",
               keepBinWidth = TRUE,
               xlim = c(0.0,0.25), ylim = c(0,80))

############################################################################
# Example 2 : run on two proteins, using Single Protein Model              #
#             compare STAN NUTS : full model vs. simple linear model       #
############################################################################


selectedProts <- c("P00441ups", "P11938")

data.selected <- 
  CPTAC1347prot[rowData(CPTAC1347prot)$protein %in% selectedProts ,]


# STEP 1 : make model
form <- ~ condition + lab + (1|peptide) 
formM <- ~ 1 + lab

myContrasts <- list()
myContrasts[[1]] <- list(contrastName = "contrast_6B-6A",
                         factor = "condition",
                         levels = c("6B", "6A"))

myContrasts[[2]] <- list(contrastName = "contrast_6C-6B",
                         factor = "condition",
                         levels = c("6C", "6B"))

myProtModel <- makeModel( modelName = "Stan NUTS - Full model",
                          data = data.selected,
                          formula = form,
                          formulaM = formM,
                          protColName = "protein",
                          peptColName = "sequence",
                          refIntensity = medianIntensity,
                          contrasts = myContrasts)

show(myProtModel)
proteinInfo(myProtModel)

# STEP 2 : run model

# prior distributions of alpha coefficients(in log2) => N(mu,sd)
intensityLow <- 0
intensityHigh <- 30
foldChangeLow <- -10
foldChangeHigh <- 10

# prior distribution of variances => ScaledInvChi2(nu,scale)
NuEpsilon = 0.002
ScaleEpsilon = 1
NuRdm = 0.002
ScaleRdm = 1

gammaInterceptMean = 0
gammaInterceptSd = 10
gammaTermsMean = 0
gammaTermsSd = 10
gammaYMean = 0
gammaYSd = 10

# model ingredients
# student-t likelihood(1) or normal(0)
t_likelihood <- TRUE
# fixed degrees of freedom for student-t
t_fixedDF <- 4
# model missingness or not
modelMiss <- TRUE
# NMAR component or not
NMAR <- TRUE
# hierarch. set-up for sigma
sigmaHierarch <- FALSE
# hierarch. set-up for sigma of random effects(if any)
sigmaRdmHierarch <- FALSE
# hierarch. set-up for gamma (if modelMiss)
gammaHierarch <- FALSE

modelSpecPars = ProtModelSpecPars( 
  .priorIntensityLow = intensityLow,
  .priorIntensityHigh = intensityHigh,
  .priorFoldChangeLow = foldChangeLow,
  .priorFoldChangeHigh = foldChangeHigh,
  .priorGammaInterceptMean = gammaInterceptMean,
  .priorGammaTermsMean = gammaTermsMean,
  .priorGammaInterceptSd = gammaInterceptSd,
  .priorGammaTermsSd = gammaTermsSd,
  .priorGammaYMean = gammaYMean,
  .priorGammaYSd = gammaYSd,
  .priorNuEpsilon = NuEpsilon,
  .priorScaleEpsilon = ScaleEpsilon,
  .priorNuRdm = NuRdm,
  .priorScaleRdm = ScaleRdm,
  .t_likelihood = t_likelihood,           
  .t_fixedDF = t_fixedDF,                 
  .modelMiss = modelMiss,                 
  .NMAR = NMAR,                           
  .sigmaHierarch = sigmaHierarch,         
  .sigmaRdmHierarch = sigmaRdmHierarch,   
  .gammaHierarch = gammaHierarch          
)

# STAN NUTS arguments
stanSamplingArgs <- list(seed = mySeed,
                         cores = 1,
                         iter = 2000,
                         chains = 4,
                         refresh = 100)#,
# control = list(
#   max_treedepth = 10,
#   adapt_delta = 0.85)
#byProt = FALSE

# Stan model script file => provided with StanProt package
StanModelScript <- file_path_sans_ext(
  system.file("Stan", 
              "FullHierarchical.stan", 
              package = "StanProt"))

# MODEL 1 : STAN NUTS with full model

myProtModel <- runModel(myProtModel,
                        providerPackage = "stan",
                        modelScript = StanModelScript,
                        runMode = "sampling",
                        modelSpecPars = modelSpecPars,
                        algoArgs = stanSamplingArgs)

getElapsedTime(ModelFit(myProtModel))
checkAllDiagnostics(ModelFit(myProtModel))

myProtModel <- 
  computeProtModelSummary(myProtModel)

(sm1 <- round(getProtModelSummary(myProtModel,
                                  pattern = c("alpha", "gamm", "beta","nu0",
                                              "s0", "mu", "sigma2", "contr"),
                                  whichProtein = c(1,2)), 2))

# MODEL 2 : STAN NUTS with simple linear model
modelSpecPars2 <- modelSpecPars
t_likelihood(modelSpecPars2) <- FALSE
modelMiss(modelSpecPars2) <- FALSE
myProtModel2 <- myProtModel
modelName(myProtModel2) <- "Stan NUTS - simple linear model"
myProtModel2 <- runModel(myProtModel2,
                         providerPackage = "stan",
                         modelScript = StanModelScript,
                         runMode = "sampling",
                         modelSpecPars = modelSpecPars2,
                         algoArgs = stanSamplingArgs)
                        
getElapsedTime(ModelFit(myProtModel2))
checkAllDiagnostics(ModelFit(myProtModel2))

myProtModel2 <- 
  computeProtModelSummary(myProtModel2)

(sm2 <- round(getProtModelSummary(myProtModel2,
                                  pattern = c("alpha", "gamm", "beta","nu0",
                                              "s0", "mu", "sigma2", "contr"),
                                  whichProtein = c(1,2)), 2))

# comparing posterior distributions of contrasts
mod1 <- myProtModel
mod2 <- myProtModel2
mod3 <- NULL


plotParamDistr(mod1, mod2, mod3,
               paramName = "contrast_6B-6A[1]",
               xlim = c(-2, 4), ylim = c(0.0, 0.8))

plotParamDistr(mod1, mod2, mod3,
               paramName = "contrast_6C-6B[1]",
               xlim = c(-1, 4), ylim = c(0.0, 1.0))

# USE ANOTHER PROTEIN for PROTEIN 2 (P02294)

plotParamDistr(mod1, mod2, mod3,
               paramName = "contrast_6B-6A[2]",
               xlim = c(-1, 1), ylim = c(0.0, 3.0))

plotParamDistr(mod1, mod2, mod3,
               paramName = "contrast_6C-6B[2]",
               xlim = c(-1.5, 1.5), ylim = c(0.0, 3.0))
