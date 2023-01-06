##' @include ProtModelFitMFVB.R
NULL

##' @include optimizeProteinRun.R
NULL

##' @include nuDistribution.R
NULL



# useful functions

blockDiag <- function(M1,M2){
  block_A_B <- rbind(cbind(M1,matrix(0,nrow=nrow(M1),ncol=ncol(M2))),
                     cbind(matrix(0,nrow=nrow(M2),ncol=ncol(M1)),M2))
  block_A_B
}

tr <- function (M){sum(diag(M))}





printIteration <- function(l, MAlpha, MBeta, MGamma, NObs, N, modelMiss, NMAR, 
                           nRdm, muqAlphaBeta, SigmaqAlphaBeta, eqinvsigma2eps, 
                           eqinvsigma2rdm, muqgammagy, Sigmaqgammagy, muqY, 
                           SigmaqY, mua, muqa,
                           ELBOValue, ELBOIndex){ 
  cat("Iteration ", l, "\n")
  cat("Alpha = ", round(muqAlphaBeta[1:MAlpha], 4), "\n")
  if(MAlpha > 1){
    cat("sigma2(Alpha) = ", round(sqrt(diag(SigmaqAlphaBeta)[1:MAlpha]), 4), "\n")
  } else {
    cat("sigma2(Alpha) = ", round(sqrt(SigmaqAlphaBeta), 4), "\n")
  }
  
  cat("eqinvsigma2eps ", round(eqinvsigma2eps, 4), "\n")
  if(nRdm > 0)
  {
    cat("Beta = ", round(muqAlphaBeta[(MAlpha+1):(MAlpha+MBeta)], 4), "\n")
    cat("eqinvsigma2rdm = ", round(eqinvsigma2rdm,4), "\n")
  }
  if(modelMiss){
    cat("Gamma = ", round(muqgammagy,4), "\n")
    if(MGamma + ifelse(NMAR, 1, 0) > 1){
      cat("sigma2(Gamma) = ", round(sqrt(diag(Sigmaqgammagy)), 4), "\n")
    } else {
      cat("sigma2(Gamma) = ", round(sqrt(Sigmaqgammagy), 4), "\n")
    }
    
    cat("Yobs average = ", round(mean(muqY[1:NObs]), 4), "\n")
    cat("Ymis average = ", round(mean(muqY[(NObs+1):N]), 4), "\n")
    cat("Trace of SigmaYmis = ", round(tr(SigmaqY[(NObs+1):N, (NObs+1):N]), 2), "\n")
    cat("mua_obs average = ", round(mean(mua[1:NObs]), 4), "\n")
    cat("mua_mis average = ", round(mean(mua[(NObs+1):N]), 4), "\n")
    cat("Aobs average = ", round(mean(muqa[1:NObs]), 4), "\n")
    cat("Amis average = ", round(mean(muqa[(NObs+1):N]), 4), "\n")
    
  }
  cat("ELBO = ", round(ELBOValue[ELBOIndex], 0), "\n")
  cat("\n")
}


printIterationHP <- function(l, gammaHierarch, sigmaHierarch, sigmaRdmHierarch,
                             modelMiss, nRdm, nuApproxi,
                             eq0_mu_gammagy, sigma2qmu0gammagy, 
                             eqinvsigma20gammagy, BqSigma20gammagy, eqNu0EpsBy2,
                             sigma2Nu0Eps, eqs02Eps, nuEpsBy2coef, BqNu0EpsBy2, 
                             CNu0EpsBy2, eqLogs02Eps, sumAllEqLogSigma2Eps, 
                             sumAllEqInvSigma2Eps, eqNu0RdmBy2, sigma2Nu0Rdm, 
                             eqs02Rdm, nuRdmBy2coef, BqNu0RdmBy2, CNu0RdmBy2,
                             GlobalELBOValue, ELBOIndex){
                               
  cat("Iteration ", l, "\n")
  if(modelMiss && gammaHierarch){
    cat("Mean(muGammaGy) = ", round(eq0_mu_gammagy,2), "\n")
    cat("Sd(muGammaGy) = ", round(sqrt(sigma2qmu0gammagy),2), "\n")
    cat("Mean(tauGammaGy) = ", round(eqinvsigma20gammagy,2),"\n")
    cat("Sd(tauGammaGy) = ", round(sqrt(eqinvsigma20gammagy/BqSigma20gammagy),2), "\n")
  }
  if(sigmaHierarch){
    cat("Mean(nu0Eps) = ", 2 * round(eqNu0EpsBy2, 2), "\n")
    if(nuApproxi == "normal"){
      cat("sigma2Nu0Eps = ", round(sigma2Nu0Eps, 2), "\n")
    } else if(nuApproxi == "gamma"){
      cat("BqNu0EpsBy2 = ", round(BqNu0EpsBy2, 2), "\n")
    } else {
      cat("CNu0EpsBy2 = ", round(CNu0EpsBy2, 2), "\n")
    }
    cat("nuEpsBy2coef = ", round(nuEpsBy2coef, 2), "\n")
    cat("Mean(s02Eps) = ", round(eqs02Eps, 2), "\n")
    cat("eqLogs02Eps = ", round(eqLogs02Eps, 2), "\n")
    cat("sumAllEqLogSigma2Eps = ", round(sumAllEqLogSigma2Eps, 2), "\n")
    cat("sumAllEqInvSigma2Eps = ", round(sumAllEqInvSigma2Eps, 2), "\n")
  }
  if(sigmaRdmHierarch && nRdm > 0){
    cat("Mean(nu0Rdm) = ", 2 * round(eqNu0RdmBy2, 2), "\n")
    cat("nuRdmBy2coef = ", round(nuRdmBy2coef, 2), "\n")
    if(nuApproxi == "normal"){
      cat("sigma2Nu0Rdm = ", round(sigma2Nu0Rdm, 2), "\n")
    } else if (nuApproxi == "gamma") {
      cat("BqNu0RdmBy2 = ", round(BqNu0RdmBy2, 2), "\n")
    } else {
      cat("CNu0RdmBy2 = ", round(CNu0RdmBy2, 2), "\n")
    }
    cat("Mean(s02Rdm) = ", round(eqs02Rdm, 2), "\n")
  }
  cat("Calculated ELBO : ", round(GlobalELBOValue[ELBOIndex], 4), "\n")
  if(ELBOIndex>1){
    cat("Previous ELBO : ", round(GlobalELBOValue[ELBOIndex-1], 4), "\n")
    cat("Relative diff : ", 
        abs(GlobalELBOValue[ELBOIndex] - GlobalELBOValue[ELBOIndex-1]) /
          abs(GlobalELBOValue[ELBOIndex]), "\n")
  }
}


# used in hierarchical version
updateOneProtParams <- function(lst,
                                modelSpecPars, 
                                algoArgs){
  
  distribStyle <- algoArgs$distribStyle
  if(is.null(distribStyle)){
    distribStyle <- "NonFactorized"
  } 
  
  taskVerbose <- algoArgs$taskVerbose
  if(is.null(taskVerbose)){
    taskVerbose <- FALSE
  }
  
  t_likelihood = t_likelihood(modelSpecPars)
  nu = t_fixedDF(modelSpecPars)
  modelMiss = modelMiss(modelSpecPars)
  NMAR = NMAR(modelSpecPars)
  sigmaHierarch = sigmaHierarch(modelSpecPars)
  sigmaRdmHierarch = sigmaRdmHierarch(modelSpecPars)
  gammaHierarch = gammaHierarch(modelSpecPars)
  
  if(taskVerbose){
    cat("*******************************************************************************\n")
    cat("Task starting updating parameters for protein ", lst$k, "(Nb obs = ", lst$N, ")\n")
    cat("*******************************************************************************\n")
  }
  
  start_time <- Sys.time() 
  lst <- within(lst,{
    if(distribStyle == "NonFactorized"){
       
      if(modelMiss){
        # Step 1 : update Alpha, Beta, yMis
        if(NMAR){
          downRightBlockSigmaqAlphaBetaYmis <- 
            diag(Sigmaqgammagy[MGamma+1, MGamma+1] + 
                   muqgammagy[MGamma+1] * muqgammagy[MGamma+1] + 
                   eqinvsigma2eps * eqinvb[(NObs+1):N])
          
        } else {
          downRightBlockSigmaqAlphaBetaYmis <- 
            diag(eqinvsigma2eps * eqinvb[(NObs+1):N])
        }
        upRightBlockSigmaqAlphaBetaYmis <- 
          - eqinvsigma2eps * t(WMis) %*% eqD[(NObs+1):N, (NObs+1):N]
        
        tUpRightBlockSigmaqAlphaBetaYmis <- t(upRightBlockSigmaqAlphaBetaYmis)
        
        invSigmaqAlphaBetaYmis[1:(MAlpha+MBeta), 1:(MAlpha+MBeta)] <-
          eqinvsigma2eps * WTEqDW + muqInvSigmaAlphaBeta
        invSigmaqAlphaBetaYmis[1:(MAlpha+MBeta), 
                               (MAlpha+MBeta+1):(MAlpha+MBeta+NMis)] <-
          upRightBlockSigmaqAlphaBetaYmis
        invSigmaqAlphaBetaYmis[(MAlpha+MBeta+1):(MAlpha+MBeta+NMis), 
                               1:(MAlpha+MBeta)] <-
          tUpRightBlockSigmaqAlphaBetaYmis
        invSigmaqAlphaBetaYmis[(MAlpha+MBeta+1):(MAlpha+MBeta+NMis),
                               (MAlpha+MBeta+1):(MAlpha+MBeta+NMis)] <-
          downRightBlockSigmaqAlphaBetaYmis
        
        SigmaqAlphaBetaYmis <- solve(invSigmaqAlphaBetaYmis)
        SigmaqAlphaBeta <- SigmaqAlphaBetaYmis[1:(MAlpha+MBeta),
                                               1:(MAlpha+MBeta)]
        
        SigmaqYmis <- SigmaqAlphaBetaYmis[(MAlpha+MBeta+1):(MAlpha+MBeta+NMis), 
                                          (MAlpha+MBeta+1):(MAlpha+MBeta+NMis)]
        SigmaqY[(NObs+1):N,(NObs+1):N] <- 
          SigmaqYmis
        
        SigmaqAlphaBetaY[1:(MAlpha+MBeta), 1:(MAlpha+MBeta)] <-
          SigmaqAlphaBeta
        SigmaqAlphaBetaY[1:(MAlpha+MBeta), 
                         (MAlpha+MBeta+NObs+1):(MAlpha+MBeta+N)] <-
          SigmaqAlphaBetaYmis[1:(MAlpha+MBeta), 
                              (MAlpha+MBeta+1):(MAlpha+MBeta+NMis)]
        SigmaqAlphaBetaY[(MAlpha+MBeta+NObs+1):(MAlpha+MBeta+N), 
                         1:(MAlpha+MBeta)] <-
          SigmaqAlphaBetaYmis[(MAlpha+MBeta+1):(MAlpha+MBeta+NMis),
                              1:(MAlpha+MBeta)]
        SigmaqAlphaBetaY[(MAlpha+MBeta+NObs+1):(MAlpha+MBeta+N),
                         (MAlpha+MBeta+NObs+1):(MAlpha+MBeta+N)] <-
          SigmaqYmis
        
        muqAlphaBetaYmis <- eqinvsigma2eps * t(WObs) %*% eqD[1:NObs,1:NObs] %*% yObs + 
          muqInvSigmaAlphaBeta %*% prior_mu_alpha_beta
        if(NMAR){
          if(MGamma > 1){
            muqAlphaBetaYmis <- 
              c(muqAlphaBetaYmis,
                muqgammagy[MGamma+1] * 
                  ( muqa[(NObs+1):N] - TT[(NObs+1):N,] %*% muqgammagy[1:MGamma] + 
                      muqgammagy[MGamma+1] * yRef))
          } else {
            muqAlphaBetaYmis <- 
              c(muqAlphaBetaYmis,
                muqgammagy[MGamma+1] * 
                  ( muqa[(NObs+1):N] - TT[(NObs+1):N,] * muqgammagy[1:MGamma] + 
                      muqgammagy[MGamma+1] * yRef))
          }
          
        } else {
          muqAlphaBetaYmis <- 
            c(muqAlphaBetaYmis,
              rep(0., NMis))
        }
        muqAlphaBetaYmis <- SigmaqAlphaBetaYmis %*% muqAlphaBetaYmis
        
        muqAlphaBeta <- muqAlphaBetaYmis[1:(MAlpha+MBeta)] 
        muqYMis <- muqAlphaBetaYmis[(MAlpha+MBeta+1):(MAlpha+MBeta+NMis)]
        muqY[(NObs+1):N] <- muqYMis
        
        
        if(NMAR){
          TmuqY <- cbind(TT, muqY - yRef)
          MatZeroMGamma[MGamma+1,MGamma+1] <- tr(SigmaqYmis)
        } # otherwise no need to update these objects as they are fixed 
        
        USigmaqAlphaBetaYUT <- U %*% SigmaqAlphaBetaY %*% t(U)
        
        # Step 2 : update gamma
        Sigmaqgammagy <- solve(t(TmuqY) %*% TmuqY + MatZeroMGamma + 
                                 eq0_inv_Sigma_gammagy)
        muqgammagy <- Sigmaqgammagy %*% 
          (t(TmuqY) %*% muqa + eq0_inv_Sigma_gammagy %*% eq0_mu_gammagy)
        
        # Step 3 : update a
        mua <- TmuqY %*% muqgammagy
        phimua <- dnorm(mua)
        Phia <- pnorm(TwoRMinus1*mua)
        muqa <- mua + TwoRMinus1*phimua/Phia
      } else { # we do not model missingness
        # Step 1 : update Alpha, Beta and NO step 2,3,4
        SigmaqAlphaBeta <- solve(eqinvsigma2eps * WTEqDW + muqInvSigmaAlphaBeta)
        
        muqAlphaBeta <- SigmaqAlphaBeta %*% 
          (eqinvsigma2eps * t(WObs) %*% eqD[1:NObs,1:NObs] %*% yObs + 
             muqInvSigmaAlphaBeta %*% prior_mu_alpha_beta)
        
        WSigmaqAlphaBetaWT <- W %*% SigmaqAlphaBeta %*% t(W)
      }
      
      #WmuqAlphaBeta <- W %*% muqAlphaBeta
      
      muqYMinusWmuqAlphaBeta <- muqY - W %*% muqAlphaBeta
      
      # Step 4 : update sigma2_epsilon
      if(sigmaHierarch){
        Aqeps <- eqNu0EpsBy2 + 0.5 * N
        prior_B_eps <- eqNu0EpsBy2 * eqs02Eps
      }
      if(modelMiss){
        Bqeps <- as.vector(
          prior_B_eps + 0.5 *
            (t(muqYMinusWmuqAlphaBeta) %*% eqD %*% muqYMinusWmuqAlphaBeta +
               tr(eqD %*% USigmaqAlphaBetaYUT)))
      } else {
        Bqeps <- as.vector(
          prior_B_eps + 0.5 *
            (t(muqYMinusWmuqAlphaBeta) %*% eqD %*% muqYMinusWmuqAlphaBeta +
               tr(eqD %*% WSigmaqAlphaBetaWT)))
      }
      # cat("SSE = ", t(muqYMinusWmuqAlphaBeta) %*% eqD %*% muqYMinusWmuqAlphaBeta, "\n")
      # cat("2*Bqeps = ", 2*(Bqeps-prior_B_eps), "\n")
      # cat("\n")
      # 
      eqinvsigma2eps <- Aqeps / Bqeps
      if(sigmaHierarch){
        eqlogsigma2eps <- log(Bqeps) - digamma(Aqeps)
      }
      
      # Step 6 : update sigma2_rdm
      if(nRdm>0){
        if(sigmaRdmHierarch){
          Aqrdm <- eqNu0RdmBy2 + 0.5 * MBetaPerRdmEffect
          prior_B_rdm <- eqNu0RdmBy2 * eqs02Rdm
        }
        for(r in 1:nRdm){
          start_index_r <- MAlpha + 1
          if(r>1) start_index_r <- start_index_r + 
              sum(MBetaPerRdmEffect[1:(r-1)])
          index_r <- seq(start_index_r, by=1, length.out = MBetaPerRdmEffect[r])
          
          Bqrdm[r] <- as.vector(prior_B_rdm[r] + 0.5 * 
                                  (t(muqAlphaBeta[index_r]) %*% muqAlphaBeta[index_r] +
                                     tr(SigmaqAlphaBeta[index_r, index_r])))
        }
        eqinvsigma2rdm <- Aqrdm / Bqrdm
        muqInvSigmaAlphaBeta <- blockDiag(prior_inv_sigma2_alpha_mat, 
                                          diag(rep(eqinvsigma2rdm, 
                                                   MBetaPerRdmEffect)))
        if(sigmaRdmHierarch){
          eqlogsigma2rdm <- log(Bqrdm) - digamma(Aqrdm)
        }
      }
      
      # Step 7 : update bi's
      if(t_likelihood){
        if(modelMiss){
          Bqb <- as.vector(
            prior_B_b + 0.5 * eqinvsigma2eps *
              (muqYMinusWmuqAlphaBeta * muqYMinusWmuqAlphaBeta + 
                 diag(USigmaqAlphaBetaYUT))
          )
        } else {
          Bqb <- as.vector(
            prior_B_b + 0.5 * eqinvsigma2eps *
              (muqYMinusWmuqAlphaBeta * muqYMinusWmuqAlphaBeta + 
                 diag(WSigmaqAlphaBetaWT))
          )
        }
        eqinvb <- Aqb/Bqb
        eqD <- diag(eqinvb)
        WTEqDW <- t(W) %*% eqD %*% W
        
      }
      # end of update parameters part for non factorized distributions style
        
    } else {
      
      # TO DO...
    }
    
  })
  
  end_time <- Sys.time() 
  
  if(taskVerbose){
    cat("*****************************************************************************************\n")
    cat("Task ending updating params for protein ", lst$k, "; elapsed = ", end_time - start_time,"\n")
    cat("*****************************************************************************************\n")
  }
  
  lst
  
}

# used in hierarchical version
updateOneProtELBO <- function(lst, modelSpecPars, algoArgs,
                              iter){
  nIterMax <- algoArgs$nIterMax
  refresh <- algoArgs$refresh
  
  distribStyle <- algoArgs$distribStyle
  if(is.null(distribStyle)){
    distribStyle <- "NonFactorized"
  } 
  
  taskVerbose <- algoArgs$taskVerbose
  if(is.null(taskVerbose)){
    taskVerbose <- FALSE
  }
  
  t_likelihood = t_likelihood(modelSpecPars)
  nu = t_fixedDF(modelSpecPars)
  modelMiss = modelMiss(modelSpecPars)
  NMAR = NMAR(modelSpecPars)
  sigmaHierarch = sigmaHierarch(modelSpecPars)
  sigmaRdmHierarch = sigmaRdmHierarch(modelSpecPars)
  gammaHierarch = gammaHierarch(modelSpecPars)
  
  ## START COMPUTING ELBO
  if(taskVerbose){
    cat("**************************************************************************\n")
    cat("Task starting computing ELBO for protein ", lst$k, "(Nb obs = ", lst$N, ")\n")
    cat("**************************************************************************\n")
  }
  start_time <- Sys.time()  
  lst <- within(lst,{
    ELBOIndex <- ELBOIndex + 1
    #ELBOIter[ELBOIndex] <- iter
    
    if(distribStyle == "NonFactorized")
    {
      
      if(modelMiss){
        ## ELBO STEP 1
        ELBODetails[ELBOIndex,1] <- 
          sum(log(pnorm(mua[1:NObs]))) +
          sum(log(rep(1,NMis)-pnorm(mua[(NObs+1):(NObs+NMis)]))) - 0.5 *
          tr(t(TmuqY) %*% TmuqY %*% Sigmaqgammagy)
        
        if(NMAR){
          ELBODetails[ELBOIndex,1] <- ELBODetails[ELBOIndex,1] - 0.5 *
            tr(SigmaqY) * ( muqgammagy[MGamma+1] * muqgammagy[MGamma+1] + 
                              Sigmaqgammagy[MGamma+1,MGamma+1] )
        } 
        
        ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] + 
          ELBODetails[ELBOIndex,1]
        
        # ELBO STEP 2
        ELBODetails[ELBOIndex,2] <- 
          0.5 * determinant(Sigmaqgammagy, logarithm = TRUE)$modulus -
          0.5 * (t(muqgammagy - eq0_mu_gammagy ) %*% 
                   eq0_inv_Sigma_gammagy %*% 
                   (muqgammagy - eq0_mu_gammagy)) -
          0.5 * tr(eq0_inv_Sigma_gammagy %*% Sigmaqgammagy)
        ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] + ELBODetails[ELBOIndex,2]
      }
      
      # ELBO STEP 3
      if(modelMiss){
        ELBODetails[ELBOIndex,3] <- - 0.5 * eqinvsigma2eps *
          (t(muqYMinusWmuqAlphaBeta) %*% eqD %*% muqYMinusWmuqAlphaBeta +
             tr(eqD %*% U %*% SigmaqAlphaBetaY %*% t(U))) - 0.5 *
          (t(muqAlphaBeta - prior_mu_alpha_beta) %*% muqInvSigmaAlphaBeta %*%
             (muqAlphaBeta - prior_mu_alpha_beta) +
             tr(muqInvSigmaAlphaBeta %*% SigmaqAlphaBeta)) +
          0.5 * determinant(SigmaqAlphaBetaYmis, logarithm = TRUE)$modulus
      } else {
        ELBODetails[ELBOIndex,3] <- - 0.5 * eqinvsigma2eps *
          (t(muqYMinusWmuqAlphaBeta) %*% eqD %*% muqYMinusWmuqAlphaBeta +
             tr(eqD %*% W %*% SigmaqAlphaBeta %*% t(W))) - 0.5 *
          (t(muqAlphaBeta - prior_mu_alpha_beta) %*% muqInvSigmaAlphaBeta %*%
             (muqAlphaBeta - prior_mu_alpha_beta) +
             tr(muqInvSigmaAlphaBeta %*% SigmaqAlphaBeta)) +
          0.5 * determinant(SigmaqAlphaBeta, logarithm = TRUE)$modulus
      }
      ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] + ELBODetails[ELBOIndex,3]
      
      # ELBO STEP 4
      if(nRdm > 0){
        if(sigmaRdmHierarch){
          ELBODetails[ELBOIndex,4] <- 
            sum(eqinvsigma2rdm * Bqrdm) -
            sum(Aqrdm * log(Bqrdm)) +
            sum(lgamma(Aqrdm))
        } else {
          ELBODetails[ELBOIndex,4] <- 
            sum(eqinvsigma2rdm * (Bqrdm - prior_B_rdm)) -
            sum(Aqrdm * log(Bqrdm))
        }
        
        ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] + ELBODetails[ELBOIndex,4]
      }
      
      # ELBO STEP 5
      if(sigmaHierarch){
        ELBODetails[ELBOIndex,5] <- 
          eqinvsigma2eps * Bqeps - 
          Aqeps * log(Bqeps) +
          lgamma(Aqeps)
        
      } else {
        ELBODetails[ELBOIndex,5] <- 
          eqinvsigma2eps * (Bqeps - prior_B_eps) -
          Aqeps * log(Bqeps)
      }
      
      ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] + ELBODetails[ELBOIndex,5]
      
      # ELBO STEP 6
      if(t_likelihood > 0){
        ELBODetails[ELBOIndex,6] <- 
          sum(eqinvb * (Bqb - prior_B_b)) -
          sum(Aqb * log(Bqb))
        ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] + ELBODetails[ELBOIndex,6]
      }
      
      ## END ELBO COMPUTATION
      
      # end of ELBO update part for non factorized distributions style
    } else {
      # TO DO...
    }
    
    # if(iter %% refresh == 0 || iter == nIterMax){
    #   printIteration(iter, MAlpha, MBeta, MGamma, NObs, N, modelMiss, NMAR,
    #                  nRdm, muqAlphaBeta, SigmaqAlphaBeta, eqinvsigma2eps, 
    #                  eqinvsigma2rdm, muqgammagy, Sigmaqgammagy, muqY, 
    #                  SigmaqY, mua, muqa,
    #                  ELBOValue, ELBOIndex)
    # }
  
  })
  
  end_time <- Sys.time()  
  
  if(taskVerbose){
    cat("****************************************************************************************\n")
    cat("Task ending computing ELBO for protein ", lst$k, "; elapsed = ", end_time - start_time,"\n")
    cat("****************************************************************************************\n")
  }
  
  lst
}

# used for hierarchical version of model
##' @export
runMFVBAllProts <- function(x, InputData, bp){
  nIterMax <- algoArgs(x)$nIterMax
  if(is.null(nIterMax)) stop("missing in algoArgs : [nIterMax]")
  relTol <- algoArgs(x)$relTol
  if(is.null(relTol)) stop("missing in algoArgs : [relTol]")
  refresh <- algoArgs(x)$refresh
  if(is.null(refresh)) stop("missing in algoArgs : [refresh]")
  ELBO_each <- algoArgs(x)$ELBO_each
  if(is.null(ELBO_each)) stop("missing in algoArgs : [ELBO_each]")
  
  distribStyle <- algoArgs(x)$distribStyle
  if(is.null(distribStyle)){
    distribStyle <- "NonFactorized"
  } else if(!distribStyle %in% c("NonFactorized","Factorized")){
    msg <- "wrong distribStyle passed to the method : authorized styles are "
    msg <- paste0(msg, "[NonFactorized], [Factorized]")
    stop(msg)
  }
  
  if(distribStyle == "Factorized")
    stop("Factorized version of MFVB not implemented with hierarchical model")
  
  # filterOutput <- algoArgs(x)$filterOutput
  # if(is.null(filterOutput)){
  #   filterOutput = FALSE
  # }
  
  initTauValue <- algoArgs(x)$initTauValue
  if(is.null(initTauValue)){
    initTauValue = 1.
  }
  
  nuApproxi <- algoArgs(x)$nuApproxi
  if(is.null(nuApproxi)){
    nuApproxi == "gamma"
  }
  if(nuApproxi != "none" && nuApproxi != "norm" && nuApproxi != "gamma"){
    msg <- paste0("Unrecognized [nuApproxi] in algoArgs : received [", nuApproxi)
    msg <- paste0(msg, "], should be in ([none], [norm], [gamma])")
    stop(msg)
  }
  
  dumpWorkspace <- algoArgs(x)$dumpWorkspace
  if(is.null(dumpWorkspace)){
    dumpWorkspace <- FALSE
  }
  
  t_likelihood = t_likelihood(ModelSpecPars(x))
  nu = t_fixedDF(ModelSpecPars(x))
  modelMiss = modelMiss(ModelSpecPars(x))
  NMAR = NMAR(ModelSpecPars(x))
  sigmaHierarch = sigmaHierarch(ModelSpecPars(x))
  sigmaRdmHierarch = sigmaRdmHierarch(ModelSpecPars(x))
  gammaHierarch = gammaHierarch(ModelSpecPars(x))
  
  nRdm <- InputData$nRdm
  
  nProt <- InputData$K
  nELBOComputationMax <- nIterMax %/% ELBO_each + 2
  
  # in case of parallel computation, try to optimize the order of proteins
  protRunOrder <- 1:nProt
  protRunIndex <- 1:nProt
  if(!is.null(bp)){
    tasks <- bp$tasks
    if(tasks == 0){
      tasks = min(c(nProt, bp$workers * 3))
    }
    res <- optimizeProteinRun(InputData$P,
                              InputData$NMis,
                              tasks
                              )
    protRunOrder <- res$runOrder
    protRunIndex <- res$runIndex
    # cat("\n protein run order : \n")
    # print(protRunOrder)
    # cat("\n protein run index : \n")
    # print(protRunIndex)
  }
  
  start_time = Sys.time()
  
  # big memory storage for results of all intermediary computation, by protein
  cat("Creating workspace for each protein results...\n")
  lst <- list()
  
  # ELBO storage at global level
  ELBOnColCommon <- 6
  ELBOnColUnique <- 6
  GlobalELBOValue <- rep(0., nELBOComputationMax)
  GlobalELBOIter <- rep(0., nELBOComputationMax)
  GlobalELBODetails <- matrix(0., nrow = nELBOComputationMax, 
                              ncol = ELBOnColCommon + ELBOnColUnique)
  
  for(k in 1:nProt){
    lst[[protRunIndex[k]]] <- list()
    lst[[protRunIndex[k]]]$k <- k
    
    lst[[protRunIndex[k]]] <- within(data = lst[[protRunIndex[k]]], {
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
      
      res <- list()
      res$paramDistr <- list()
      res$fitResults <- list()
      
      # init timing results
      res$fitResults$elapsed <- 0.
      
      # init ELBO storage
      ELBOValue <- rep(0., nELBOComputationMax)
      #ELBOIter <- rep(0, nELBOComputationMax)
      ELBODetails <- matrix(0.,
                            nrow = nELBOComputationMax,
                            ncol = ELBOnColCommon)
      
      # convert data into the right format
      NObs <- InputData$NObs[k]
      if(NObs<2){
        msg <- paste0("Number of observations for protein ", k, " too low (<2)")
        stop(msg)
      }
      NMis <- InputData$NMis[k]
      N <- ifelse(modelMiss, NObs+NMis, NObs)
      MAlpha <- InputData$MAlpha[k]
      nRdm = InputData$nRdm
      MBeta <- InputData$MBeta[k]
      MBetaPerRdmEffect <- 
        InputData$MBetaPerRdmEffect[spos:(spos+nRdm-1)]
      MGamma <- InputData$MGamma[k]
      yObs <- InputData$yObs[seq(from = yobspos,
                                 length.out = NObs)]
      XObs <- matrix(
        InputData$XObs[seq(from = xobspos,
                           length.out = NObs * MAlpha)],
        nrow = NObs, byrow = FALSE)
      
      if(nRdm > 0){
        ZObs <- matrix(
          InputData$ZObs[seq(from = zobspos,
                             length.out = NObs * MBeta)],
          nrow = NObs, byrow = FALSE)
      }
      
      TObs <- matrix(
        InputData$TObs[seq(from = tobspos,
                           length.out = NObs * MGamma)],
        nrow = NObs, byrow = FALSE)
      XMis <- matrix(
        InputData$XMis[seq(from = xmispos,
                           length.out = NMis * MAlpha)],
        nrow = NMis, byrow = FALSE)
      if(nRdm > 0){
        ZMis <- matrix(
          InputData$ZMis[seq(from = zmispos,
                             length.out = NMis * MBeta)],
          nrow = NMis, byrow = FALSE)
        
      }
      
      TMis <- matrix(
        InputData$TMis[seq(from = tmispos,
                           length.out = NMis * MGamma)],
        nrow = NMis, byrow = FALSE)
      
      yRef <- InputData$refIntensity[k]
      if(modelMiss){
        X <- rbind(XObs, XMis)
        TT <- rbind(TObs, TMis)
        if(nRdm > 0){
          Z <- rbind(ZObs, ZMis)
        }
      } else {
        X <- XObs
        TT <- TObs
        if(nRdm > 0){
          Z <- ZObs
        }
      }
      
      if(nRdm > 0){
        W <- cbind(X, Z)
        WObs <- cbind(XObs, ZObs)
        WMis <- cbind(XMis, ZMis)
      } else {
        W <- X
        WObs <- XObs
        WMis <- XMis
      }
      if(modelMiss){
        U <- cbind(-W,
                   diag(c(rep(0., NObs),
                          rep(1., NMis))))
      } else {
        U <- cbind(-W,
                   diag(rep(0., NObs)))
      }
      
      
      # init priors
      prior_mu_alpha <- c(priorIntensityMean(ModelSpecPars(x)),
                          rep(priorFoldChangeMean(ModelSpecPars(x)),
                              MAlpha-1))
      prior_mu_alpha_beta <- c(prior_mu_alpha,
                               rep(0., MBeta))
      prior_sigma_alpha <- c(priorIntensitySd(ModelSpecPars(x)),
                             rep(priorFoldChangeSd(ModelSpecPars(x)),
                                 MAlpha-1))
      prior_inv_sigma2_alpha <- 1/
        (prior_sigma_alpha * prior_sigma_alpha)
      if( MAlpha > 1){
        prior_inv_sigma2_alpha_mat <-
          diag(prior_inv_sigma2_alpha)
      } else {
        prior_inv_sigma2_alpha_mat <-
          matrix(prior_inv_sigma2_alpha)
      }
      
      prior_mu_gamma <- priorGammaInterceptMean(ModelSpecPars(x))
      prior_sigma_gamma <- priorGammaInterceptSd(ModelSpecPars(x))
      prior_mu_gamma_y <- priorGammaYMean(ModelSpecPars(x))
      prior_sigma_gamma_y <- priorGammaYSd(ModelSpecPars(x))
      prior_inv_sigma2_gamma_y <-
        1/(prior_sigma_gamma_y * prior_sigma_gamma_y)
      if(MGamma > 1){
        # prior_mu_gamma <- c(prior_mu_gamma,
        #                     rep(priorGammaTermsMean(ModelSpecPars(x)), MGamma-1))
        # prior_sigma_gamma <- c(prior_sigma_gamma,
        #                        rep(priorGammaTermsSd(ModelSpecPars(x)), MGamma-1))
        prior_mu_gamma <- c(prior_mu_gamma,
                            priorGammaTermsMean(ModelSpecPars(x)))
        prior_sigma_gamma <- c(prior_sigma_gamma,
                               priorGammaTermsSd(ModelSpecPars(x)))
      }
      prior_mu_gammagy <- prior_mu_gamma
      prior_inv_sigma2_gammagy <- 1/
        (prior_sigma_gamma * prior_sigma_gamma)
      if(NMAR){
        prior_mu_gammagy <- c(prior_mu_gammagy,
                              prior_mu_gamma_y)
        prior_inv_sigma2_gammagy <-
          c(prior_inv_sigma2_gammagy,
            1/(prior_sigma_gamma_y * prior_sigma_gamma_y))
      }
      if(MGamma + NMAR > 1){
        prior_inv_Sigma_gammagy <- diag(prior_inv_sigma2_gammagy)
      } else {
        prior_inv_Sigma_gammagy <-
          matrix(prior_inv_sigma2_gammagy, ncol = 1, nrow = 1)
      }
      
      prior_nu_eps <- priorNuEpsilon(ModelSpecPars(x))
      prior_scale_eps <- priorScaleEpsilon(ModelSpecPars(x))
      prior_A_eps <- 0.5 * prior_nu_eps
      prior_B_eps <- 0.5 * prior_nu_eps *
        prior_scale_eps * prior_scale_eps
      
      prior_nu_rdm <- priorNuRdm(ModelSpecPars(x))
      prior_scale_rdm <- priorScaleRdm(ModelSpecPars(x))
      prior_A_rdm <- 0.5 * prior_nu_rdm
      prior_B_rdm <-
        rep(0.5 * prior_nu_rdm * prior_scale_rdm *
              prior_scale_rdm,
            nRdm)
      
      # buid missingness data vector
      
      R <- c(rep(1., NObs), rep(0., NMis))
      TwoRMinus1 <- 2*R-1
      
      # initialization step
      muqgammagy <- vector()
      Sigmaqgammagy <- matrix()
      
      if(NMAR){
        muqgammagy <- rep(0., MGamma+1)
        Sigmaqgammagy <-
          matrix(rep(0., (MGamma+1)*(MGamma+1)),
                 ncol = MGamma+1)
        Sigmaqgammagy[MGamma+1, MGamma+1] <- 1.
      } else {
        muqgammagy <- rep(0., MGamma)
        Sigmaqgammagy <- matrix(rep(0., MGamma * MGamma,
                                    ncol = MGamma))
      }
      
      mua <- rep(0., N)
      muqa <- 2 * TwoRMinus1 * dnorm(0.)
      
      eqinvb <- rep(1., N)
      eqD <- diag(eqinvb)
      
      # note we initialized tau here!
      eqinvsigma2eps <- initTauValue
      eqlogsigma2eps <- log(eqinvsigma2eps)
      eqinvsigma2rdm <- rep(1., nRdm)
      eqlogsigma2rdm <- log(eqinvsigma2rdm)
      muqInvSigmaAlphaBeta <-
        blockDiag(prior_inv_sigma2_alpha_mat,
                  diag(rep(eqinvsigma2rdm, MBetaPerRdmEffect)))
      
      # hyperparameters
      eqNu0EpsBy2 <- sqrt(nProt) #arbitrary
      eqs02Eps <- 1.
      eqNu0RdmBy2 <- rep(sqrt(nProt), nRdm)
      eqs02Rdm <- rep(1., nRdm)
      if(gammaHierarch){
        if(MGamma + NMAR > 1){
          eq0_inv_Sigma_gammagy <- diag(rep(1., MGamma + NMAR))
        } else {
          eq0_inv_Sigma_gammagy <-
            matrix(1., ncol = 1, nrow = 1)
        }
        eq0_mu_gammagy <- rep(0., MGamma + NMAR)
      } else {
        eq0_inv_Sigma_gammagy <- prior_inv_Sigma_gammagy
        eq0_mu_gammagy <- prior_mu_gammagy
      }
      
      
      # (following init part is mainly for memory allocation, as the
      # so initialized quantities will be overwritten straight away. However some
      # zeros initialized in the var/covar matrices  will still remain and be used
      muqAlpha <- prior_mu_alpha
      muqBeta <- rep(0., MBeta)
      muqAlphaBeta <- c(muqAlpha, muqBeta)
      
      muqYMis <- rep(0, NMis)
      muqY <- vector()
      # sigma2qYMis <- rep(0, NMis)
      # sigma2qY <- vector()
      if(modelMiss){
        muqY <- c(yObs, muqYMis)
        #sigma2qY <- c(rep(0, NObs), sigma2qYMis)
      } else {
        muqY <- yObs
        #sigma2qY <- rep(0, NObs)
      }
      if(NMAR){
        TmuqY <- cbind(TT, muqY)
      } else {
        TmuqY <- TT
      }
      
      muqAlphaBetaYmis <- c(muqAlphaBeta, muqYMis)
      muqAlphaBetaY <- c(muqAlphaBeta, yObs,
                         muqYMis)
      
      SigmaqAlphaBeta <-
        matrix(rep(0.,(MAlpha+MBeta)*
                     (MAlpha+MBeta)),
               ncol = (MAlpha+MBeta))
      SigmaqAlphaBetaYmis <-
        matrix(rep(0., (MAlpha+MBeta+NMis) *
                     (MAlpha+MBeta+NMis)),
               nrow = MAlpha + MBeta + NMis,
               ncol = MAlpha + MBeta + NMis)
      invSigmaqAlphaBetaYmis <-
        matrix(rep(0., (MAlpha+MBeta+NMis) *
                     (MAlpha+MBeta+NMis)),
               nrow = MAlpha + MBeta + NMis,
               ncol = MAlpha + MBeta + NMis)
      SigmaqAlphaBetaY <-
        matrix(rep(0., (MAlpha+MBeta+N) *
                     (MAlpha+MBeta+N)),
               nrow = MAlpha + MBeta + N,
               ncol = MAlpha + MBeta + N)
      SigmaqYmis <-
        matrix(rep(0., NMis * NMis),
               nrow = NMis, ncol = NMis)
      
      SigmaqY <-
        matrix(rep(0., N * N),
               nrow = N, ncol = N)
      
      
      # other quantities precomputation or memory allocation
      Aqeps <- prior_A_eps + 0.5*N
      Aqrdm <- prior_A_rdm + 0.5*MBetaPerRdmEffect
      Bqrdm <- rep(0., nRdm)
      prior_A_b <- nu/2; prior_B_b <- nu/2
      Aqb <- prior_A_b + 0.5
      Bqb <- rep(0., N)
      
      phimua <- rep(0., N)
      Phia <- rep(0., N)
      
      WmuqAlphaBeta <- W %*% muqAlphaBeta
      upRightBlockSigmaqAlphaBetaYmis <-
        matrix(nrow = MAlpha + MBeta,
               ncol = NMis)
      tUpRightBlockSigmaqAlphaBetaYmis <-
        matrix(nrow = NMis,
               ncol = MAlpha + MBeta)
      downRightBlockSigmaqAlphaBetaYmis <-
        matrix(nrow = NMis, ncol = NMis)
      WSigmaqAlphaBetaWT <-
        W %*% SigmaqAlphaBeta %*% t(W)
      WTEqDW <- t(W) %*% eqD %*% W
      
      MatZeroGamma <- matrix()
      if(NMAR){
        MatZeroMGamma <-
          matrix(rep(0., (MGamma+1) * (MGamma+1)),
                 ncol = MGamma+1)
      } else {
        MatZeroMGamma <-
          matrix(rep(0., MGamma * MGamma),
                 ncol = MGamma)
      }
      
      muqYMinusWmuqAlphaBeta <- rep(0., N)
      muqaMinusTmuqYmuqGamma <- rep(0., N)
      USigmaqAlphaBetaYUT <-
        U %*% SigmaqAlphaBetaY %*% t(U)
      
      # prepare coefficients for contrasts
      nContrast <- InputData$nContrast[k]
      contrastCoefs <- list()
      if(nContrast > 0){
        for(i in 1:nContrast){
          contrastCoefs[[i]] <- 
            InputData$contrastCoefs[(conpos+(i-1)*MAlpha):(conpos+i*MAlpha-1)]
        }
      }
      
    })
  } # end loop on proteins (k)
  
  cat("Creating workspace for each protein results... DONE!\n")
  
  # hyperparameters - global storage (if hierarchical set-up applies)
  cat("Creating workspace for hyperparameters... \n")
  HPLst <- list()
  HPLst <- within(data = HPLst, expr = {
    ELBOValue <- rep(0., nELBOComputationMax)
    ELBODetails <- matrix(0.,
                          nrow = nELBOComputationMax,
                          ncol = ELBOnColCommon + ELBOnColUnique)
    
    
    MGamma <- InputData$MGamma[1]
    nRdm <- InputData$nRdm
    eqNu0EpsBy2 <- sqrt(nProt) #arbitrary
    nuEpsBy2coef <- 0.
    CNu0EpsBy2 <- 0.
    sigma2Nu0EpsBy2 <- 0.
    AqNu0EpsBy2 <- 0.
    BqNu0EpsBy2 <- 0.
    eqs02Eps <- 1.
    eqLogs02Eps <- 0.
    Aqs02Eps <- 0.
    Bqs02Eps <- 0.
    
    eqNu0RdmBy2 <- rep(sqrt(nProt), nRdm)
    nuRdmBy2coef <- rep(0., nRdm)
    CNu0RdmBy2 <- rep(0., nRdm)
    sigma2Nu0RdmBy2 <- rep(0., nRdm)
    AqNu0RdmBy2 <- rep(0., nRdm)
    BqNu0RdmBy2 <- rep(0., nRdm)
    eqs02Rdm <- rep(1., nRdm)
    eqLogs02Rdm <- rep(0., nRdm)
    Aqs02Rdm <- rep(0., nRdm)
    Bqs02Rdm <- rep(0., nRdm)
    if(modelMiss && gammaHierarch){
      if(MGamma + NMAR > 1){
        eq0_inv_Sigma_gammagy <- diag(rep(1., MGamma + NMAR))
      } else {
        eq0_inv_Sigma_gammagy <-
          matrix(1., ncol = 1, nrow = 1)
      }
      eq0_mu_gammagy <- rep(0., MGamma + NMAR)
    } 
    sigma2qmu0gammagy <- rep(1., MGamma + NMAR)
    eqinvsigma20gammagy <- rep(1., MGamma + NMAR)
    BqSigma20gammagy <- rep(0., MGamma + NMAR)
    
    AllMuqGammagy <- matrix(0., nrow = nProt, ncol = MGamma + NMAR)
    AllSigma2qGammagy <- matrix(0., nrow = nProt, ncol = MGamma + NMAR)
    AllEqInvSigma2Eps <- rep(0., nProt)
    AllEqLogSigma2Eps <- rep(0., nProt)
    if(nRdm > 0){
      AllEqInvSigma2Rdm <- matrix(0., nrow = nProt, ncol = nRdm)
      AllEqLogSigma2Rdm <- matrix(0., nrow = nProt, ncol = nRdm)
    }
  })
  
  
  cat("Creating workspace for hyperparameters... DONE!\n")
  
  # loop niter iterations
  # if(k == 10){
  #   bbb <- 1
  # }
  #browser()
  converged <- FALSE
  ELBOIndex <- 0
  for(k in 1:nProt){
    lst[[k]]$ELBOIndex <- 0
  }
  
  start_time <- Sys.time()
  tryCatch(
    error = function(cnd) {
      print(conditionMessage(cnd))
    },
    {
      for (l in 1:nIterMax){
  
        cat("\n")
        cat("Starting iteration ", l, "...\n")
        cat("*****************************\n\n")
        
        cat("Updating parameters for all proteins ...\n")
        
        if(!is.null(bp)){
          lst <- bpmapply(FUN = updateOneProtParams,
                          lst = lst,
                          MoreArgs = list(modelSpecPars = ModelSpecPars(x),
                                          algoArgs = algoArgs(x)),
                          SIMPLIFY = FALSE,
                          BPPARAM = bp)
        } else {
          lst <- mapply(FUN = updateOneProtParams, 
                        lst = lst,
                        MoreArgs = list(modelSpecPars = ModelSpecPars(x),
                                        algoArgs = algoArgs(x)),
                        SIMPLIFY = FALSE)
        }
        
        cat("Updating parameters for all proteins ... DONE!\n")
        
        
        cat("Updating hyperparameters ...\n")
        
        # create local storage for all protein relevant parameters for ease of computation
        for(k in 1:nProt){
          if(modelMiss & gammaHierarch){
            HPLst$AllMuqGammagy[k, ] <- lst[[k]]$muqgammagy
            if(HPLst$MGamma + NMAR > 1){
              HPLst$AllSigma2qGammagy[k,] <- diag(lst[[k]]$Sigmaqgammagy)
            } else {
              HPLst$AllSigma2qGammagy[k,] <- lst[[k]]$Sigmaqgammagy
            }
          }
          if(sigmaHierarch){
            HPLst$AllEqInvSigma2Eps[k] <- lst[[k]]$eqinvsigma2eps
            HPLst$AllEqLogSigma2Eps[k] <- lst[[k]]$eqlogsigma2eps
          }
          if(nRdm > 0 && sigmaRdmHierarch){
            HPLst$AllEqInvSigma2Rdm[k,] <- lst[[k]]$eqinvsigma2rdm
            HPLst$AllEqLogSigma2Rdm[k,] <- lst[[k]]$eqlogsigma2rdm
          }
        }
        
        # hyperparameters update (per se)
        HPLst <- within(data = HPLst, expr = {
          if(modelMiss && gammaHierarch){
            # STEP 7 : update mu0 of gamma coefficients
            eq0_mu_gammagy <- colMeans(AllMuqGammagy)
            sigma2qmu0gammagy <- 1/(nProt * eqinvsigma20gammagy)
            if(MGamma + NMAR > 1){
              eq0_inv_Sigma_gammagy <- diag(sigma2qmu0gammagy)
            } else {
              eq0_inv_Sigma_gammagy <- matrix(sigma2qmu0gammagy, nrow = 1, ncol = 1)
            }
            # STEP 8 : update sigma20 of gamma coefficients
            for(i in 1:(MGamma+NMAR)){
              BqSigma20gammagy[i] <- 0.5 * 
                (sum((AllMuqGammagy[,i] - eq0_mu_gammagy[i])^2) +
                   sum(AllSigma2qGammagy[,i]) + nProt * sigma2qmu0gammagy[i])
            }
            eqinvsigma20gammagy <- (0.5*nProt - 1) / BqSigma20gammagy
          }
          
          if(sigmaHierarch){
            # STEP 9 : update s20 epsilon
            Aqs02Eps <- nProt * eqNu0EpsBy2 + 1
            Bqs02Eps <- eqNu0EpsBy2 * sum(AllEqInvSigma2Eps)
            eqs02Eps <- Aqs02Eps / Bqs02Eps
            eqLogs02Eps <- digamma(Aqs02Eps) - log(Bqs02Eps)
            
            # STEP 10 : update nu0 epsilon
            sumAllEqLogSigma2Eps <- sum(AllEqLogSigma2Eps)
            sumAllEqInvSigma2Eps <- sum(AllEqInvSigma2Eps)
            nuEpsBy2coef <- nProt * eqLogs02Eps - sumAllEqLogSigma2Eps -
              eqs02Eps * sumAllEqInvSigma2Eps
            if(nuApproxi == "norm" || nuApproxi == "gamma"){
              eqNu0EpsBy2 <- modeANu(A = 1, 
                                     xlogxcoef = nProt,
                                     lgammaxcoef = -nProt,
                                     xcoef = nuEpsBy2coef,
                                     const = 0.)
              sigma2Nu0EpsBy2 <- varianceANuApproxi(A = 1,
                                                    xlogxcoef = nProt,
                                                    lgammaxcoef = -nProt,
                                                    xcoef = nuEpsBy2coef,
                                                    const = 0.)
              
              if (nuApproxi == "gamma"){
                #AqNu0EpsBy2 <- nProt/2 + 1
                AqNu0EpsBy2 <- 1 + eqNu0EpsBy2^2/sigma2Nu0EpsBy2
                #BqNu0EpsBy2 <- -nuEpsBy2coef - nProt
                BqNu0EpsBy2 <- eqNu0EpsBy2 / sigma2Nu0EpsBy2
                eqNu0EpsBy2 <- AqNu0EpsBy2 / BqNu0EpsBy2
              } 
            } else {
              multConstant <- integrate(f = dNu, lower = 0, upper = Inf,
                                        xlogxcoef = nProt,
                                        lgammaxcoef = -nProt,
                                        xcoef = nuEpsBy2coef,
                                        const = 0)$value
              CNu0EpsBy2 <- -log(multConstant)
              eqNu0EpsBy2 <- integrate(f = nudNu, lower = 0, upper = Inf, 
                                       xlogxcoef = nProt,
                                       lgammaxcoef = -nProt,
                                       xcoef = nuEpsBy2coef,
                                       const = CNu0EpsBy2)$value
            }
            
          }
          
          
          if(nRdm > 0 && sigmaRdmHierarch){
            # STEP 11 : update s20 rdm
            Aqs02Rdm <- nProt * eqNu0RdmBy2 + 1
            Bqs02Rdm <- eqNu0RdmBy2 * colSums(AllEqInvSigma2Rdm)
            eqs02Rdm <- Aqs02Rdm / Bqs02Rdm
            eqLogs02Rdm <- digamma(Aqs02Rdm) - log(Bqs02Rdm)
            
            # STEP 12 : update nu0 rdm
            sumAllEqLogSigma2Rdm <- colSums(AllEqLogSigma2Rdm)
            sumAllEqInvSigma2Rdm <- colSums(AllEqInvSigma2Rdm)
            for(r in 1:nRdm){
              nuRdmBy2coef[r] <- nProt * eqLogs02Rdm[r] - sumAllEqLogSigma2Rdm[r] -
                eqs02Rdm[r] * sumAllEqInvSigma2Rdm[r]
              if(nuApproxi == "norm" || nuApproxi == "gamma"){
                eqNu0RdmBy2[r] <- modeANu(A = 1, 
                                          xlogxcoef = nProt,
                                          lgammaxcoef = -nProt,
                                          xcoef = nuRdmBy2coef[r],
                                          const = 0.)
                sigma2Nu0RdmBy2[r] <- varianceANuApproxi(A = 1,
                                                        xlogxcoef = nProt,
                                                        lgammaxcoef = -nProt,
                                                        xcoef = nuRdmBy2coef[r],
                                                        const = 0.)
                if(nuApproxi == "gamma"){
                  #AqNu0RdmBy2[r] <- nProt/2 + 1
                  AqNu0RdmBy2[r] <- 1 + eqNu0RdmBy2[r]^2 / sigma2Nu0RdmBy2[r]
                  #BqNu0RdmBy2[r] <- -nuRdmBy2coef[r] - nProt
                  BqNu0RdmBy2[r] <- eqNu0RdmBy2[r] / sigma2Nu0RdmBy2[r]
                  eqNu0RdmBy2[r] <- AqNu0RdmBy2[r] / BqNu0RdmBy2[r]
                } 
              } else {
                multConstant <- integrate(f = dNu, lower = 0, upper = Inf,
                                          xlogxcoef = nProt,
                                          lgammaxcoef = -nProt,
                                          xcoef = nuRdmBy2coef[r],
                                          const = 0)$value
                CNu0RdmBy2[r] <- -log(multConstant)
                eqNu0RdmBy2[r] <- integrate(f = nudNu, lower = 0, upper = Inf,
                                            xlogxcoef = nProt,
                                            lgammaxcoef = -nProt,
                                            xcoef = nuRdmBy2coef[r],
                                            const = CNu0RdmBy2[r])$value
              }
              
            }
          }
        })
        
        # send back new value of hyperparameters to all protein workspaces
        for(k in 1:nProt){
          if(modelMiss && gammaHierarch){
            lst[[k]]$eq0_mu_gammagy <- HPLst$eq0_mu_gammagy
            lst[[k]]$eq0_inv_Sigma_gammagy <- HPLst$eq0_inv_Sigma_gammagy
          }
          if(sigmaHierarch){
            lst[[k]]$eqNu0EpsBy2 <- HPLst$eqNu0EpsBy2
            lst[[k]]$eqs02Eps <- HPLst$eqs02Eps
            
          }
          if(sigmaRdmHierarch && nRdm > 0){
            lst[[k]]$eqNu0RdmBy2 <- HPLst$eqNu0RdmBy2
            lst[[k]]$eqs02Rdm <- HPLst$eqs02Rdm
          }
        }
        
        cat("Updating hyperparameters ... DONE!\n")
        
        if(l == 1 || (l %% ELBO_each) == 0 || l == nIterMax){
    
          cat("Computing ELBO for all proteins ...\n")
          ELBOIndex <- ELBOIndex+1
        
          GlobalELBOIter[ELBOIndex] <- l
          
          if(!is.null(bp)){
            lst <- bpmapply(FUN = updateOneProtELBO,
                            lst = lst,
                            MoreArgs = list(modelSpecPars = ModelSpecPars(x),
                                            algoArgs = algoArgs(x),
                                            iter = l),
                            SIMPLIFY = FALSE,
                            BPPARAM = bp)
          } else {
            lst <- mapply(FUN = updateOneProtELBO, 
                          lst = lst,
                          MoreArgs = list(modelSpecPars = ModelSpecPars(x),
                                          algoArgs = algoArgs(x),
                                          iter = l),
                          SIMPLIFY = FALSE)
          }
          
          for(k in 1:nProt){
            GlobalELBOValue[ELBOIndex] <- GlobalELBOValue[ELBOIndex] +
              lst[[k]]$ELBOValue[ELBOIndex]
            for(i in 1:(ELBOnColCommon)){
              GlobalELBODetails[ELBOIndex,i] <- GlobalELBODetails[ELBOIndex,i] +
                lst[[k]]$ELBODetails[ELBOIndex,i]
            }
          }
          
          cat("Computing ELBO for all proteins ... DONE!\n")
          
          cat("Computing ELBO for hyper parameters part ...\n")
          
          HPLst <- within(data = HPLst, expr = {
            if(modelMiss && gammaHierarch){
              # STEP 7 : mu0gammagy
              ELBODetails[ELBOIndex,7] <- -0.5 *
                sum(log(sigma2qmu0gammagy))
              ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] +
                ELBODetails[ELBOIndex,7]

              # STEP 8 : sigma20gammagy
              ELBODetails[ELBOIndex,8] <-
                sum(eqinvsigma20gammagy * BqSigma20gammagy) -
                (0.5*nProt-1) * sum(log(BqSigma20gammagy))
              ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] +
                ELBODetails[ELBOIndex,8]
            }

            if(sigmaHierarch){
              # STEP 9 : nu0 epsilon
              if(nuApproxi == "norm"){
                ELBODetails[ELBOIndex,9] <-
                  0.5 * log(sigma2Nu0EpsBy2) +
                  (nProt/2) * (log(eqNu0EpsBy2) -
                                 0.5*sigma2Nu0EpsBy2/eqNu0EpsBy2/eqNu0EpsBy2) +
                  nProt * eqNu0EpsBy2 * eqLogs02Eps -
                  eqNu0EpsBy2 * eqs02Eps * sumAllEqInvSigma2Eps
              } else if(nuApproxi == "gamma"){
                ELBODetails[ELBOIndex,9] <-
                  eqNu0EpsBy2 * sumAllEqLogSigma2Eps -
                  AqNu0EpsBy2 * log(BqNu0EpsBy2) + 
                  lgamma(AqNu0EpsBy2)
              } else {
                ELBODetails[ELBOIndex,9] <-
                  eqNu0EpsBy2 * sumAllEqLogSigma2Eps - CNu0EpsBy2
              }

              ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] +
                ELBODetails[ELBOIndex,9]

              # STEP 10 : so2 epsilon
              ELBODetails[ELBOIndex,10] <-
                - (Aqs02Eps * log(Bqs02Eps) - lgamma(Aqs02Eps) +
                     (Aqs02Eps-1) * eqLogs02Eps - Bqs02Eps * eqs02Eps)
              ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] +
                ELBODetails[ELBOIndex,10]
            }

            if(sigmaRdmHierarch && nRdm > 0){
              # STEP 11 : nu0 rdm
              if(nuApproxi == "norm"){
                ELBODetails[ELBOIndex,11] <-
                  0.5 * sum(log(sigma2Nu0RdmBy2)) +
                  (nProt/2) * sum(log(eqNu0RdmBy2) -
                                 0.5*sigma2Nu0RdmBy2/eqNu0RdmBy2/eqNu0RdmBy2) +
                  nProt * sum(eqNu0RdmBy2 * eqLogs02Rdm) -
                  sum(eqNu0RdmBy2 * eqs02Rdm * sumAllEqInvSigma2Rdm)
              } else if(nuApproxi == "gamma"){
                ELBODetails[ELBOIndex,11] <-
                  sum(eqNu0RdmBy2 * sumAllEqLogSigma2Rdm) -
                  sum(AqNu0RdmBy2 * log(BqNu0RdmBy2)) + 
                  sum(lgamma(AqNu0RdmBy2))
              } else {
                ELBODetails[ELBOIndex,11] <- sum(
                  eqNu0RdmBy2 * sumAllEqLogSigma2Rdm - CNu0RdmBy2)
              }

              ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] +
                ELBODetails[ELBOIndex,11]

              # STEP 12 : so2 rdm
              ELBODetails[ELBOIndex,12] <-
                - (Aqs02Rdm * log(Bqs02Rdm) - lgamma(Aqs02Rdm) +
                     (Aqs02Rdm-1) * eqLogs02Rdm - Bqs02Rdm * eqs02Rdm)
              ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] +
                ELBODetails[ELBOIndex,12]
            }
          })

          GlobalELBOValue[ELBOIndex] <- GlobalELBOValue[ELBOIndex] + 
            HPLst$ELBOValue[ELBOIndex]
          for(i in (ELBOnColCommon+1):(ELBOnColCommon + ELBOnColUnique)){
            GlobalELBODetails[ELBOIndex,i] <- HPLst$ELBODetails[ELBOIndex,i] 
          }
        
          
          cat("Computing ELBO for hyper parameters part ... DONE!\n")
          
          if(l %% refresh == 0 || l == nIterMax){
            printIterationHP(l, gammaHierarch, sigmaHierarch, sigmaRdmHierarch,
                             modelMiss, nRdm, nuApproxi,
                             HPLst$eq0_mu_gammagy, HPLst$sigma2qmu0gammagy, 
                             HPLst$eqinvsigma20gammagy, HPLst$BqSigma20gammagy, HPLst$eqNu0EpsBy2,
                             HPLst$sigma2Nu0Eps, HPLst$eqs02Eps, HPLst$nuEpsBy2coef, HPLst$BqNu0EpsBy2, 
                             HPLst$CNu0EpsBy2, HPLst$eqLogs02Eps, HPLst$sumAllEqLogSigma2Eps, 
                             HPLst$sumAllEqInvSigma2Eps, HPLst$eqNu0RdmBy2, HPLst$sigma2Nu0Rdm, 
                             HPLst$eqs02Rdm, HPLst$nuRdmBy2coef, HPLst$BqNu0RdmBy2, HPLst$CNu0RdmBy2,
                             GlobalELBOValue, ELBOIndex)
            if(l %% refresh == 0 && dumpWorkspace){
              cat("Dumping full workspace... \n")
              saveRDS(object <- lst, 
                      file <- paste0("ProtWorkspaceDump_It", l,".rds"))
              saveRDS(object <- HPLst, 
                      file <- paste0("HPWorkspaceDump_It", l,".rds"))
              cat("Dumping full workspace... DONE !\n")
            }
          }
          
          if(ELBOIndex>1){
            if(abs(GlobalELBOValue[ELBOIndex] - GlobalELBOValue[ELBOIndex-1]) <
               relTol * abs(GlobalELBOValue[ELBOIndex])){
              converged <- TRUE
              break
            }
          }
        } # end computing ELBO
    
      }# end loop on l (iteration)
    }) # body of tryCatch
    
  end_time = Sys.time()

  # store general fit results
  x@fitResults <- list( nIter = l,
                        converged = converged,
                        nELBO = ELBOIndex,
                        ELBOIter = GlobalELBOIter[1:ELBOIndex],
                        ELBOValue = GlobalELBOValue[1:ELBOIndex],
                        ELBODetails = GlobalELBODetails[1:ELBOIndex,],
                        elapsed = end_time - start_time)
  
  

  # update marginal distributions with current protein results
  

  # start with per protein parameters
  for(k in 1:nProt){
    lst[[k]] <- within(lst[[k]],{
      paramDistr <- list()
      # normal distributions for alpha
      newParamDistrib <- list()
      for(i in 1:MAlpha){
        mui <- muqAlphaBeta[i]
        if(MAlpha + MBeta > 1){
          sdi <- sqrt(SigmaqAlphaBeta[i,i])
        } else {
          sdi <- sqrt(SigmaqAlphaBeta[1])
        }
        newParamDistrib[[i]] <- list(main = TRUE,
                                     mean = mui,
                                     sd = sdi,
                                     type = "norm",
                                     args = list(mean = mui,
                                                 sd = sdi))
      }

      names(newParamDistrib) <- paste0("alpha[",
                                       seq(from = apos, to = apos+MAlpha-1),
                                       "]")
      paramDistr <- newParamDistrib
      
      # normal distributions for contrasts on alpha's
      if(nContrast > 0){
        newParamDistrib <- list()
        for(i in 1:nContrast){
          mui <- sum(muqAlphaBeta[1:MAlpha]*contrastCoefs[[i]])
          if(MAlpha + MBeta > 1){
            vari <-  t(contrastCoefs[[i]]) %*% SigmaqAlphaBeta[1:MAlpha,1:MAlpha] %*% 
              contrastCoefs[[i]]
          } else {
            vari <- SigmaqAlphaBeta[1] * contrastCoefs[[i]][1] * contrastCoefs[[i]][1]
          }
          sdi <- sqrt(vari)
          newParamDistrib[[i]] <- list(main = TRUE,
                                       mean = mui,
                                       sd = sdi,
                                       type = "norm",
                                       args = list(mean = mui,
                                                   sd = sdi))
        }
        
        if(k == 1){
          start <- 1
        } else {
          start <- sum(InputData$nContrast[1:(k-1)]) + 1
        }
        names(newParamDistrib) <- 
          paste0("contrast[", 
                 seq(from = start, to = start+nContrast-1),
                 "]")
        
        paramDistr <- c(paramDistr, newParamDistrib)
      }
      
      # normal distributions for beta
      if(nRdm > 0){
        newParamDistrib <- list()
        for(i in 1:MBeta){
          mui <- muqAlphaBeta[MAlpha+i]
          if(MAlpha + MBeta > 1){
            sdi <- sqrt(SigmaqAlphaBeta[MAlpha+i, MAlpha+i])
          } else {
            sdi <- sqrt(SigmaqAlphaBeta[1])
          }
          newParamDistrib[[i]] <- list(main = TRUE,
                                       mean = mui,
                                       sd = sdi,
                                       type = "norm",
                                       args = list(mean = mui,
                                                   sd = sdi))
        }

        names(newParamDistrib) <- paste0("beta[",
                                         seq(from = bpos, to = bpos+MBeta-1),
                                         "]")
        paramDistr <- c(paramDistr, newParamDistrib)
      }

      # normal distributions for gamma and gammaY
      if(modelMiss){
        newParamDistrib <- list()
        for(i in 1:MGamma){
          mui <- muqgammagy[i]
          if(MGamma > 1){
            sdi <- sqrt(Sigmaqgammagy[i, i])
          } else {
            sdi <- sqrt(Sigmaqgammagy[1])
          }
          newParamDistrib[[i]] <- list(main = TRUE,
                                       mean = mui,
                                       sd = sdi,
                                       type = "norm",
                                       args = list(mean = mui,
                                                   sd = sdi))
        }

        distrNames <- paste0("gamma[",
                             seq(from = cpos, to = cpos+MGamma-1),
                             "]")

        if(NMAR){
          gammay_name <- paste0("gammY[", k, "]")
          mui <- muqgammagy[MGamma+1]
          sdi <- sqrt(Sigmaqgammagy[MGamma+1,MGamma+1])
          newParamDistrib[[MGamma+1]] <- list(main = TRUE,
                                              mean = mui,
                                              sd = sdi,
                                              type = "norm",
                                              args = list(mean = mui,
                                                          sd = sdi))
          distrNames <- c(distrNames, gammay_name)
        }
        names(newParamDistrib) <- distrNames
        paramDistr <- c(paramDistr, newParamDistrib)

      }

      # normal distributions for Ymis
      #if(!filterOutput && modelMiss){
      if(modelMiss)
      {
        newParamDistrib <- list()
        for(i in 1:NMis){
          mui <- muqY[NObs+i]
          sdi <- sqrt(SigmaqY[NObs+i,NObs+i])

          newParamDistrib[[i]] <- list(main = FALSE,
                                       mean = mui,
                                       sd = sdi,
                                       type = "norm",
                                       args = list(mean = mui,
                                                   sd = sdi))
        }

        names(newParamDistrib) <- paste0("yMis[",
                                         seq(from = ymispos, to = ymispos+NMis-1),
                                         "]")

        paramDistr <- c(paramDistr, newParamDistrib)
      }


      # gamma distributions for tau(eps)
      tau_name <- paste0("tauEps[", k, "]")
      newParamDistrib <- list()
      newParamDistrib[[1]] <- list(main = TRUE,
                                   mean = Aqeps/Bqeps,
                                   sd = sqrt(Aqeps/Bqeps/Bqeps),
                                   type = "gamma",
                                   args = list(shape = Aqeps,
                                               rate = Bqeps))
      names(newParamDistrib) <- c(tau_name)
      paramDistr <- c(paramDistr, newParamDistrib)
      
      # inverse gamma distributions for sigma2(eps)
      sigma2_name <- paste0("sigma2Eps[", k, "]")
      mui <- Bqeps/(Aqeps-1)
      sdi <- sqrt(Bqeps*Bqeps/((Aqeps-1)^2*(Aqeps-2)))
      newParamDistrib <- list()
      newParamDistrib[[1]] <- list(main = TRUE,
                                   mean = mui,
                                   sd = sdi,
                                   type = "invgamma",
                                   args = list(alpha = Aqeps,
                                               beta = Bqeps))
      names(newParamDistrib) <- c(sigma2_name)
      paramDistr <- c(paramDistr, newParamDistrib)

      # gamma distributions for tau_rdm and inverse gamma distributions for 
      # sigma2_rdm
      if(nRdm > 0){
        newParamDistrib <- list()
        for(i in 1:nRdm){
          newParamDistrib[[i]] <- list(main = TRUE,
                                       mean = Aqrdm[i]/Bqrdm[i],
                                       sd = sqrt(Aqrdm[i]/Bqrdm[i]/Bqrdm[i]),
                                       type = "gamma",
                                       args = list(shape = Aqrdm[i],
                                                   rate = Bqrdm[i]))
        }

        names(newParamDistrib) <- paste0("tauRdm[",
                                         seq(from = spos, to = spos+nRdm-1),
                                         "]")

        paramDistr <- c(paramDistr, newParamDistrib)
        
        newParamDistrib <- list()
        for(i in 1:nRdm){
          alpha <- Aqrdm[i]
          beta <- Bqrdm[i]
          mui <- beta/(alpha-1)
          sdi <- sqrt(beta*beta/((alpha-1)^2*(alpha-2)))
          newParamDistrib[[i]] <- list(main = TRUE,
                                       mean = mui,
                                       sd = sdi,
                                       type = "invgamma",
                                       args = list(alpha = alpha,
                                                   beta = beta))
        }
        
        names(newParamDistrib) <- paste0("sigma2Rdm[",
                                         seq(from = spos, to = spos+nRdm-1),
                                         "]")
        
        paramDistr <- c(paramDistr, newParamDistrib)
      }
      
      

      # gamma distribution for inv(b)
      #if(!filterOutput && t_likelihood > 0){
      if(t_likelihood > 0){
        # inv(b) parameters for observed data
        newParamDistrib <- list()
        for(i in 1:NObs){
          newParamDistrib[[i]] <- list(main = FALSE,
                                       mean = Aqb/Bqb[i],
                                       sd = sqrt(Aqb/Bqb[i]/Bqb[i]),
                                       type = "gamma",
                                       args = list(shape = Aqb,
                                                   rate = Bqb[i]))
        }
        names(newParamDistrib) <- paste0("invBObs[",
                                         seq(from = yobspos, to = yobspos+NObs-1),
                                         "]")

        paramDistr <- c(paramDistr, newParamDistrib)

        # inv(b) parameters for missing data
        if(modelMiss){
          newParamDistrib <- list()
          for(i in 1:NMis){
            newParamDistrib[[i]] <- list(main = FALSE,
                                         mean = Aqb/Bqb[NObs+i],
                                         sd = sqrt(Aqb/Bqb[NObs+i]/Bqb[NObs+i]),
                                         type = "gamma",
                                         args = list(shape = Aqb,
                                                     rate = Bqb[NObs+i]))
          }
          names(newParamDistrib) <- paste0("invBMis[",
                                           seq(from = ymispos, to = ymispos+NMis-1),
                                           "]")

          paramDistr <- c(paramDistr, newParamDistrib)
        }

      }

      # for 'a' parameter : truncated normals
      #if(!filterOutput && modelMiss){
      if(modelMiss){  
        # for observed responses
        newParamDistrib <- list()
        for(i in 1:NObs){
          ratioi <- dnorm(mua[i])/pnorm(mua[i])
          sdi <- sqrt(1 - mua[i]*ratioi - ratioi^2)
          newParamDistrib[[i]] <- list(main = FALSE,
                                       mean = muqa[i],
                                       sd = sdi,
                                       type = "normTrunc",
                                       args = list(mean = mua[i],
                                                   sd = 1,
                                                   min = 0.))
        }
        names(newParamDistrib) <- paste0("aObs[",
                                         seq(from = yobspos, to = yobspos+NObs-1),
                                         "]")

        paramDistr <- c(paramDistr, newParamDistrib)

        # for missing responses
        newParamDistrib <- list()
        for(i in 1:NMis){
          ratioi <- dnorm(mua[NObs+i])/(1-pnorm(mua[NObs+i]))
          sdi <- sqrt(1 + mua[NObs+i]*ratioi - ratioi^2)
          newParamDistrib[[i]] <- list(main = FALSE,
                                       mean = muqa[NObs+i],
                                       sd = sdi,
                                       type = "normTrunc",
                                       args = list(mean = mua[NObs+i],
                                                   sd = 1,
                                                   max = 0.))
        }
        names(newParamDistrib) <- paste0("aMis[",
                                         seq(from = ymispos, to = ymispos+NMis-1),
                                         "]")

        paramDistr <- c(paramDistr, newParamDistrib)
      }

    })
    # end storage of distributions for protein k
  }

  # now store all distributions in correct order
  for(k in 1:nProt){
    x@paramDistr <- c(x@paramDistr, lst[[protRunIndex[k]]]$paramDistr)
  }
  
  
  # now add hyperparameters
  HPLst <- within(data = HPLst, expr = {
    paramDistr = list()
    if(modelMiss && gammaHierarch){
      # normal distributions for mu0gammagy
      newParamDistrib <- list()
      for(i in 1:(MGamma+NMAR)){
        mui <- eq0_mu_gammagy[i]
        sdi <- sqrt(sigma2qmu0gammagy[i])
        newParamDistrib[[i]] <- list(main = TRUE,
                                     mean = mui,
                                     sd = sdi,
                                     type = "norm",
                                     args = list(mean = mui,
                                                 sd = sdi))
      }
      names(newParamDistrib) <- paste0("mu0gm[",
                                       seq(from = 1, to = MGamma+NMAR),
                                       "]")
      paramDistr <- c(paramDistr, newParamDistrib)

      # gamma distribution for tau0gammagy

      newParamDistrib <- list()
      for(i in 1:(MGamma+NMAR)){
        mui <- eqinvsigma20gammagy[i]
        sdi <- sqrt(eqinvsigma20gammagy[i]/BqSigma20gammagy[i])
        newParamDistrib[[i]] <- list(main = TRUE,
                                     mean = mui,
                                     sd = sdi,
                                     type = "gamma",
                                     args = list(shape = 0.5*nProt - 1,
                                                 rate = BqSigma20gammagy[i]))
      }
      names(newParamDistrib) <- paste0("tau0gm[",
                                       seq(from = 1, to = MGamma+NMAR),
                                       "]")
      paramDistr <- c(paramDistr, newParamDistrib)
      
      # inverse gamma distribution for sigma20gammagy
      
      newParamDistrib <- list()
      for(i in 1:(MGamma+NMAR)){
        alpha <- 0.5*nProt - 1
        beta <- BqSigma20gammagy[i]
        mui <- beta/(alpha-1)
        sdi <- sqrt(beta*beta/((alpha-1)^2*(alpha-2)))
        newParamDistrib[[i]] <- list(main = TRUE,
                                     mean = mui,
                                     sd = sdi,
                                     type = "invgamma",
                                     args = list(alpha = alpha,
                                                 beta = beta))
      }
      names(newParamDistrib) <- paste0("sigma20gm[",
                                       seq(from = 1, to = MGamma+NMAR),
                                       "]")
      paramDistr <- c(paramDistr, newParamDistrib)
    }

    if(sigmaHierarch){
      # gamma distribution for s0Eps2
      newParamDistrib <- list()
      newParamDistrib[[1]] <- list(main = TRUE,
                                   mean = Aqs02Eps/Bqs02Eps,
                                   sd = sqrt(Aqs02Eps/Bqs02Eps/Bqs02Eps),
                                   type = "gamma",
                                   args = list(shape = Aqs02Eps,
                                               rate = Bqs02Eps))
      names(newParamDistrib) <- c("s0Eps2[1]")
      paramDistr <- c(paramDistr, newParamDistrib)

      # nu special distribution for nu0Eps, or normal approxi, or gamma approxi
      newParamDistrib <- list()

      if(nuApproxi == "norm"){
        munu <- 2 * eqNu0EpsBy2
        sdnu <- 2 * sqrt(sigma2Nu0EpsBy2)
        newParamDistrib[[1]] <- list(main = TRUE,
                                     mean = munu,
                                     sd = sdnu,
                                     type = "norm",
                                     args = list(mean = munu,
                                                 sd = sdnu))
      } else if(nuApproxi == "gamma"){
        munu <- AqNu0EpsBy2/(BqNu0EpsBy2/2)
        sdnu <- sqrt(munu/(BqNu0EpsBy2/2))
        newParamDistrib[[1]] <- list(main = TRUE,
                                     mean = munu,
                                     sd = sdnu,
                                     type = "gamma",
                                     args = list(shape = AqNu0EpsBy2,
                                                 rate = BqNu0EpsBy2/2))
      } else {
        munu <- 2 * eqNu0EpsBy2
        sdnu <- 2 * sqrt( integrate(f = nu2dNu,
                                    lower = 0,
                                    upper = Inf,
                                    xlogxcoef = nProt,
                                    lgammaxcoef = -nProt,
                                    xcoef = nuEpsBy2coef,
                                    const = CNu0EpsBy2)$value
                          - eqNu0EpsBy2 * eqNu0EpsBy2)

        newParamDistrib[[1]] <- list(main = TRUE,
                                     mean = munu,
                                     sd = sdnu,
                                     type = "ANu",
                                     args = list(A = 2,
                                                 xlogxcoef = nProt,
                                                 lgammaxcoef = -nProt,
                                                 xcoef = nuEpsBy2coef,
                                                 const = CNu0EpsBy2))
      }


      names(newParamDistrib) <- c("nu0Eps[1]")
      paramDistr <- c(paramDistr, newParamDistrib)
    }

    if(sigmaRdmHierarch && nRdm > 0){
      # gamma distribution for s0Rdm2
      newParamDistrib <- list()
      for(r in 1:nRdm){

        newParamDistrib[[r]] <- list(main = TRUE,
                                     mean = Aqs02Rdm[r]/Bqs02Rdm[r],
                                     sd = sqrt(Aqs02Rdm[r]/Bqs02Rdm[r]/Bqs02Rdm[r]),
                                     type = "gamma",
                                     args = list(shape = Aqs02Rdm[r],
                                                 rate = Bqs02Rdm[r]))
      }

      names(newParamDistrib) <- paste0("s0Rdm2[",
                                       seq(from = 1, to = nRdm),
                                       "]")

      paramDistr <- c(paramDistr, newParamDistrib)

      # nu special distribution for nu0Rdm, or normal approxi, or gamma approxi
      newParamDistrib <- list()
      for(r in 1:nRdm){
        if(nuApproxi == "norm"){
          munu <- 2 * eqNu0RdmBy2[r]
          sdnu <- 2 * sqrt(sigma2Nu0RdmBy2[r])
          newParamDistrib[[r]] <- list(main = TRUE,
                                       mean = munu,
                                       sd = sdnu,
                                       type = "norm",
                                       args = list(mean = munu,
                                                   sd = sdnu))
        } else if(nuApproxi == "gamma"){
          munu <- AqNu0RdmBy2[r]/(BqNu0RdmBy2[r]/2)
          sdnu <- sqrt(munu/(BqNu0RdmBy2[r]/2))
          newParamDistrib[[r]] <- list(main = TRUE,
                                       mean = munu,
                                       sd = sdnu,
                                       type = "gamma",
                                       args = list(shape = AqNu0RdmBy2[r],
                                                   rate = BqNu0RdmBy2[r]/2))
        } else {
          munu <- 2 * eqNu0RdmBy2[r]
          sdnu <- 2 * sqrt( integrate(f = nu2dNu,
                                      lower = 0,
                                      upper = Inf,
                                      xlogxcoef = nProt,
                                      lgammaxcoef = -nProt,
                                      xcoef = nuRdmBy2coef[r],
                                      const = CNu0RdmBy2[r])$value
                            - eqNu0RdmBy2[r] * eqNu0RdmBy2[r])

          newParamDistrib[[r]] <- list(main = TRUE,
                                       mean = munu,
                                       sd = sdnu,
                                       type = "ANu",
                                       args = list(A = 2,
                                                   xlogxcoef = nProt,
                                                   lgammaxcoef = -nProt,
                                                   xcoef = nuRdmBy2coef[r],
                                                   const = CNu0RdmBy2[r]))
        }
      }

      names(newParamDistrib) <- paste0("nu0Rdm[",
                                       seq(from = 1, to = nRdm),
                                       "]")

      paramDistr <- c(paramDistr, newParamDistrib)
    }
  })

  x@paramDistr <- c(x@paramDistr, HPLst$paramDistr)
  
  x
  
}

runMFVBOneProt <- function(x, InputData, protIndex){
  nIterMax <- algoArgs(x)$nIterMax
  if(is.null(nIterMax)) stop("missing in algoArgs : [nIterMax]")
  relTol <- algoArgs(x)$relTol
  if(is.null(relTol)) stop("missing in algoArgs : [relTol]")
  refresh <- algoArgs(x)$refresh
  if(is.null(refresh)) stop("missing in algoArgs : [refresh]")
  ELBO_each <- algoArgs(x)$ELBO_each
  if(is.null(ELBO_each)) stop("missing in algoArgs : [ELBO_each]")
  
  distribStyle <- algoArgs(x)$distribStyle
  if(is.null(distribStyle)){
    distribStyle <- "NonFactorized"
  } else if(!distribStyle %in% c("NonFactorized","Factorized")){
    msg <- "wrong distribStyle passed to the method : authorized styles are "
    msg <- paste0(msg, "[NonFactorized], [Factorized]")
    stop(msg)
  }
  
  # filterOutput <- algoArgs(x)$filterOutput
  # if(is.null(filterOutput)){
  #   filterOutput = FALSE
  # }
  
  initTauValue <- algoArgs(x)$initTauValue
  if(is.null(initTauValue)){
    initTauValue = 1.
  }
  
  t_likelihood = t_likelihood(ModelSpecPars(x))
  nu = t_fixedDF(ModelSpecPars(x))
  modelMiss = modelMiss(ModelSpecPars(x))
  NMAR = NMAR(ModelSpecPars(x))
  sigmaHierarch = sigmaHierarch(ModelSpecPars(x))
  sigmaRdmHierarch = sigmaRdmHierarch(ModelSpecPars(x))
  
  
  nProt <- InputData$K
  
  k <- protIndex
  
  
  nELBOComputationMax <- nIterMax %/% ELBO_each + 2
  
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
  
  
  res <- list()
  res$paramDistr <- list()
  res$fitResults <- list()
  
  
  cat("Starting MFVB run for protein ", k, "\n")
  cat("********************************\n\n")
  
  start_time <- Sys.time()
  
  # init ELBO storage
  ELBOValue <- rep(0., nELBOComputationMax)
  ELBOIter <- rep(0., nELBOComputationMax)
  ELBODetails <- matrix(0., nrow = nELBOComputationMax, ncol = 7)
  
  # convert data into the right format
  
  NObs <- InputData$NObs[k]
  if(NObs<2){
    msg <- paste0("Number of observations for protein ", k, " too low (<2)")
    stop(msg)
  }
  NMis <- InputData$NMis[k]
  N <- ifelse(modelMiss, NObs+NMis, NObs)
  MAlpha <- InputData$MAlpha[k]
  nRdm = InputData$nRdm
  MBeta <- InputData$MBeta[k]
  MBetaPerRdmEffect <- InputData$MBetaPerRdmEffect[spos:(spos+nRdm-1)]
  MGamma <- InputData$MGamma[k]
  yObs <- InputData$yObs[seq(from = yobspos, length.out = NObs)]
  XObs <- matrix(InputData$XObs[seq(from = xobspos, length.out = NObs * MAlpha)],
                 nrow = NObs, byrow = FALSE)
  if(nRdm > 0){
    ZObs <- matrix(InputData$ZObs[seq(from = zobspos, length.out = NObs * MBeta)],
                   nrow = NObs, byrow = FALSE)
  }
  
  TObs <- matrix(InputData$TObs[seq(from = tobspos, length.out = NObs * MGamma)],
                 nrow = NObs, byrow = FALSE)
  XMis <- matrix(InputData$XMis[seq(from = xmispos, length.out = NMis * MAlpha)],
                 nrow = NMis, byrow = FALSE)
  if(nRdm > 0){
    ZMis <- matrix(InputData$ZMis[seq(from = zmispos, length.out = NMis * MBeta)],
                   nrow = NMis, byrow = FALSE)
    
  }
  
  TMis <- matrix(InputData$TMis[seq(from = tmispos, length.out = NMis * MGamma)],
                 nrow = NMis, byrow = FALSE)
  
  yRef <- InputData$refIntensity[k]
  if(modelMiss){
    X <- rbind(XObs, XMis)
    TT <- rbind(TObs, TMis)
    if(nRdm > 0){
      Z <- rbind(ZObs, ZMis)
    }
  } else {
    X <- XObs
    TT <- TObs
    if(nRdm > 0){
      Z <- ZObs
    }
  }
  
  if(nRdm > 0){
    W <- cbind(X,Z)
    WObs <- cbind(XObs, ZObs)
    WMis <- cbind(XMis, ZMis)
  } else {
    W <- X
    WObs <- XObs
    WMis <- XMis
  }
  if(modelMiss){
    U <- cbind(-W, diag(c(rep(0., NObs), rep(1., NMis))))
  } else {
    U <- cbind(-W, diag(rep(0., NObs)))
  }
  
  
  # init priors
  prior_mu_alpha <- c(priorIntensityMean(ModelSpecPars(x)), 
                      rep(priorFoldChangeMean(ModelSpecPars(x)), MAlpha-1))
  prior_mu_alpha_beta <- c(prior_mu_alpha, rep(0., MBeta))
  prior_sigma_alpha <- c(priorIntensitySd(ModelSpecPars(x)),
                         rep(priorFoldChangeSd(ModelSpecPars(x)), MAlpha-1))
  prior_inv_sigma2_alpha <- 1/(prior_sigma_alpha*prior_sigma_alpha)
  if( MAlpha > 1){
    prior_inv_sigma2_alpha_mat <- diag(prior_inv_sigma2_alpha)
  } else {
    prior_inv_sigma2_alpha_mat <- matrix(prior_inv_sigma2_alpha)
  }
  
  prior_mu_gamma <- priorGammaInterceptMean(ModelSpecPars(x))
  prior_sigma_gamma <- priorGammaInterceptSd(ModelSpecPars(x))
  prior_mu_gamma_y <- priorGammaYMean(ModelSpecPars(x))
  prior_sigma_gamma_y <- priorGammaYSd(ModelSpecPars(x))
  prior_inv_sigma2_gamma_y <- 1/(prior_sigma_gamma_y*prior_sigma_gamma_y)
  if(MGamma > 1){
    # prior_mu_gamma <- c(prior_mu_gamma, 
    #                     rep(priorGammaTermsMean(ModelSpecPars(x)), MGamma-1))
    # prior_sigma_gamma <- c(prior_sigma_gamma, 
    #                        rep(priorGammaTermsSd(ModelSpecPars(x)), MGamma-1))
    prior_mu_gamma <- c(prior_mu_gamma, 
                        priorGammaTermsMean(ModelSpecPars(x)))
    prior_sigma_gamma <- c(prior_sigma_gamma, 
                           priorGammaTermsSd(ModelSpecPars(x)))
  }
  prior_mu_gammagy <- prior_mu_gamma
  prior_inv_sigma2_gammagy <- 1/(prior_sigma_gamma * prior_sigma_gamma)
  if(NMAR){
    prior_mu_gammagy <- c(prior_mu_gammagy, prior_mu_gamma_y)
    prior_inv_sigma2_gammagy <- c(prior_inv_sigma2_gammagy, 
                                  1/(prior_sigma_gamma_y*prior_sigma_gamma_y))
  }
  if(MGamma + NMAR > 1){
    prior_inv_Sigma_gammagy <- diag(prior_inv_sigma2_gammagy)
  } else {
    prior_inv_Sigma_gammagy <- matrix(prior_inv_sigma2_gammagy,
                                      ncol = 1,
                                      nrow = 1)
  }
  
  prior_nu_eps <- priorNuEpsilon(ModelSpecPars(x))
  prior_scale_eps <- priorScaleEpsilon(ModelSpecPars(x))
  prior_A_eps <- 0.5 * prior_nu_eps
  prior_B_eps <- 0.5 * prior_nu_eps * prior_scale_eps * prior_scale_eps
  
  prior_nu_rdm <- priorNuRdm(ModelSpecPars(x))
  prior_scale_rdm <- priorScaleRdm(ModelSpecPars(x))
  prior_A_rdm <- 0.5 * prior_nu_rdm
  prior_B_rdm <- rep(0.5 * prior_nu_rdm * prior_scale_rdm * prior_scale_rdm,
                     nRdm)
  
  # buid missingness data vector
  
  R <- c(rep(1., NObs), rep(0., NMis))
  TwoRMinus1 <- 2*R-1
  
  # initialization step
  muqgammagy <- vector()
  Sigmaqgammagy <- matrix()
  
  if(NMAR){
    muqgammagy <- rep(0., MGamma+1)
    Sigmaqgammagy <- matrix(rep(0., (MGamma+1)*(MGamma+1)), 
                            ncol = MGamma+1)
    Sigmaqgammagy[MGamma+1, MGamma+1] <- 1.
  } else {
    muqgammagy <- rep(0., MGamma)
    Sigmaqgammagy <- matrix(rep(0., MGamma * MGamma, 
                                ncol = MGamma))
  }
  
  mua <- rep(0., N)
  muqa <- 2 * TwoRMinus1 * dnorm(0.)
  
  eqinvb <- rep(1., N)
  eqD <- diag(eqinvb)
  
  # change init value here!
  eqinvsigma2eps <- initTauValue
  eqinvsigma2rdm <- rep(1., nRdm)
  muqInvSigmaAlphaBeta <- blockDiag(prior_inv_sigma2_alpha_mat, 
                                    diag(rep(eqinvsigma2rdm, MBetaPerRdmEffect)))
  
  
  # (following init part is mainly for memory allocation, as the 
  # so initialized quantities will be overwritten straight away. However some
  # zeros initialized in the var/covar matrices  will still remain and be used
  muqAlpha <- prior_mu_alpha
  muqBeta <- rep(0., MBeta)
  muqAlphaBeta <- c(muqAlpha, muqBeta)
  
  muqYMis <- rep(0, NMis)
  muqY <- vector()
  # sigma2qYMis <- rep(0, NMis)
  # sigma2qY <- vector()
  if(modelMiss){
    muqY <- c(yObs, muqYMis)
    #sigma2qY <- c(rep(0, NObs), sigma2qYMis)
  } else {
    muqY <- yObs
    #sigma2qY <- rep(0, NObs)
  }
  if(NMAR){
    TmuqY <- cbind(TT, muqY)
  } else {
    TmuqY <- TT
  }
  
  muqAlphaBetaYmis <- c(muqAlphaBeta, muqYMis)
  muqAlphaBetaY <- c(muqAlphaBeta, yObs, muqYMis)
  
  SigmaqAlphaBeta <- matrix(rep(0.,(MAlpha+MBeta)*(MAlpha+MBeta)),
                            ncol = (MAlpha+MBeta))
  SigmaqAlphaBetaYmis <- 
    matrix(rep(0., (MAlpha+MBeta+NMis) * (MAlpha+MBeta+NMis)),
           nrow = MAlpha + MBeta + NMis,
           ncol = MAlpha + MBeta + NMis)
  invSigmaqAlphaBetaYmis <- 
    matrix(rep(0., (MAlpha+MBeta+NMis) * (MAlpha+MBeta+NMis)),
           nrow = MAlpha + MBeta + NMis,
           ncol = MAlpha + MBeta + NMis)
  SigmaqAlphaBetaY <- 
    matrix(rep(0., (MAlpha+MBeta+N) * (MAlpha+MBeta+N)),
           nrow = MAlpha + MBeta + N,
           ncol = MAlpha + MBeta + N)
  SigmaqYmis <- 
    matrix(rep(0., NMis * NMis),
           nrow = NMis, ncol = NMis)
  
  SigmaqY <- 
    matrix(rep(0., N * N),
           nrow = N, ncol = N)
  
  
  # other quantities precomputation or memory allocation
  Aqeps <- prior_A_eps + 0.5*N
  Aqrdm <- prior_A_rdm + 0.5*MBetaPerRdmEffect
  Bqrdm <- rep(0., nRdm)
  prior_A_b <- nu/2; prior_B_b <- nu/2
  Aqb <- prior_A_b + 0.5
  Bqb <- rep(0., N)
  
  phimua <- rep(0., N)
  Phia <- rep(0., N)
  
  WmuqAlphaBeta <- W %*% muqAlphaBeta
  upRightBlockSigmaqAlphaBetaYmis <- matrix(nrow = MAlpha + MBeta, 
                                            ncol = NMis)
  tUpRightBlockSigmaqAlphaBetaYmis <- matrix(nrow = NMis,
                                             ncol = MAlpha + MBeta)
  downRightBlockSigmaqAlphaBetaYmis <- matrix(nrow = NMis, ncol = NMis)
  WSigmaqAlphaBetaWT <- W %*% SigmaqAlphaBeta %*% t(W)
  WTEqDW <- t(W) %*% eqD %*% W
  
  MatZeroGamma <- matrix()
  if(NMAR){
    MatZeroMGamma <- matrix(rep(0., (MGamma+1) * (MGamma+1)), ncol = MGamma+1)
  } else {
    MatZeroMGamma <- matrix(rep(0., MGamma * MGamma), ncol = MGamma)
  }
  
  muqYMinusWmuqAlphaBeta <- rep(0., N)
  muqaMinusTmuqYmuqGamma <- rep(0., N)
  USigmaqAlphaBetaYUT <- U %*% SigmaqAlphaBetaY %*% t(U)
  
  # prepare coefficients for contrasts
  nContrast <- InputData$nContrast[k]
  contrastCoefs <- list()
  if(nContrast > 0){
    for(i in 1:nContrast){
      contrastCoefs[[i]] <- 
        InputData$contrastCoefs[(conpos+(i-1)*MAlpha):(conpos+i*MAlpha-1)]
    }
  }
  
  
  # loop niter iterations
  # if(k == 10){
  #   bbb <- 1
  # }
  #browser()
  converged <- FALSE
  ELBOIndex <- 0
  for (l in 1:nIterMax){
    if(distribStyle == "NonFactorized"){
      if(modelMiss){
        # Step 1 : update Alpha, Beta, yMis
        if(NMAR){
          downRightBlockSigmaqAlphaBetaYmis <- 
            diag(Sigmaqgammagy[MGamma+1, MGamma+1] + 
                   muqgammagy[MGamma+1] * muqgammagy[MGamma+1] + 
                   eqinvsigma2eps * eqinvb[(NObs+1):N])
          
        } else {
          downRightBlockSigmaqAlphaBetaYmis <- 
            diag(eqinvsigma2eps * eqinvb[(NObs+1):N])
        }
        upRightBlockSigmaqAlphaBetaYmis <- 
          - eqinvsigma2eps * t(WMis) %*% eqD[(NObs+1):N, (NObs+1):N]
        
        tUpRightBlockSigmaqAlphaBetaYmis <- t(upRightBlockSigmaqAlphaBetaYmis)
        
        invSigmaqAlphaBetaYmis[1:(MAlpha+MBeta), 1:(MAlpha+MBeta)] <-
          eqinvsigma2eps * WTEqDW + muqInvSigmaAlphaBeta
        invSigmaqAlphaBetaYmis[1:(MAlpha+MBeta), 
                               (MAlpha+MBeta+1):(MAlpha+MBeta+NMis)] <-
          upRightBlockSigmaqAlphaBetaYmis
        invSigmaqAlphaBetaYmis[(MAlpha+MBeta+1):(MAlpha+MBeta+NMis), 
                               1:(MAlpha+MBeta)] <-
          tUpRightBlockSigmaqAlphaBetaYmis
        invSigmaqAlphaBetaYmis[(MAlpha+MBeta+1):(MAlpha+MBeta+NMis),
                               (MAlpha+MBeta+1):(MAlpha+MBeta+NMis)] <-
          downRightBlockSigmaqAlphaBetaYmis
        
        SigmaqAlphaBetaYmis <- solve(invSigmaqAlphaBetaYmis)
        SigmaqAlphaBeta <- SigmaqAlphaBetaYmis[1:(MAlpha+MBeta),
                                               1:(MAlpha+MBeta)]
        
        SigmaqYmis <- SigmaqAlphaBetaYmis[(MAlpha+MBeta+1):(MAlpha+MBeta+NMis), 
                                          (MAlpha+MBeta+1):(MAlpha+MBeta+NMis)]
        SigmaqY[(NObs+1):N,(NObs+1):N] <- 
          SigmaqYmis
        
        SigmaqAlphaBetaY[1:(MAlpha+MBeta), 1:(MAlpha+MBeta)] <-
          SigmaqAlphaBeta
        SigmaqAlphaBetaY[1:(MAlpha+MBeta), 
                         (MAlpha+MBeta+NObs+1):(MAlpha+MBeta+N)] <-
          SigmaqAlphaBetaYmis[1:(MAlpha+MBeta), 
                              (MAlpha+MBeta+1):(MAlpha+MBeta+NMis)]
        SigmaqAlphaBetaY[(MAlpha+MBeta+NObs+1):(MAlpha+MBeta+N), 
                         1:(MAlpha+MBeta)] <-
          SigmaqAlphaBetaYmis[(MAlpha+MBeta+1):(MAlpha+MBeta+NMis),
                              1:(MAlpha+MBeta)]
        SigmaqAlphaBetaY[(MAlpha+MBeta+NObs+1):(MAlpha+MBeta+N),
                         (MAlpha+MBeta+NObs+1):(MAlpha+MBeta+N)] <-
          SigmaqYmis
        
        muqAlphaBetaYmis <- eqinvsigma2eps * t(WObs) %*% eqD[1:NObs,1:NObs] %*% yObs + 
          muqInvSigmaAlphaBeta %*% prior_mu_alpha_beta
        if(NMAR){
          if(MGamma > 1){
            muqAlphaBetaYmis <- 
              c(muqAlphaBetaYmis,
                muqgammagy[MGamma+1] * 
                  ( muqa[(NObs+1):N] - TT[(NObs+1):N,] %*% muqgammagy[1:MGamma] + 
                      muqgammagy[MGamma+1] * yRef))
          } else {
            muqAlphaBetaYmis <- 
              c(muqAlphaBetaYmis,
                muqgammagy[MGamma+1] * 
                  ( muqa[(NObs+1):N] - TT[(NObs+1):N,] * muqgammagy[1:MGamma] + 
                      muqgammagy[MGamma+1] * yRef))
          }
          
        } else {
          muqAlphaBetaYmis <- 
            c(muqAlphaBetaYmis,
              rep(0., NMis))
        }
        muqAlphaBetaYmis <- SigmaqAlphaBetaYmis %*% muqAlphaBetaYmis
        
        muqAlphaBeta <- muqAlphaBetaYmis[1:(MAlpha+MBeta)] 
        muqYMis <- muqAlphaBetaYmis[(MAlpha+MBeta+1):(MAlpha+MBeta+NMis)]
        muqY[(NObs+1):N] <- muqYMis
        
        
        if(NMAR){
          TmuqY <- cbind(TT, muqY - yRef)
          MatZeroMGamma[MGamma+1,MGamma+1] <- tr(SigmaqYmis)
        } # otherwise no need to update these objects as they are fixed 
        
        USigmaqAlphaBetaYUT <- U %*% SigmaqAlphaBetaY %*% t(U)
        
        # Step 2 : update gamma
        Sigmaqgammagy <- solve(t(TmuqY) %*% TmuqY + MatZeroMGamma + 
                                 prior_inv_Sigma_gammagy)
        muqgammagy <- Sigmaqgammagy %*% 
          (t(TmuqY) %*% muqa + prior_inv_Sigma_gammagy %*% prior_mu_gammagy)
        
        # Step 3 : update a
        mua <- TmuqY %*% muqgammagy
        phimua <- dnorm(mua)
        Phia <- pnorm(TwoRMinus1*mua)
        muqa <- mua + TwoRMinus1*phimua/Phia
        
        
      } else { # we do not model missingness
        # Step 1 : update Alpha, Beta and NO step 2,3,4
        SigmaqAlphaBeta <- solve(eqinvsigma2eps * WTEqDW + muqInvSigmaAlphaBeta)
        
        muqAlphaBeta <- SigmaqAlphaBeta %*% 
          (eqinvsigma2eps * t(WObs) %*% eqD[1:NObs,1:NObs] %*% yObs + 
             muqInvSigmaAlphaBeta %*% prior_mu_alpha_beta)
        
        WSigmaqAlphaBetaWT <- W %*% SigmaqAlphaBeta %*% t(W)
      }
      
      #WmuqAlphaBeta <- W %*% muqAlphaBeta
      
      muqYMinusWmuqAlphaBeta <- muqY - W %*% muqAlphaBeta
      
      # Step 4 : update sigma2_epsilon
      if(modelMiss){
        Bqeps <- as.vector(
          prior_B_eps + 0.5 *
            (t(muqYMinusWmuqAlphaBeta) %*% eqD %*% muqYMinusWmuqAlphaBeta +
               tr(eqD %*% USigmaqAlphaBetaYUT)))
      } else {
        Bqeps <- as.vector(
          prior_B_eps + 0.5 *
            (t(muqYMinusWmuqAlphaBeta) %*% eqD %*% muqYMinusWmuqAlphaBeta +
               tr(eqD %*% WSigmaqAlphaBetaWT)))
      }
      # cat("SSE = ", t(muqYMinusWmuqAlphaBeta) %*% eqD %*% muqYMinusWmuqAlphaBeta, "\n")
      # cat("2*Bqeps = ", 2*(Bqeps-prior_B_eps), "\n")
      # cat("\n")
      # 
      eqinvsigma2eps <- Aqeps / Bqeps
      
      # Step 6 : update sigma2_rdm
      if(nRdm>0){
        for(r in 1:nRdm){
          start_index_r <- MAlpha + 1
          if(r>1) start_index_r <- start_index_r + 
              sum(MBetaPerRdmEffect[1:(r-1)])
          index_r <- seq(start_index_r, by=1, length.out = MBetaPerRdmEffect[r])
          
          Bqrdm[r] <- as.vector(prior_B_rdm[r] + 0.5 * 
                                  (t(muqAlphaBeta[index_r]) %*% muqAlphaBeta[index_r] +
                                     tr(SigmaqAlphaBeta[index_r, index_r])))
        }
        eqinvsigma2rdm <- Aqrdm / Bqrdm
        muqInvSigmaAlphaBeta <- blockDiag(prior_inv_sigma2_alpha_mat, 
                                          diag(rep(eqinvsigma2rdm, 
                                                   MBetaPerRdmEffect)))
      }
      
      # Step 7 : update bi's
      if(t_likelihood){
        if(modelMiss){
          Bqb <- as.vector(
            prior_B_b + 0.5 * eqinvsigma2eps *
              (muqYMinusWmuqAlphaBeta * muqYMinusWmuqAlphaBeta + 
                 diag(USigmaqAlphaBetaYUT))
          )
        } else {
          Bqb <- as.vector(
            prior_B_b + 0.5 * eqinvsigma2eps *
              (muqYMinusWmuqAlphaBeta * muqYMinusWmuqAlphaBeta + 
                 diag(WSigmaqAlphaBetaWT))
          )
        }
        eqinvb <- Aqb/Bqb
        eqD <- diag(eqinvb)
        WTEqDW <- t(W) %*% eqD %*% W
        
      }
      # end of update parameters part for non factorized distributions style
    } else {
      if(modelMiss){
        # Step 1 : update yMis
        if(NMAR){
          
          sigma2qYMis <- 1/
            (Sigmaqgammagy[MGamma+1, MGamma+1] + 
               muqgammagy[MGamma+1] * muqgammagy[MGamma+1] + 
               eqinvsigma2eps * eqinvb[(NObs+1):N])
          if(MGamma > 1){
            if(MAlpha+MBeta > 1){
              muqYMis <- sigma2qYMis * 
                (muqgammagy[MGamma+1] *
                   ( muqa[(NObs+1):N] - TT[(NObs+1):N,] %*% muqgammagy[1:MGamma] + 
                       muqgammagy[MGamma+1] * yRef) + 
                   eqinvsigma2eps * eqinvb[(NObs+1):N] * 
                   (W[(NObs+1):N,] %*% muqAlphaBeta))
            } else {
              muqYMis <- sigma2qYMis * 
                (muqgammagy[MGamma+1] *
                   ( muqa[(NObs+1):N] - TT[(NObs+1):N,] %*% muqgammagy[1:MGamma] + 
                       muqgammagy[MGamma+1] * yRef) + 
                   eqinvsigma2eps * eqinvb[(NObs+1):N] * 
                   (W[(NObs+1):N,] * c(muqAlphaBeta)))
            }
            
          } else {
            if(MAlpha+MBeta > 1){
              muqYMis <- sigma2qYMis * 
                (muqgammagy[MGamma+1] * 
                   ( muqa[(NObs+1):N] - TT[(NObs+1):N,] * muqgammagy[1:MGamma] + 
                       muqgammagy[MGamma+1] * yRef) + 
                   eqinvsigma2eps * eqinvb[(NObs+1):N] * 
                   (W[(NObs+1):N,] %*% muqAlphaBeta))
            } else 
            {
              muqYMis <- sigma2qYMis * 
                (muqgammagy[MGamma+1] * 
                   ( muqa[(NObs+1):N] - TT[(NObs+1):N,] * muqgammagy[1:MGamma] + 
                       muqgammagy[MGamma+1] * yRef) + 
                   eqinvsigma2eps * eqinvb[(NObs+1):N] * 
                   (W[(NObs+1):N,] * c(muqAlphaBeta)))
            }
            
          }
          
        } else {
          if(MAlpha + MBeta > 1){
            sigma2qYMis <- 1/(eqinvsigma2eps * eqinvb[(NObs+1):N])
            muqYMis <- sigma2qYMis * eqinvsigma2eps * eqinvb[(NObs+1):N] * 
              (W[(NObs+1):N,] %*% muqAlphaBeta)
          } else {
            sigma2qYMis <- 1/(eqinvsigma2eps * eqinvb[(NObs+1):N])
            muqYMis <- sigma2qYMis * eqinvsigma2eps * eqinvb[(NObs+1):N] * 
              (W[(NObs+1):N,] * c(muqAlphaBeta))
          }
          
        }
        
        muqY <- c(yObs, muqYMis)
        #sigma2qY[(NObs+1):N] <- sigma2qYMis
        SigmaqY <- diag(c(rep(0.,NObs), sigma2qYMis))
        
        if(NMAR){
          TmuqY <- cbind(TT, muqY - yRef)
          MatZeroMGamma[MGamma+1,MGamma+1] <- sum(sigma2qYMis)
        } # otherwise no need to update these objects
        
        # Step 2 : update a
        mua <- TmuqY %*% muqgammagy
        phimua <- dnorm(mua)
        Phia <- pnorm(TwoRMinus1*mua)
        muqa <- mua + TwoRMinus1*phimua/Phia
        
        # Step 3 : update gamma
        Sigmaqgammagy <- solve(t(TmuqY) %*% TmuqY + MatZeroMGamma + 
                                 prior_inv_Sigma_gammagy)
        muqgammagy <- Sigmaqgammagy %*% 
          (t(TmuqY) %*% muqa + prior_inv_Sigma_gammagy %*% prior_mu_gammagy)
      }
      
      # Step 4 : update Alpha and Beta
      SigmaqAlphaBeta <-solve(eqinvsigma2eps * WTEqDW + 
                                muqInvSigmaAlphaBeta)
      muqAlphaBeta <- SigmaqAlphaBeta %*% 
        (eqinvsigma2eps * t(W) %*% eqD %*% muqY +
           muqInvSigmaAlphaBeta %*% prior_mu_alpha_beta)
      
      WmuqAlphaBeta <- W %*% muqAlphaBeta
      WSigmaqAlphaBetaWT <- W %*% SigmaqAlphaBeta %*% t(W)
      
      # Step 5 : update sigma2_epsilon
      Bqeps <- as.vector(prior_B_eps + 0.5 *
                           (t(muqY) %*% eqD %*% muqY + tr(eqD %*% SigmaqY) - 
                              2 * t(muqAlphaBeta) %*% t(W) %*% eqD %*% muqY + 
                              t(muqAlphaBeta) %*% WTEqDW %*% muqAlphaBeta +
                              tr(eqD %*% W %*% SigmaqAlphaBeta %*% t(W))))
      # cat("SSE = ", t(muqY) %*% eqD %*% muqY - 
      #       2 * t(muqAlphaBeta) %*% t(W) %*% eqD %*% muqY + 
      #       t(muqAlphaBeta) %*% WTEqDW %*% muqAlphaBeta, "\n")
      # cat("2*Bqeps = ", 2*(Bqeps-prior_B_eps), "\n")
      # cat("\n")
      
      eqinvsigma2eps <- Aqeps / Bqeps
      
      # Step 6 : update sigma2_rdm
      if(nRdm>0){
        for(r in 1:nRdm){
          start_index_r <- MAlpha + 1
          if(r>1) start_index_r <- start_index_r + 
              sum(MBetaPerRdmEffect[1:(r-1)])
          index_r <- seq(start_index_r, by=1, length.out = MBetaPerRdmEffect[r])
          
          Bqrdm[r] <- as.vector(prior_B_rdm[r] + 0.5 * 
                                  (t(muqAlphaBeta[index_r]) %*% muqAlphaBeta[index_r] +
                                     tr(SigmaqAlphaBeta[index_r, index_r])))
        }
        eqinvsigma2rdm <- Aqrdm / Bqrdm
        muqInvSigmaAlphaBeta <- blockDiag(prior_inv_sigma2_alpha_mat, 
                                          diag(rep(eqinvsigma2rdm, 
                                                   MBetaPerRdmEffect)))
      }
      
      # Step 7 : update bi's
      if(t_likelihood){
        Bqb <- mapply(FUN = function(muqY, SigmaqY, WmuqAlphaBeta, WSigmaqAlphaBetaWT){
          x <- muqY * muqY + SigmaqY - 2 * WmuqAlphaBeta * muqY +
            WmuqAlphaBeta * WmuqAlphaBeta + WSigmaqAlphaBetaWT
          x},
          muqY = muqY,
          SigmaqY = diag(SigmaqY),
          WmuqAlphaBeta = WmuqAlphaBeta,
          WSigmaqAlphaBetaWT = diag(WSigmaqAlphaBetaWT))
        Bqb <- as.vector(prior_B_b + 0.5 * eqinvsigma2eps *  Bqb)
        eqinvb <- Aqb/Bqb
        eqD <- diag(eqinvb)
        WTEqDW <- t(W) %*% eqD %*% W
      }
      
      # end of update parameters part for factorized distributions style
    }
    
    if(l == 1 || (l %% ELBO_each) == 0 || l == nIterMax){
      
      ## START COMPUTING ELBO
      ELBOIndex <- ELBOIndex + 1
      ELBOIter[ELBOIndex] <- l
      
      if(distribStyle == "NonFactorized"){
        if(modelMiss){
          ## ELBO STEP 1
          ELBODetails[ELBOIndex,1] <- 
            sum(log(pnorm(mua[1:NObs]))) +
            sum(log(rep(1,NMis)-pnorm(mua[(NObs+1):(NObs+NMis)]))) - 0.5 *
            tr(t(TmuqY) %*% TmuqY %*% Sigmaqgammagy)
          
          if(NMAR){
            ELBODetails[ELBOIndex,1] <- ELBODetails[ELBOIndex,1] - 0.5 *
              tr(SigmaqY) * ( muqgammagy[MGamma+1] * muqgammagy[MGamma+1] + 
                                Sigmaqgammagy[MGamma+1,MGamma+1] )
          } 
          
          ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] + 
            ELBODetails[ELBOIndex,1]
          
          # ELBO STEP 2
          ELBODetails[ELBOIndex,2] <- 
            0.5 * determinant(Sigmaqgammagy, logarithm = TRUE)$modulus -
            0.5 * (t(muqgammagy - prior_mu_gammagy ) %*% 
                     prior_inv_Sigma_gammagy %*% 
                     (muqgammagy - prior_mu_gammagy)) -
            0.5 * tr(prior_inv_Sigma_gammagy %*% Sigmaqgammagy)
          ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] + ELBODetails[ELBOIndex,2]
        }
        
        # ELBO STEP 3
        if(modelMiss){
          ELBODetails[ELBOIndex,3] <- - 0.5 * eqinvsigma2eps *
            (t(muqYMinusWmuqAlphaBeta) %*% eqD %*% muqYMinusWmuqAlphaBeta +
               tr(eqD %*% U %*% SigmaqAlphaBetaY %*% t(U))) - 0.5 *
            (t(muqAlphaBeta - prior_mu_alpha_beta) %*% muqInvSigmaAlphaBeta %*%
               (muqAlphaBeta - prior_mu_alpha_beta) +
               tr(muqInvSigmaAlphaBeta %*% SigmaqAlphaBeta)) +
            0.5 * determinant(SigmaqAlphaBetaYmis, logarithm = TRUE)$modulus
        } else {
          ELBODetails[ELBOIndex,3] <- - 0.5 * eqinvsigma2eps *
            (t(muqYMinusWmuqAlphaBeta) %*% eqD %*% muqYMinusWmuqAlphaBeta +
               tr(eqD %*% W %*% SigmaqAlphaBeta %*% t(W))) - 0.5 *
            (t(muqAlphaBeta - prior_mu_alpha_beta) %*% muqInvSigmaAlphaBeta %*%
               (muqAlphaBeta - prior_mu_alpha_beta) +
               tr(muqInvSigmaAlphaBeta %*% SigmaqAlphaBeta)) +
            0.5 * determinant(SigmaqAlphaBeta, logarithm = TRUE)$modulus
        }
        ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] + ELBODetails[ELBOIndex,3]
        
        # ELBO STEP 4
        if(nRdm > 0){
          ELBODetails[ELBOIndex,4] <- 
            sum(eqinvsigma2rdm * (Bqrdm - prior_B_rdm)) -
            sum(Aqrdm * log(Bqrdm))
          ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] + ELBODetails[ELBOIndex,4]
        }
        
        # ELBO STEP 5
        ELBODetails[ELBOIndex,5] <- 
          eqinvsigma2eps * (Bqeps - prior_B_eps) -
          Aqeps * log(Bqeps)
        ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] + ELBODetails[ELBOIndex,5]
        
        # ELBO STEP 6
        if(t_likelihood > 0){
          ELBODetails[ELBOIndex,6] <- 
            sum(eqinvb * (Bqb - prior_B_b)) -
            sum(Aqb * log(Bqb))
          ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] + ELBODetails[ELBOIndex,6]
        }
        
        ## END ELBO COMPUTATION
        
        # end of ELBO update part for non factorized distributions style
      } else {
        if(modelMiss){
          ## ELBO STEP 1
          ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] +
            sum(log(pnorm(mua[1:NObs]))) +
            sum(log(rep(1,NMis)-pnorm(mua[(NObs+1):(NObs+NMis)]))) - 0.5 *
            tr(t(TmuqY) %*% TmuqY %*% Sigmaqgammagy)
          
          if(NMAR){
            ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] - 0.5 *
              tr(SigmaqY) * ( muqgammagy[MGamma+1] * muqgammagy[MGamma+1] + 
                                Sigmaqgammagy[MGamma+1,MGamma+1])
          } 
          
          # ELBO STEP 2
          ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] +
            0.5 * determinant(Sigmaqgammagy, logarithm = TRUE)$modulus - 
            0.5 * (t(muqgammagy - prior_mu_gammagy ) %*% 
                     prior_inv_Sigma_gammagy %*%
                     (muqgammagy - prior_mu_gammagy)) -
            0.5 * tr(prior_inv_Sigma_gammagy %*% Sigmaqgammagy)
        }
        
        # ELBO STEP 3
        ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] - 0.5 * eqinvsigma2eps * 
          (t(muqY - WmuqAlphaBeta) %*% eqD %*% (muqY - WmuqAlphaBeta) +
             tr(eqD %*% (SigmaqY + W %*% SigmaqAlphaBeta %*% t(W))))
        if(modelMiss){
          ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] +
            0.5 * determinant(SigmaqY[(NObs+1):N,(NObs+1):N], logarithm = TRUE)$modulus
        }
        
        
        # ELBO STEP 4
        ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] - 0.5 *
          (t(muqAlphaBeta - prior_mu_alpha_beta) %*% muqInvSigmaAlphaBeta %*%
             (muqAlphaBeta - prior_mu_alpha_beta) + 
             tr(muqInvSigmaAlphaBeta %*% SigmaqAlphaBeta)) +
          0.5 * determinant(SigmaqAlphaBeta, logarithm = TRUE)$modulus
        
        #ELBO STEP 5  
        if(nRdm > 0){
          ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] + 
            sum(eqinvsigma2rdm * (Bqrdm - prior_B_rdm)) - 
            sum(Aqrdm * log(Bqrdm))
        }  
        
        #ELBO STEP 6
        ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] +
          eqinvsigma2eps * (Bqeps - prior_B_eps) - 
          Aqeps * log(Bqeps)
        
        #ELBO STEP 7
        if(t_likelihood > 0){
          ELBOValue[ELBOIndex] <- ELBOValue[ELBOIndex] + 
            sum(eqinvb * (Bqb - prior_B_b)) - 
            sum(Aqb * log(Bqb))
        }  
        # end of ELBO update part for factorized distributions style
      }
      
    }
    if(ELBOIndex>1){
      if(abs(ELBOValue[ELBOIndex] - ELBOValue[ELBOIndex-1]) < 
         relTol * abs(ELBOValue[ELBOIndex])){
        printIteration(l, MAlpha, MBeta, MGamma, NObs, N, modelMiss, NMAR,
                       nRdm, muqAlphaBeta, SigmaqAlphaBeta, eqinvsigma2eps, 
                       eqinvsigma2rdm, muqgammagy, Sigmaqgammagy, muqY, 
                       SigmaqY, mua, muqa,
                       ELBOValue, ELBOIndex)
        converged <- TRUE
        break
      }
    }
    if(l %% refresh == 0 || l == nIterMax){
      printIteration(l, MAlpha, MBeta, MGamma, NObs, N, modelMiss, NMAR,
                     nRdm, muqAlphaBeta, SigmaqAlphaBeta, eqinvsigma2eps, 
                     eqinvsigma2rdm, muqgammagy, Sigmaqgammagy, muqY, 
                     SigmaqY, mua, muqa,
                     ELBOValue, ELBOIndex)
    }
    
  }# end loop on l (iteration)
  
  end_time <- Sys.time()
  
  # store general fit results
  res$fitResults <- list( nIter = l, 
                          converged = converged,
                          nELBO = ELBOIndex,
                          ELBOIter = ELBOIter[1:ELBOIndex],
                          ELBOValue = ELBOValue[1:ELBOIndex],
                          ELBODetails = ELBODetails[1:ELBOIndex,],
                          elapsed = end_time - start_time)
  
  
  # update marginal distributions with current protein results
  
  # normal distributions for alpha
  newParamDistrib <- list()
  for(i in 1:MAlpha){
    mui <- muqAlphaBeta[i]
    if(MAlpha + MBeta > 1){
      sdi <- sqrt(SigmaqAlphaBeta[i,i])
    } else {
      sdi <- sqrt(SigmaqAlphaBeta[1])
    }
    newParamDistrib[[i]] <- list(main = TRUE,
                                 mean = mui,
                                 sd = sdi,
                                 type = "norm",
                                 args = list(mean = mui,
                                             sd = sdi))
  }
  #names(newParamDistrib) <- res$coefficients$alphaNames[apos:(apos+MAlpha-1)]
  names(newParamDistrib) <- paste0("alpha[", 
                                   seq(from = apos, to = apos+MAlpha-1),
                                   "]")
  res$paramDistr <- newParamDistrib
  
  # normal distributions for contrasts on alpha's
  if(nContrast > 0){
    newParamDistrib <- list()
    for(i in 1:nContrast){
      mui <- sum(muqAlphaBeta[1:MAlpha]*contrastCoefs[[i]])
      if(MAlpha + MBeta > 1){
        vari <-  t(contrastCoefs[[i]]) %*% SigmaqAlphaBeta[1:MAlpha,1:MAlpha] %*% 
          contrastCoefs[[i]]
      } else {
        vari <- SigmaqAlphaBeta[1] * contrastCoefs[[i]][1] * contrastCoefs[[i]][1]
      }
      sdi <- sqrt(vari)
      newParamDistrib[[i]] <- list(main = TRUE,
                                   mean = mui,
                                   sd = sdi,
                                   type = "norm",
                                   args = list(mean = mui,
                                               sd = sdi))
    }
    
    if(k == 1){
      start <- 1
    } else {
      start <- sum(InputData$nContrast[1:(k-1)]) + 1
    }
    names(newParamDistrib) <- 
      paste0("contrast[", 
             seq(from = start, to = start+nContrast-1),
             "]")
    
    res$paramDistr <- c(res$paramDistr, newParamDistrib)
  }
  
  
  # normal distributions for beta
  if(nRdm > 0){
    newParamDistrib <- list()
    for(i in 1:MBeta){
      mui <- muqAlphaBeta[MAlpha+i]
      if(MAlpha + MBeta > 1){
        sdi <- sqrt(SigmaqAlphaBeta[MAlpha+i, MAlpha+i])
      } else {
        sdi <- sqrt(SigmaqAlphaBeta[1])
      }
      newParamDistrib[[i]] <- list(main = TRUE,
                                   mean = mui,
                                   sd = sdi,
                                   type = "norm",
                                   args = list(mean = mui,
                                               sd = sdi))
    }
    #names(newParamDistrib) <- res$coefficients$betaNames[bpos:(bpos+MBeta-1)]
    names(newParamDistrib) <- paste0("beta[", 
                                     seq(from = bpos, to = bpos+MBeta-1),
                                     "]")
    res$paramDistr <- c(res$paramDistr, newParamDistrib)
  }
  
  # normal distributions for gamma and gammaY
  if(modelMiss){
    newParamDistrib <- list()
    for(i in 1:MGamma){
      mui <- muqgammagy[i]
      if(MGamma > 1){
        sdi <- sqrt(Sigmaqgammagy[i, i])
      } else {
        sdi <- sqrt(Sigmaqgammagy[1])
      }
      newParamDistrib[[i]] <- list(main = TRUE,
                                   mean = mui,
                                   sd = sdi,
                                   type = "norm",
                                   args = list(mean = mui,
                                               sd = sdi))
    }
    #distrNames <- res$coefficients$gammaNames[cpos:(cpos+MGamma-1)]
    distrNames <- paste0("gamma[", 
                         seq(from = cpos, to = cpos+MGamma-1),
                         "]")
    
    if(NMAR){
      gammay_name <- paste0("gammY[", k, "]")
      mui <- muqgammagy[MGamma+1]
      sdi <- sqrt(Sigmaqgammagy[MGamma+1,MGamma+1])
      newParamDistrib[[MGamma+1]] <- list(main = TRUE,
                                          mean = mui,
                                          sd = sdi,
                                          type = "norm",
                                          args = list(mean = mui,
                                                      sd = sdi))
      distrNames <- c(distrNames, gammay_name)
    }
    names(newParamDistrib) <- distrNames
    res$paramDistr <- c(res$paramDistr, newParamDistrib)
    
  }
  
  # normal distributions for Ymis
  #if(!filterOutput && modelMiss){
  if(modelMiss){
    newParamDistrib <- list()
    for(i in 1:NMis){
      mui <- muqY[NObs+i]
      sdi <- sqrt(SigmaqY[NObs+i,NObs+i])
      
      newParamDistrib[[i]] <- list(main = FALSE,
                                   mean = mui,
                                   sd = sdi,
                                   type = "norm",
                                   args = list(mean = mui,
                                               sd = sdi))
    }
    #names(newParamDistrib) <- res$coefficients$yMisNames[ymispos:(ymispos+NMis-1)]
    names(newParamDistrib) <- paste0("yMis[", 
                                     seq(from = ymispos, to = ymispos+NMis-1),
                                     "]")
    
    res$paramDistr <- c(res$paramDistr, newParamDistrib)
  }
  
  
  # gamma distributions for tau(eps)
  tau_name <- paste0("tauEps[", k, "]")
  newParamDistrib <- list()
  newParamDistrib[[1]] <- list(main = TRUE,
                               mean = Aqeps/Bqeps,
                               sd = sqrt(Aqeps/Bqeps/Bqeps),
                               type = "gamma",
                               args = list(shape = Aqeps,
                                           rate = Bqeps))
  names(newParamDistrib) <- c(tau_name)
  res$paramDistr <- c(res$paramDistr, newParamDistrib)
  
  # inverse gamma distributions for sigma2(eps)
  sigma2_name <- paste0("sigma2Eps[", k, "]")
  mui <- Bqeps/(Aqeps-1)
  sdi <- sqrt(Bqeps*Bqeps/((Aqeps-1)^2*(Aqeps-2)))
  newParamDistrib <- list()
  newParamDistrib[[1]] <- list(main = TRUE,
                               mean = mui,
                               sd = sdi,
                               type = "invgamma",
                               args = list(alpha = Aqeps,
                                           beta = Bqeps))
  names(newParamDistrib) <- c(sigma2_name)
  res$paramDistr <- c(res$paramDistr, newParamDistrib)
  
  # gamma distributions for tau_rdm and inverse gamma distributions for 
  # sigma2_rdm
  if(nRdm > 0){
    newParamDistrib <- list()
    for(i in 1:nRdm){
      newParamDistrib[[i]] <- list(main = TRUE,
                                   mean = Aqrdm[i]/Bqrdm[i],
                                   sd = sqrt(Aqrdm[i]/Bqrdm[i]/Bqrdm[i]),
                                   type = "gamma",
                                   args = list(shape = Aqrdm[i],
                                               rate = Bqrdm[i]))
    }
    
    names(newParamDistrib) <- paste0("tauRdm[",
                                     seq(from = spos, to = spos+nRdm-1),
                                     "]")
    
    res$paramDistr <- c(res$paramDistr, newParamDistrib)
    
    newParamDistrib <- list()
    for(i in 1:nRdm){
      alpha <- Aqrdm[i]
      beta <- Bqrdm[i]
      mui <- beta/(alpha-1)
      sdi <- sqrt(beta*beta/((alpha-1)^2*(alpha-2)))
      newParamDistrib[[i]] <- list(main = TRUE,
                                   mean = mui,
                                   sd = sdi,
                                   type = "invgamma",
                                   args = list(alpha = alpha,
                                               beta = beta))
    }
    
    names(newParamDistrib) <- paste0("sigma2Rdm[",
                                     seq(from = spos, to = spos+nRdm-1),
                                     "]")
    
    res$paramDistr <- c(res$paramDistr, newParamDistrib)
  }
  
  # gamma distribution for inv(b)
  #if(!filterOutput && t_likelihood > 0){
  if(t_likelihood > 0){
    # inv(b) parameters for observed data
    newParamDistrib <- list()
    for(i in 1:NObs){
      newParamDistrib[[i]] <- list(main = FALSE,
                                   mean = Aqb/Bqb[i],
                                   sd = sqrt(Aqb/Bqb[i]/Bqb[i]),
                                   type = "gamma",
                                   args = list(shape = Aqb,
                                               rate = Bqb[i]))
    }
    names(newParamDistrib) <- paste0("invBObs[", 
                                     seq(from = yobspos, to = yobspos+NObs-1),
                                     "]")
    
    res$paramDistr <- c(res$paramDistr, newParamDistrib)
    
    # inv(b) parameters for missing data => only if we don't filter output
    if(modelMiss){
      newParamDistrib <- list()
      for(i in 1:NMis){
        newParamDistrib[[i]] <- list(main = FALSE,
                                     mean = Aqb/Bqb[NObs+i],
                                     sd = sqrt(Aqb/Bqb[NObs+i]/Bqb[NObs+i]),
                                     type = "gamma",
                                     args = list(shape = Aqb,
                                                 rate = Bqb[NObs+i]))
      }
      names(newParamDistrib) <- paste0("invBMis[", 
                                       seq(from = ymispos, to = ymispos+NMis-1),
                                       "]")
      
      res$paramDistr <- c(res$paramDistr, newParamDistrib)
    }
    
  }
  
  # for 'a' parameter : truncated normals
  #if(!filterOutput && modelMiss){
  if(modelMiss){
    # for observed responses
    newParamDistrib <- list()
    for(i in 1:NObs){
      ratioi <- dnorm(mua[i])/pnorm(mua[i])
      sdi <- sqrt(1 - mua[i]*ratioi - ratioi^2)
      newParamDistrib[[i]] <- list(main = FALSE,
                                   mean = muqa[i],
                                   sd = sdi,
                                   type = "normTrunc",
                                   args = list(mean = mua[i],
                                               sd = 1,
                                               min = 0.))
    }
    names(newParamDistrib) <- paste0("aObs[", 
                                     seq(from = yobspos, to = yobspos+NObs-1),
                                     "]")
    
    res$paramDistr <- c(res$paramDistr, newParamDistrib)
    
    # for missing responses
    newParamDistrib <- list()
    for(i in 1:NMis){
      ratioi <- dnorm(mua[NObs+i])/(1-pnorm(mua[NObs+i]))
      sdi <- sqrt(1 + mua[NObs+i]*ratioi - ratioi^2)
      newParamDistrib[[i]] <- list(main = FALSE,
                                   mean = muqa[NObs+i],
                                   sd = sdi,
                                   type = "normTrunc",
                                   args = list(mean = mua[NObs+i],
                                               sd = 1,
                                               max = 0.))
    }
    names(newParamDistrib) <- paste0("aMis[", 
                                     seq(from = ymispos, to = ymispos+NMis-1),
                                     "]")
    
    res$paramDistr <- c(res$paramDistr, newParamDistrib)
  }
    
  
  return(res)
  
}

##' @import BiocParallel
##'
##' @export
runMFVBByProt <- function(x, InputData)
{
  cores <- algoArgs(x)$cores
  if(is.null(cores)){
    cores = 1
  }
  
  nProt <- InputData$K
  protIndexes <- c(1:nProt)
  
  if(cores > 1 && nProt > 1){
    # take registered bpparam if suitable for parallel computation
    # if not, create SnowParam by default as suitable for any platform
    used_bp <- bpparam() 
    if(is.null(used_bp) || inherits(used_bp, "SerialParam")){
      used_bp <- SnowParam(workers = cores,
                           tasks = nProt,
                           progressbar = TRUE,
                           type = "SOCK",
                           #RNGseed = my.seed,
                           timeout = 30L * 24L * 60L * 60L,
                           exportglobals = TRUE,
                           log = FALSE,
                           logdir = NA_character_,
                           jobname = "MFVB")
    }
    
    res <- bpmapply(FUN = runMFVBOneProt, 
                    protIndex = protIndexes,
                    MoreArgs = list(x = x, 
                                    InputData = InputData),
                    SIMPLIFY = FALSE,
                    BPPARAM = used_bp)
    
  } else {
    res <- mapply(FUN = runMFVBOneProt, 
                  protIndex = protIndexes,
                  MoreArgs = list(x = x, InputData = InputData),
                  SIMPLIFY = FALSE)
  }
  
  # store everything back into object x
  x@fitResults <- list()
  x@paramDistr <- list()
  for(i in 1:nProt){
    x@fitResults[[i]] <- res[[i]]$fitResults
    x@paramDistr <- c(x@paramDistr, res[[i]]$paramDistr)
  }
  x
}

##' @import BiocParallel
##'
##' @export
runMFVBAllProtsAtOnce <- function(x, InputData)
{
  cores <- algoArgs(x)$cores
  if(is.null(cores)){
    cores = 1
  }
  
  nProt <- InputData$K
  protIndexes <- c(1:nProt)
  used_bp = NULL;
  
  if(cores > 1 && nProt > 1){
    # take registered bpparam if suitable for parallel computation
    # if not, create SnowParam by default as suitable for any platform
    used_bp <- bpparam() 
    if(is.null(used_bp) || inherits(used_bp, "SerialParam")){
      used_bp <- SnowParam(workers = cores,
                           tasks = 0L,
                           progressbar = TRUE,
                           type = "SOCK",
                           #RNGseed = my.seed,
                           timeout = 30L * 24L * 60L * 60L,
                           exportglobals = TRUE,
                           log = FALSE,
                           logdir = NA_character_,
                           jobname = "MFVB")
    }
  }
  
  x <- runMFVBAllProts(x = x, 
                       InputData = InputData, 
                       bp = used_bp)
    
  x
}

##' Take an initialised `ProtModelFitMFVB` object, runs MFVB approximation
##' algorithm and returns an updated `ProtModelFitMFVB` object containing the 
##' estimated parameters distributions
##' 
##' @title Run a Bayes model using variational approximation
##' 
##' @param x A `ProtModelFitMFVB` object
##'
##' @return A updated `ProtModelFitMFVB` object.
##' 
##' @importFrom stats dnorm qnorm pnorm rnorm dgamma qgamma rgamma density quantile sd
##' @importFrom EnvStats dnormTrunc qnormTrunc rnormTrunc
##' @importFrom extraDistr dinvgamma qinvgamma rinvgamma
##' @import BiocParallel
##'
##' @export
##' 
##' @author Philippe Hauchamps  
runProtModelEngine.ProtModelFitMFVB <- function(x, 
                                                InputData,
                                                initialValues)
{
  stopifnot(inherits(x, "ProtModelFitMFVB"))
  
  if( length(InputData) == 0) 
    stop("InputData list is empty")
  
  if(is.null(ModelSpecPars(x))){
    stop("Model Specs not initialized!")
  }
  
  if(is.null(algoArgs(x))){
    stop("Algorithm arguments not initialized!")
  }
  
  hierarchicalModel <- sigmaHierarch(ModelSpecPars(x)) || 
    sigmaRdmHierarch(ModelSpecPars(x)) ||
    gammaHierarch(ModelSpecPars(x))
  
  byProt <- algoArgs(x)$byProt
  if(is.null(byProt)){
    if(hierarchicalModel){
      byProt = FALSE
    } else {
      byProt = TRUE
    }
    
  }
  
  if(byProt && hierarchicalModel){
    msg <- "detected inconsistency in AlgoArgs : byProt cannot be [TRUE] when "
    msg <- paste0(msg, "model is hierarchical")
    stop(msg)
  }
  
  if(byProt){
    x <- runMFVBByProt(x = x, InputData = InputData)
  } else {
    x <- runMFVBAllProtsAtOnce(x = x, InputData = InputData)
  }
  
  x
  
}
  
  
setMethod("runProtModelEngine", "ProtModelFitMFVB",
          runProtModelEngine.ProtModelFitMFVB)

