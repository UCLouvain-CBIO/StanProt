// saved as FullHierarchicalLaplaceApproxi2.stan

data {
  int<lower=1> K;                  // number of proteins
  int<lower=0> nRdm;               // number of random effects
  int<lower=1> NObs[K];            // number of observations per protein
  int<lower=1> NMis[K];            // number of missing data per protein
  int<lower=1> MAlpha[K];          // number of alpha's (fixed effects coefs) per protein
  int<lower=0> MBeta[K];           // number of beta's (random effects coefs) per protein
  int<lower=1> MBetaPerRdmEffect[nRdm*K]; // same as before, per random effect
  int<lower=1> MGamma[K];          // number of gamma coeffs (missingness model)
  int<lower=1> NObsSum;            // total number of observations (=sum(NObs))
  int<lower=1> NMisSum;            // total number of missing data (=sum(NMis))
  int<lower=1> XObsLen;            // dimension of the flat vector containing all model matrices 
                                   // for the fixed effects of the observed responses
  int<lower=1> XMisLen;            // dimension of the flat vector containing all model matrices 
                                   // for the fixed effects of the missing responses
  int<lower=0> ZObsLen;            // dimension of the flat vector containing all model matrices 
                                   // for the random effects of the observed responses
  int<lower=0> ZMisLen;            // dimension of the flat vector containing all model matrices 
                                   // for the random effects of the missing responses
  int<lower=1> TObsLen;            // dimension of the flat vector containing all model matrices 
                                   // for the fixed effects of the missing indicators equal 0
  int<lower=1> TMisLen;            // dimension of the flat vector containing all model matrices 
                                   // for the fixed effects of the missing indicators equal 1
  
  row_vector[XObsLen] XObs;        // concatenation of model matrices in a flat row vector format, 
                                   // for the fixed effects of the observed responses
  row_vector[XMisLen] XMis;        // concatenation of model matrices in a flat row vector format, 
                                   // for the fixed effects of the missing responses
  vector[NObsSum] yObs;            // observed responses (log-2 intensities of peptides)
  row_vector[ZObsLen] ZObs;        // concatenation of model matrices in a flat row vector format, 
                                   // for the random effects of the observed responses
  row_vector[ZMisLen] ZMis;        // concatenation of model matrices in a flat row vector format, 
                                   // for the random effects of the missing responses
  row_vector[TObsLen] TObs;        // concatenation of model matrices in a flat row vector format, 
                                   // for the fixed effects of the missing indicators equal 0
  row_vector[TMisLen] TMis;        // concatenation of model matrices in a flat row vector format, 
                                   // for the fixed effects of the missing indicators equal 1
  
  // general model settings
  int<lower=0, upper=1> t_likelihood;      // use student likelihood (robust version), 
                                           // otherwise normal likelihood
  int<lower=1> t_fixedDF;                  // the fixed df parameter for the student likelihood
  int<lower=0, upper=1> modelMiss;         // model missingness (1) or not (0)
  int<lower=0, upper=1> NMAR;              // include NMAR component, yes (1) or no (0)
  vector[K] refIntensity;                  // ref intensity set in the NMAR term 
                                           // of the missingness model
  int<lower=0, upper=1> sigmaHierarch;     // hierarchical set-up for sigma (1) or not (0)
  int<lower=0, upper=1> sigmaRdmHierarch;  // hierarchical set-up for sigma of random effects (1) 
                                           // or not (0)
  int<lower=0, upper=1> gammaHierarch;     // hierarchical set-up for gamma (1) or not (0)
  
  // priors for intensity model
  real priorIntensityMean;            // prior distribution for intercept : mean
  real priorIntensitySd;              // prior distribution for intercept : std dev
  real priorFoldChangeMean;           // prior distribution for fold changes : mean
  real priorFoldChangeSd;             // prior distribution for fold changes : std dev
  real priorNuEpsilon;                // prior inv chi2 distribution for error variances : 
                                      // nu parameter
  real priorScaleEpsilon;             // prior inv chi2 distribution for error variances : 
                                      // scale parameter
  real priorNuRdm;                    // prior inv chi2 distribution for random effects variances : 
                                      // nu parameter
  real priorScaleRdm;                 // prior inv chi2 distribution for random effects variances : 
                                      // scale parameter
  
  // priors for missingness model
  real priorGammaInterceptMean;       // prior normal distribution for intercept 
                                      // in missingness model : mean
  real priorGammaInterceptSd;         // prior normal distribution for intercept 
                                      // in missingness model : std dev
  real priorGammaTermsMean[max(MGamma[1]-1,1)];// prior normal distribution for the other MAR coefficients
                                      // in missingness model : mean
  real priorGammaTermsSd[max(MGamma[1]-1,1)];  // prior normal distribution for the other MAR coefficients 
                                      // in missingness model : std dev
  real priorGammaYMean;               // prior normal distribution for NMAR coefficient 
                                      // in missingness model : mean
  real priorGammaYSd;                 // prior normal distribution for NMAR coefficient 
                                      // in missingness model : std dev
  
  // contrasts
  int<lower=0> nContrast[K];          // nb of contrasts per protein
  int<lower=0> contrastCoefLen;       // dimension of the flat vector containing all contrast coefficients 
  vector[contrastCoefLen] contrastCoefs;// concatenation of all constrast coefficient vectors

}


transformed data {
  int NMAlphaObs[K];
  int NMBetaObs[K];
  int NMGammaObs[K];
  int NMAlphaMis[modelMiss?K:0];
  int NMBetaMis[modelMiss?K:0];
  int NMGammaMis[modelMiss?K:0];
  int MAlphaBeta[K];
  int MMAlphaBeta[K];
  int MAlphaSum;
  int MBetaSum;
  int MGammaSum;
  int MGammaMax;
  int MAlphaBetaSum;
  int MMAlphaBetaSum;
  int NContrastSum;

  // store N*M for all k
  for(k in 1:K){
    NMAlphaObs[k] = NObs[k] * MAlpha[k];
    NMBetaObs[k] = NObs[k] * MBeta[k];
    NMGammaObs[k] = NObs[k] * MGamma[k];
    MAlphaBeta[k] = MAlpha[k]+MBeta[k];
    MMAlphaBeta[k] = (MAlpha[k]+MBeta[k]) * (MAlpha[k]+MBeta[k]);
    if(modelMiss){
      NMAlphaMis[k] = NMis[k] * MAlpha[k];
      NMBetaMis[k] = NMis[k] * MBeta[k];
      NMGammaMis[k] = NMis[k] * MGamma[k];
    } 
  }
  MAlphaSum = sum(MAlpha);
  MBetaSum = sum(MBeta);
  MGammaSum = sum(MGamma);
  MGammaMax = max(MGamma);
  MAlphaBetaSum = sum(MAlphaBeta);
  MMAlphaBetaSum = sum(MMAlphaBeta);
  NContrastSum = sum(nContrast);
}

parameters {
  // parameters of intensity model
  // alpha, beta are here fixed, depending on the value of the other parameters
  // variance parameters are parametrized in log()
  vector[K] logSigma;                             // log of error std dev
  vector[nRdm*K] logSigmaRdm;                     // log of std dev of random effects
  vector[t_likelihood*NObsSum] logBObs;           // log of b auxiliary variables 
                                                  // (for student t likelihood)
  vector[t_likelihood*modelMiss*NMisSum] logBMis; // log of b auxiliary variables 
                                                  // (for student t likelihood)
  
  // hyperparameters of hierarchical intensity model
  real<lower=1> nu0Eps[sigmaHierarch?1:0];        // degree of freedom hyperparameter 
                                                  // for error variances
  real<lower=0> s0Eps2[sigmaHierarch?1:0];        // scale hyperparameter for error variances
  real<lower=1> nu0Rdm[nRdm*sigmaRdmHierarch];    // DoF for variances of rdm effects
  real<lower=0> s0Rdm2[nRdm*sigmaRdmHierarch];    // scale for variances of rdm effects 
  
  // * parameters of missingness model *
  
  // coefficients for predictors in probit regression of missingness (MAR)
  vector[modelMiss?MGammaSum:0] gamma;            
  // coefficients for Y predictor in probit regression of missingness (NMAR)
  vector[modelMiss*NMAR?K:0] gammY; 
  // missing responses
  vector[modelMiss?NMisSum:0] yMis;  
  
  // * hyperparameters of hierarchical missingness model *
  
  // mean hyperparameters for gamma parameters
  vector [modelMiss*gammaHierarch?(MGammaMax+NMAR):0] mu0gm;  
  // scale hyperparameters for gamma parameters
  vector<lower=0> [modelMiss*gammaHierarch?(MGammaMax+NMAR):0] sigma20gm; 
}

transformed parameters{
  
  vector[MAlphaBetaSum] AllMuqAlBet;
  vector[MMAlphaBetaSum] AllInvSigqAlBet;
  
  {
    int yobspos;
    int zobspos;
    int xobspos;
    int tobspos;
    int ymispos;
    int zmispos;
    int xmispos;
    int tmispos;
    int apos;
    int bpos;
    int spos;
  
    int muAlphaBetapos;
    int invSigmaAlphaBetapos;
  
    yobspos = 1;
    xobspos = 1;
    zobspos = 1;
    tobspos = 1;
    ymispos = 1;
    xmispos = 1;
    zmispos = 1;
    tmispos = 1;
    apos = 1;
    bpos = 1;
    spos = 1;
  
    muAlphaBetapos = 1;
    invSigmaAlphaBetapos = 1;
  
  
    for(k in 1:K){
      vector[NObs[k]] YO;
      matrix[NObs[k], MAlpha[k]] XO;
      matrix[NObs[k], MBeta[k]] ZO;
      matrix[modelMiss?NObs[k]:0, modelMiss?MGamma[k]:0] TO;
    
      vector[modelMiss?NMis[k]:0] YM;
      matrix[modelMiss?NMis[k]:0, MBeta[k]] ZM;
      matrix[modelMiss?NMis[k]:0, modelMiss?MAlpha[k]:0] XM;
      matrix[modelMiss?NMis[k]:0, modelMiss?MGamma[k]:0] TM;
    
      matrix[modelMiss?(NObs[k]+NMis[k]):NObs[k], MAlpha[k] + MBeta[k]] W;
      matrix[NObs[k], MAlpha[k] + MBeta[k]] WO;
      matrix[modelMiss?NMis[k]:0, MAlpha[k] + MBeta[k]] WM;
    
      vector[NObs[k]] invBO;
      vector[modelMiss?NMis[k]:0] invBM;
      matrix[modelMiss?(NObs[k]+NMis[k]):NObs[k], modelMiss?(NObs[k]+NMis[k]):NObs[k]] D;
      matrix[NObs[k], NObs[k]]  DO;
      matrix[modelMiss?NMis[k]:0, modelMiss?NMis[k]:0] DM;
    
      vector[modelMiss?NObs[k]:0] aO;
      vector[modelMiss?NMis[k]:0] aM;
      vector[modelMiss?(NObs[k]+NMis[k]):0] a;
    
      real invSigma2Eps;
      real sigma2Eps;
      
      vector[nRdm] invSigma2Rdm;
      vector[nRdm] currentLogSigmaRdm;
      vector[nRdm] sigma2Rdm;
    
    
      vector[MAlpha[k]] currentAlpha;
      vector[MBeta[k]] currentBeta;
      vector[modelMiss?MGamma[k]:0] currentGamma;
      real currentGammY;
    
      vector[MAlpha[k]] alphaInvVarPriors;
      matrix[MAlpha[k]+MBeta[k], MAlpha[k]+MBeta[k]] invSigmaAlphaBeta;

      matrix[MAlphaBeta[k], MAlphaBeta[k]] invSigmaqAlphaBeta;

      vector[MAlpha[k]] alphaMeanPriors;
      vector[MAlpha[k] + MBeta[k]] muAlphaBeta;
      vector[MAlphaBeta[k]] muqAlphaBeta;
    
      // initialize vectors and matrices
      YO = segment(yObs, yobspos, NObs[k]);
      XO = to_matrix( segment(XObs, xobspos, NMAlphaObs[k]), NObs[k], MAlpha[k]);
    
      if(modelMiss){
        YM = segment(yMis, ymispos, NMis[k]);
        XM = to_matrix( segment(XMis, xmispos, NMAlphaMis[k]), NMis[k], MAlpha[k]);
        TO = to_matrix( segment(TObs, tobspos, NMGammaObs[k]), NObs[k], MGamma[k]);
        TM = to_matrix( segment(TMis, tmispos, NMGammaMis[k]), NMis[k], MGamma[k]);
      }
    
      if(nRdm > 0){
        ZO = to_matrix( segment(ZObs, zobspos, NMBetaObs[k]), NObs[k], MBeta[k]);
        WO = append_col(XO, ZO);
        if(modelMiss){
          ZM = to_matrix( segment(ZMis, zmispos, NMBetaMis[k]), NMis[k], MBeta[k]);
          WM = append_col(XM, ZM);
          W = append_row(WO, WM);
        }
        else{
          W = append_col(XO, ZO);
        }
      } else {
        WO = XO;
        if(modelMiss){
          WM = XM;
          W = append_row(WO, WM);
        } else {
          W = WO;
        }
      }
    
      invSigma2Eps = exp(-2*logSigma[k]);
      if(t_likelihood){
        invBO = exp(-segment(logBObs, yobspos, NObs[k]));
        if(modelMiss){
          invBM = exp(-segment(logBMis, ymispos, NMis[k]));
        }
      } else {
        invBO = rep_vector(1., NObs[k]);
        if(modelMiss){
          invBM = rep_vector(1., NMis[k]);
        }
      }
      DO = diag_matrix( invBO );
      if(modelMiss){
        DM = diag_matrix( invBM );
        D = diag_matrix(append_row(invBO, invBM));
      } else {
        D = DO;
      }
    
      if(nRdm > 0){
        currentLogSigmaRdm = segment(logSigmaRdm, spos, nRdm);
        invSigma2Rdm = exp(-2*currentLogSigmaRdm);
        sigma2Rdm = inv(invSigma2Rdm);
      }
    
      // reorder priors for alpha and beta => compute related matrices and vectors
      if(MAlpha[k] > 1){
        alphaInvVarPriors = append_row(1/(priorIntensitySd*priorIntensitySd),
        rep_vector(1/(priorFoldChangeSd*priorFoldChangeSd), MAlpha[k]-1));
        alphaMeanPriors = append_row(priorIntensityMean,
        rep_vector(priorFoldChangeMean, MAlpha[k]-1));
      } else {
        alphaInvVarPriors = rep_vector(1/(priorIntensitySd*priorIntensitySd), 1);
        alphaMeanPriors = rep_vector(priorIntensityMean, 1);
      }
      if(nRdm > 0){
        vector[MBeta[k]] betaInvVarPriors;
        int ipos;
        ipos = 1;
        for(r in 1:nRdm){
          for(i in ipos:(ipos+MBetaPerRdmEffect[spos+r-1]-1)){
            betaInvVarPriors[i] = invSigma2Rdm[r];
          }
          ipos += MBetaPerRdmEffect[spos+r-1];
        }
        invSigmaAlphaBeta = diag_matrix(append_row(alphaInvVarPriors, betaInvVarPriors));
        muAlphaBeta = append_row(alphaMeanPriors, rep_vector(0., MBeta[k]));
      } else {
        invSigmaAlphaBeta = diag_matrix(alphaInvVarPriors);
        muAlphaBeta = alphaMeanPriors;
      }
    
      // STEP 1 : 
      // precompute fixed values for mu and sigma of (alpha, beta, YMis)
      // and calculate logPost term linked to determinant()
      invSigmaqAlphaBeta = invSigma2Eps * W'* D * W + invSigmaAlphaBeta;
      muqAlphaBeta = invSigma2Eps * W'* D * append_row(YO, YM) + invSigmaAlphaBeta * muAlphaBeta;
      muqAlphaBeta = invSigmaqAlphaBeta \ muqAlphaBeta; // A \ b = inverse(A) * b

      // store mu and Sigma for alpha, beta, YMis
      AllMuqAlBet[muAlphaBetapos:(muAlphaBetapos+MAlphaBeta[k]-1)] = 
        muqAlphaBeta;
      AllInvSigqAlBet[invSigmaAlphaBetapos:(invSigmaAlphaBetapos+MMAlphaBeta[k]-1)] =
        to_vector(invSigmaqAlphaBeta);
    
      // update postion indexes  
      apos += MAlpha[k];
      bpos += MBeta[k];
      yobspos += NObs[k];
      xobspos += NMAlphaObs[k];
      zobspos += NMBetaObs[k];
      spos += nRdm;
    
      if(modelMiss){
        ymispos += NMis[k];
        xmispos += NMAlphaMis[k];
        zmispos += NMBetaMis[k];
        tobspos += NMGammaObs[k];
        tmispos += NMGammaMis[k];
      }
    
      muAlphaBetapos += MAlphaBeta[k];
      invSigmaAlphaBetapos += MMAlphaBeta[k];
    }
  }
}
  
model {
  int yobspos;
  int zobspos;
  int xobspos;
  int tobspos;
  int ymispos;
  int zmispos;
  int xmispos;
  int tmispos;
  int apos;
  int bpos;
  int cpos;
  int spos;
  int muAlphaBetapos;
  int invSigmaAlphaBetapos;
  
  
  
  yobspos = 1;
  xobspos = 1;
  zobspos = 1;
  tobspos = 1;
  ymispos = 1;
  xmispos = 1;
  zmispos = 1;
  tmispos = 1;
  apos = 1;
  bpos = 1;
  cpos = 1;
  spos = 1;
  
  muAlphaBetapos = 1;
  invSigmaAlphaBetapos = 1;
  
  
  for(k in 1:K){
    vector[NObs[k]] YO;
    matrix[NObs[k], MAlpha[k]] XO;
    matrix[NObs[k], MBeta[k]] ZO;
    matrix[modelMiss?NObs[k]:0, modelMiss?MGamma[k]:0] TO;
    
    vector[modelMiss?NMis[k]:0] YM;
    matrix[modelMiss?NMis[k]:0, MBeta[k]] ZM;
    matrix[modelMiss?NMis[k]:0, modelMiss?MAlpha[k]:0] XM;
    matrix[modelMiss?NMis[k]:0, modelMiss?MGamma[k]:0] TM;
    
    vector[NObs[k]] invBO;
    vector[modelMiss?NMis[k]:0] invBM;
    
    vector[NObs[k]] BO;
    vector[modelMiss?NMis[k]:0] BM;
    
    real invSigma2Eps;
    real sigma2Eps;
    vector[NObs[k]] sigmaEpsO;
    vector[modelMiss?NMis[k]:0] sigmaEpsM;
    vector[nRdm] invSigma2Rdm;
    vector[nRdm] currentLogSigmaRdm;
    vector[nRdm] sigma2Rdm;
    
    
    vector[MAlpha[k]] currentAlpha;
    vector[MBeta[k]] currentBeta;
    vector[modelMiss?MGamma[k]:0] currentGamma;
    real currentGammY;

    matrix[MAlphaBeta[k], MAlphaBeta[k]] invSigmaqAlphaBeta;

    vector[MAlphaBeta[k]] muqAlphaBeta;
    
    // initialize vectors and matrices
    YO = segment(yObs, yobspos, NObs[k]);
    XO = to_matrix( segment(XObs, xobspos, NMAlphaObs[k]), NObs[k], MAlpha[k]);
    
    if(modelMiss){
      YM = segment(yMis, ymispos, NMis[k]);
      XM = to_matrix( segment(XMis, xmispos, NMAlphaMis[k]), NMis[k], MAlpha[k]);
      TO = to_matrix( segment(TObs, tobspos, NMGammaObs[k]), NObs[k], MGamma[k]);
      TM = to_matrix( segment(TMis, tmispos, NMGammaMis[k]), NMis[k], MGamma[k]);
    }
    
    if(nRdm > 0){
      ZO = to_matrix( segment(ZObs, zobspos, NMBetaObs[k]), NObs[k], MBeta[k]);
      if(modelMiss){
        ZM = to_matrix( segment(ZMis, zmispos, NMBetaMis[k]), NMis[k], MBeta[k]);
      }
    } 
    
    invSigma2Eps = exp(-2*logSigma[k]);
    if(t_likelihood){
      invBO = exp(-segment(logBObs, yobspos, NObs[k]));
      BO = inv(invBO);
      if(modelMiss){
        invBM = exp(-segment(logBMis, ymispos, NMis[k]));
        BM = inv(invBM);
      }
    } else {
      invBO = rep_vector(1., NObs[k]);
      if(modelMiss){
        invBM = rep_vector(1., NMis[k]);
      }
    }
    
    sigmaEpsO = sqrt(inv(invBO * invSigma2Eps));
    if(modelMiss){
      sigmaEpsM = sqrt(inv(invBM * invSigma2Eps));
    }
    sigma2Eps = inv(invSigma2Eps);
    
    if(nRdm > 0){
      currentLogSigmaRdm = segment(logSigmaRdm, spos, nRdm);
      invSigma2Rdm = exp(-2*currentLogSigmaRdm);
      sigma2Rdm = inv(invSigma2Rdm);
    }
    
    if(modelMiss){
      currentGamma = segment(gamma, cpos, MGamma[k]);
      if(NMAR){
        currentGammY = gammY[k];
      } else {
        currentGammY = 0.;
      }
    }
    
    muqAlphaBeta = AllMuqAlBet[muAlphaBetapos:(muAlphaBetapos+MAlphaBeta[k]-1)];
    invSigmaqAlphaBeta = to_matrix(
        AllInvSigqAlBet[invSigmaAlphaBetapos:(invSigmaAlphaBetapos+MMAlphaBeta[k]-1)],
        MAlphaBeta[k], MAlphaBeta[k]);
    
    currentAlpha = segment(muqAlphaBeta, 1, MAlpha[k]);
    if(nRdm > 0){
      currentBeta = segment(muqAlphaBeta, MAlpha[k]+1, MBeta[k]);
    }
    
    target += -0.5 * log(determinant(invSigmaqAlphaBeta));
  
    // STEP 2 : now calculate the remaining terms of the log posterior 
    
    // likelihood of observed responses
    
    
    if(nRdm > 0){
      
      
      int rpos = 1;
      
      YO ~ normal( XO * currentAlpha + ZO * currentBeta, sigmaEpsO );
      
      // likelihood for current Beta
      for(r in 1:nRdm){
        int MCol = MBetaPerRdmEffect[(k-1)*nRdm+r];
        segment(currentBeta, rpos, MCol) ~ 
            normal(0, sqrt(sigma2Rdm[r]));
        rpos += MCol;
      }
      
    } else {
      YO ~ normal( XO * currentAlpha, sigmaEpsO );
    }
    
    
    
    
    if (modelMiss == 1){
      
      int Obs[NObs[k]];
      int Mis[NMis[k]];
      
      Obs = rep_array(1, NObs[k]);
      Mis = rep_array(0, NMis[k]);
    
      // missing responses distributions
      
      if(nRdm > 0){
        YM ~ normal( XM * currentAlpha + ZM * currentBeta, sigmaEpsM);
        
        // NOTE : currentBeta already has a prior defined (in likelihood section of observed data)
      } else {
        YM ~ normal( XM * currentAlpha, sigmaEpsM);
      }
      
      // likelihood of non missingness of observed responses
      
      { vector[NObs[k]] temp;
        if(NMAR){
          temp = Phi_approx( TO * currentGamma + (YO - refIntensity[k]) * gammY[k] );
        } else {
          temp = Phi_approx( TO * currentGamma );
        }
        for(n in 1:NObs[k]){
          Obs[n] ~ bernoulli(fmin(fmax(temp[n], 0.000001), 0.999999));
        }
      }
      
    
      // likelihood of missingness of missing responses
      { vector[NMis[k]] temp;
        if(NMAR){
          temp = Phi_approx( TM * currentGamma + (YM - refIntensity[k]) * gammY[k] );
        } else {
          temp = Phi_approx( TM * currentGamma );
        }
        for(n in 1:NMis[k]){
          Mis[n] ~ bernoulli(fmin(fmax(temp[n], 0.000001), 0.999999));
        }
      }
      
      
      // priors for gamma and gammaY
      if(gammaHierarch == 1){
        segment(gamma, cpos, MGamma[k]) ~ normal( mu0gm[1:MGamma[k]], sqrt(sigma20gm[1:MGamma[k]]) );
        
        // NMAR coefficient
        if(NMAR){
          gammY[k] ~ normal(mu0gm[MGamma[k]+1], sqrt(sigma20gm[MGamma[k]+1]));
        }
      } else {
        // intercept
        gamma[cpos] ~ normal(priorGammaInterceptMean, priorGammaInterceptSd);
        // other MAR coefficients
        if(MGamma[k] > 1){
          segment(gamma, cpos+1, MGamma[k]-1) ~ normal( priorGammaTermsMean, priorGammaTermsSd);
        }
        // NMAR coefficient
        if(NMAR){
          gammY[k] ~ normal(priorGammaYMean, priorGammaYSd);
        }
      }
      
    } // end if(modelMiss)
    
    // priors for sigma
    if(sigmaHierarch == 1){
      sigma2Eps ~ scaled_inv_chi_square(nu0Eps[1], sqrt(s0Eps2[1]));
    } else {
      sigma2Eps ~ scaled_inv_chi_square(priorNuEpsilon, priorScaleEpsilon);
    }
    
    if(nRdm > 0){
      if(sigmaRdmHierarch == 1){
        vector[nRdm] currentSigma2Rdm;
        spos = 1;
        sigma2Rdm ~ scaled_inv_chi_square(nu0Rdm, sqrt(s0Rdm2));
      } else {
        sigma2Rdm ~ scaled_inv_chi_square(priorNuRdm, priorScaleRdm);
      }
    }
    
    // priors for b
    if(t_likelihood){
      BO ~ inv_gamma(t_fixedDF/2., t_fixedDF/2.);
      if(modelMiss){
        BM ~ inv_gamma(t_fixedDF/2., t_fixedDF/2.);
      }
    }
    
    // STEP 3 : add additional terms due to change of variables :
    // log(sigmaEps), log(sigmaRdm), log(b)
    target += log(sigma2Eps);
    if(nRdm > 0){
      target += sum(log(sigma2Rdm));
    }
    
    if(t_likelihood){
      target += - sum(log(invBO));
      if(modelMiss){
        target += - sum(log(invBM));
      }
    }
    
    
    // STEP 4 : update position indicators
    apos += MAlpha[k];
    bpos += MBeta[k];
    cpos += MGamma[k];
    yobspos += NObs[k];
    xobspos += NMAlphaObs[k];
    zobspos += NMBetaObs[k];
    spos += nRdm;
    
    if(modelMiss){
      ymispos += NMis[k];
      xmispos += NMAlphaMis[k];
      zmispos += NMBetaMis[k];
      tobspos += NMGammaObs[k];
      tmispos += NMGammaMis[k];
    }
    
    muAlphaBetapos += MAlphaBeta[k];
    invSigmaAlphaBetapos += MMAlphaBeta[k];
  } // end loop on proteins
  
  // hyper parameters priors => no priors !!
  
}

generated quantities {
  
  vector[MAlphaSum] alpha;                           // coefficients for fixed effects
  vector[nRdm>0?MBetaSum:0] beta;                    // coefficients for random effects
  
  vector[K] sigmaEps;
  vector[K] tauEps;
  vector[K] sigma2Eps;
  vector [nRdm*K] sigmaRdm;
  vector [nRdm*K] tauRdm;
  vector [nRdm*K] sigma2Rdm;
  vector [t_likelihood*NObsSum] invBObs;
  vector [t_likelihood*modelMiss*NMisSum] invBMis;
  
  // contrasts
  vector [NContrastSum] contrast;
  
  
  // generate samples for alpha, beta, yMis
  {
    
    int apos;
    int bpos;
    int ymispos;
    int muAlphaBetapos;
    int invSigmaAlphaBetapos;
    
    apos = 1;
    bpos = 1;
    ymispos = 1;
    
    muAlphaBetapos = 1;
    invSigmaAlphaBetapos = 1;
    
    for(k in 1:K){
      vector[MAlphaBeta[k]] muqAlphaBeta; 
      vector[MAlphaBeta[k]] alphaBeta;
      matrix[MAlphaBeta[k], MAlphaBeta[k]] invSigmaqAlphaBeta;
      
      muqAlphaBeta = AllMuqAlBet[muAlphaBetapos:(muAlphaBetapos+MAlphaBeta[k]-1)];
      
      invSigmaqAlphaBeta = to_matrix(
        AllInvSigqAlBet[invSigmaAlphaBetapos:(invSigmaAlphaBetapos+MMAlphaBeta[k]-1)],
        MAlphaBeta[k], MAlphaBeta[k]);
      
      alphaBeta = multi_normal_rng(muqAlphaBeta, inverse(invSigmaqAlphaBeta));
      //alphaBeta = muqAlphaBeta;
      alpha[apos:(apos+MAlpha[k]-1)] = alphaBeta[1:MAlpha[k]];
      if(nRdm > 0){
        beta[bpos:(bpos+MBeta[k]-1)] = alphaBeta[(MAlpha[k]+1):(MAlpha[k]+MBeta[k])];
      }
      
      apos += MAlpha[k];
      bpos += MBeta[k];
      if(modelMiss){
        ymispos += NMis[k];
      }
      muAlphaBetapos += MAlphaBeta[k];
      invSigmaAlphaBetapos += MMAlphaBeta[k];
    }
  }
  
  sigmaEps = exp(logSigma);
  sigma2Eps = sigmaEps .* sigmaEps;
  tauEps = inv(sigma2Eps);
  if(nRdm>0){
    sigmaRdm = exp(logSigmaRdm);
    sigma2Rdm = sigmaRdm .* sigmaRdm;
    tauRdm = inv(sigma2Rdm);
  }
  if(t_likelihood){
    invBObs = exp(-logBObs);
    if(modelMiss){
      invBMis = exp(-logBMis);
    }
    
  }
  
  {
    int aPos;
    int ctrPos;
    int ctrCoefPos;
    aPos = 1;
    ctrCoefPos = 1;
    ctrPos = 1;
    for(k in 1:K){
      if(nContrast[k]>0){
        for(c in 1:nContrast[k]){
          contrast[ctrPos] =
            sum(alpha[aPos:(aPos+MAlpha[k]-1)] .* 
                contrastCoefs[ctrCoefPos:(ctrCoefPos+MAlpha[k]-1)]);
          ctrCoefPos += MAlpha[k];
          ctrPos += 1;
        }
      }
      
      aPos += MAlpha[k];
    }
  }
  
}
