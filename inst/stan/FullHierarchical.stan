// saved as FullHierarchical.stan

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
  int<lower=0, upper=1> NMAR;              // include NMAR component, yes (1) or not (0)
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
  real priorGammaTermsMean[max(MGamma[1]-1,1)];// prior normal distribution for the other MAR 
                                      // coefficients in missingness model : mean
  real priorGammaTermsSd[max(MGamma[1]-1,1)];  // prior normal distribution for the other MAR  
                                      // coefficients in missingness model : std dev
  real priorGammaYMean;               // prior normal distribution for NMAR coefficient 
                                      // in missingness model : mean
  real priorGammaYSd;                 // prior normal distribution for NMAR coefficient 
                                      // in missingness model : std dev
  
  // contrasts
  int<lower=0> nContrast[K];          // nb of contrasts per protein
  int<lower=0> contrastCoefLen;       // dimension of the flat vector containing all contrast 
                                      // coefficients 
  vector[contrastCoefLen] contrastCoefs;// concatenation of all constrast coefficient vectors
}


transformed data {
  int NMAlphaObs[K];
  int NMBetaObs[K];
  int NMGammaObs[K];
  int NMAlphaMis[modelMiss?K:0];
  int NMBetaMis[modelMiss?K:0];
  int NMGammaMis[modelMiss?K:0];
  int MAlphaSum;
  int MBetaSum;
  int MGammaSum;
  int MGammaMax;
  int NContrastSum;

  // store N*M for all k
  for(k in 1:K){
    NMAlphaObs[k] = NObs[k] * MAlpha[k];
    NMBetaObs[k] = NObs[k] * MBeta[k];
    NMGammaObs[k] = NObs[k] * MGamma[k];
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
  NContrastSum = sum(nContrast);
}

parameters {
  
  // * parameters of intensity model *
  
  // coefficients for fixed effects
  vector[MAlphaSum] alpha;    
  // coefficients for random effects
  vector[nRdm>0?MBetaSum:0] beta;  
  // error variances
  real<lower=0> sigma2Eps[K];   
  // variances of random effects
  vector<lower=0> [nRdm*K] sigma2Rdm;  
         
  
  // * hyperparameters of hierarchical intensity model *
  
  // degrees of freedom hyperparameter for error variances
  real<lower=1> nu0Eps[sigmaHierarch?1:0];     
  // scale hyperparameter for error variances
  real<lower=0> s0Eps2[sigmaHierarch?1:0];    
  // degrees of freedom hyperparameter for variances of rdm effects
  real<lower=1> nu0Rdm[nRdm*sigmaRdmHierarch];   
  // scale hyperparameter for variances of rdm effects 
  real<lower=0> s0Rdm2[nRdm*sigmaRdmHierarch];    
  
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
  
  
  for(k in 1:K){
    vector[MAlpha[k]] currentAlpha;
    matrix[NObs[k], MAlpha[k]] XO;
    
    vector[NObs[k]] YO;
    
    currentAlpha = segment(alpha, apos, MAlpha[k]);
    
    // likelihood of observed responses
    
    YO = segment(yObs, yobspos, NObs[k]);
    XO = to_matrix( segment(XObs, xobspos, NMAlphaObs[k]), NObs[k], MAlpha[k]);
    
    if(nRdm > 0){
      vector[MBeta[k]] currentBeta;
      vector[nRdm] currentSigma2Rdm;
      matrix[NObs[k], MBeta[k]] ZO;
      int rpos = 1;
    
      currentBeta = segment(beta, bpos, MBeta[k]);
      currentSigma2Rdm = segment(sigma2Rdm, spos, nRdm);
      ZO = to_matrix( segment(ZObs, zobspos, NMBetaObs[k]), NObs[k], MBeta[k]);
      
      if(t_likelihood){
        YO ~ student_t(t_fixedDF, XO * currentAlpha + ZO * currentBeta, sqrt(sigma2Eps[k]) );
      } else {
        YO ~ normal( XO * currentAlpha + ZO * currentBeta, sqrt(sigma2Eps[k]) );
      }
      
      
      for(r in 1:nRdm){
        int MCol = MBetaPerRdmEffect[(k-1)*nRdm+r];
        segment(currentBeta, rpos, MCol) ~ 
            normal(0, sqrt(currentSigma2Rdm[r]));
        rpos += MCol;
      }
      
    } else {
      if(t_likelihood){
        YO ~ student_t(t_fixedDF, XO * currentAlpha, sqrt(sigma2Eps[k]) );
      } else {
        YO ~ normal( XO * currentAlpha, sqrt(sigma2Eps[k]) );
      }
      
    }
    
    // priors for Alpha's
    // intercept
    alpha[apos] ~ normal(priorIntensityMean, priorIntensitySd);
    // logFoldChanges all
    if(MAlpha[k] > 1){
      segment(alpha, apos+1, MAlpha[k]-1) ~ normal( priorFoldChangeMean, priorFoldChangeSd);
    }
    
    yobspos += NObs[k];
    xobspos += NMAlphaObs[k];
    zobspos += NMBetaObs[k];
    
    
    if (modelMiss == 1){
      
      matrix[modelMiss?NMis[k]:0, modelMiss?MAlpha[k]:0] XM;
      matrix[modelMiss?NObs[k]:0, modelMiss?MGamma[k]:0] TO;
      matrix[modelMiss?NMis[k]:0, modelMiss?MGamma[k]:0] TM;
      vector[MGamma[k]] currentGamma;
      
      vector[NMis[k]] YM;
      
      int Obs[NObs[k]];
      int Mis[NMis[k]];
    
      // missing responses distributions
      YM = segment(yMis, ymispos, NMis[k]);
      XM = to_matrix( segment(XMis, xmispos, NMAlphaMis[k]), NMis[k], MAlpha[k]);
      
      if(nRdm > 0){
        vector[MBeta[k]] currentBeta;
        matrix[NMis[k], MBeta[k]] ZM;
        
        currentBeta = segment(beta, bpos, MBeta[k]);
        
        ZM = to_matrix( segment(ZMis, zmispos, NMBetaMis[k]), NMis[k], MBeta[k]);
    
        if(t_likelihood){
          YM ~ student_t( t_fixedDF, XM * currentAlpha + ZM * currentBeta, sqrt(sigma2Eps[k]) );
        } else {
          YM ~ normal( XM * currentAlpha + ZM * currentBeta, sqrt(sigma2Eps[k]) );
        }
        
        // NOTE : currentBeta already has a prior defined (in likelihood section of observed data)
      } else {
        if(t_likelihood){
          YM ~ student_t( t_fixedDF, XM * currentAlpha, sqrt(sigma2Eps[k]));
        } else {
          YM ~ normal( XM * currentAlpha, sqrt(sigma2Eps[k]));
        }
        
      }
      
      // likelihood of non missingness of observed responses
      currentGamma = segment(gamma, cpos, MGamma[k]);
      TO = to_matrix( segment(TObs, tobspos, NMGammaObs[k]), NObs[k], MGamma[k]);
      Obs = rep_array(1, NObs[k]);

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
      TM = to_matrix( segment(TMis, tmispos, NMGammaMis[k]), NMis[k], MGamma[k]);
      Mis = rep_array(0, NMis[k]);
      
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
      
      
      cpos += MGamma[k];
      ymispos += NMis[k];
      xmispos += NMAlphaMis[k];
      zmispos += NMBetaMis[k];
      tobspos += NMGammaObs[k];
      tmispos += NMGammaMis[k];
      
    } // end if(modelMiss)
    
    apos += MAlpha[k];
    bpos += MBeta[k];
    spos += nRdm;
  } // end loop on proteins
  
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
      for(k in 1:K){
        currentSigma2Rdm = segment(sigma2Rdm, spos, nRdm);
        currentSigma2Rdm ~ scaled_inv_chi_square(nu0Rdm, sqrt(s0Rdm2));
        spos += nRdm;
      }
      
    } else {
      sigma2Rdm ~ scaled_inv_chi_square(priorNuRdm, priorScaleRdm);
    }
  }
  
  
  // hyper parameters priors => no priors !!
  

}

generated quantities {
 
    //other ways to express sigma2 parameters
    real sigmaEps[K];
    real tauEps[K];
    vector [nRdm*K] sigmaRdm;
    vector [nRdm*K] tauRdm;
    vector [modelMiss*gammaHierarch?(MGammaMax+NMAR):0] tau0gm;
    vector [modelMiss*gammaHierarch?(MGammaMax+NMAR):0] sigma0gm;
    
    // contrasts
    vector [NContrastSum] contrast;
    
    
    sigmaEps = sqrt(sigma2Eps);
    tauEps = inv(sigma2Eps);
    if(nRdm>0){
      sigmaRdm = sqrt(sigma2Rdm);
      tauRdm = inv(sigma2Rdm);
    }
    if(modelMiss * gammaHierarch){
      tau0gm = inv(sigma20gm);
      sigma0gm = sqrt(sigma20gm);
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
