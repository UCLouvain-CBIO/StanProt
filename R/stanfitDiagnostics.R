


##' @title Check transitions that ended with a divergence
##' 
##' @param object a `stanfit` object
##' 
##' @import rstan
##' 
##' @export
check_div <- function(object) {
  stopifnot(inherits(object, "stanfit"))
  sampler_params <- 
    get_sampler_params(object, inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  n = sum(divergent)
  N = length(divergent)
  
  msg <- sprintf('%s of %s iterations ended with a divergence (%s%%)',
                 n, N, 100 * n / N)
  
  if (n > 0){
    msg <- paste0(msg, "\n",
                  "  Try running with larger adapt_delta to remove the divergences")
  }
  
  return(list(isProb = (n > 0),
              msg = msg))
}

##' @title Check transitions that ended prematurely due to maximum tree depth limit
##' 
##' @param object a `stanfit` object
##' 
##' @import rstan
##' 
##' @export
check_treedepth <- function(object) {
  stopifnot(inherits(object, "stanfit"))
  
  max_depth <-
    attributes(object)$stan_args[[1]]$control$max_treedepth
  if(is.null(max_depth))max_depth <- 10 # by default
  sampler_params <- get_sampler_params(object, inc_warmup=FALSE)
  treedepths <- do.call(rbind, sampler_params)[,'treedepth__']
  n = length(treedepths[sapply(treedepths, function(x) x == max_depth)])
  N = length(treedepths)
  
  msg <- sprintf('%s of %s iterations saturated the maximum tree depth of %s (%s%%)',
                 n, N, max_depth, 100 * n / N)
  
  if (n > 0){
    msg <- paste0(msg, "\n",
                  "  Run again with max_depth set to a larger value to avoid saturation")
  }
    
  return(list(isProb = (n > 0),
              msg = msg))
  
}

##' @title Check energy Bayesian fraction of missing information (E-BFMI)
##' 
##' @param object a `stanfit` object
##' 
##' @import rstan
##' 
##' @export
check_energy <- function(object) {
  stopifnot(inherits(object, "stanfit"))
  
  sampler_params <- get_sampler_params(object,
                                       inc_warmup=FALSE)
  no_warning <- TRUE
  msg <- ""
  for (n in 1:length(sampler_params)) {
    currentMsg <- ""
    energies = sampler_params[n][[1]][,'energy__']
    numer = sum(diff(energies)**2) / length(energies)
    denom = stats::var(energies)
    if (numer / denom < 0.2) {
      currentMsg <- sprintf('Chain %s: E-BFMI = %s', n, numer / denom)
      no_warning <- FALSE
    }
    if(currentMsg != ""){
      if(msg != ""){
        msg <- paste0(msg,"\n")
      }
      msg <- paste0(msg, currentMsg)
    }
    
  }
  if (no_warning){
    msg <- "E-BFMI indicated no pathological behavior"
  }
  else{
    currentMsg <- "  E-BFMI below 0.2 indicates you may need to reparameterize your model"
    msg <- paste0(msg, "\n", currentMsg)
  }
      
  return(list(isProb = !no_warning,
              msg = msg))
}


##' @title Check effective sample size per iteration
##' 
##' @param object a `stanfit` object
##' 
##' @import rstan
##' 
##' @export
check_n_eff <- function(object) {
  stopifnot(inherits(object, "stanfit"))
  
  # extract median (probs = 0.5) of all parameters (use cache for future use)
  fit_summary <- summary(object, probs = c(0.5))$summary
  N <- dim(fit_summary)[[1]]
  nchains <- length(attr(object,"stan_args"))
  n_eff_thresh <- 20 * nchains; # min 10 * 2 splits otherwise Rhat diagnostic is useless
  
  iter <- dim(extract(object)[[1]])[[1]]
  
  nb_warning_ratio <- 0
  nb_warning_absolute <- 0
  msg <- ""
  for (n in 1:N) {
    n_eff <- fit_summary[,5][n]
    ratio <- n_eff / iter
    currentMsg <- ""
    if (ratio < 0.001) {
      currentMsg <- 
        sprintf("n_eff / iter for parameter %s is %s!",
                rownames(fit_summary)[n], ratio)
      nb_warning_ratio <- nb_warning_ratio + 1
    }
    if (n_eff < n_eff_thresh){
      if(currentMsg != ""){
        currentMsg <- paste0(currentMsg,";")
      }
      currentMsg <- paste0(currentMsg,
                           sprintf("n_eff for parameter %s is %s!",
                                   rownames(fit_summary)[n], n_eff))
      nb_warning_absolute <- nb_warning_absolute + 1
    }
    if( currentMsg != ""){
      if(msg != ""){
        msg <- paste0(msg, "\n")
      }
      msg <- paste0(msg, currentMsg)
    }
  }
  if (nb_warning_ratio == 0){
    if(msg != ""){
      msg <- paste0(msg, "\n")
    }
    msg <- paste0(msg,
                  "n_eff / iter looks reasonable for all parameters")
  } else {
    msg <- 
      paste0(msg, 
             "  n_eff / iter below 0.001 indicates that the effective sample ",
             "size has likely been overestimated")
  }
    
  if (nb_warning_absolute == 0){
    if(msg != ""){
      msg <- paste0(msg, "\n")
    }
    msg <- paste0(msg,
                  "n_eff looks reasonable for all parameters")
  } else {
    msg <- 
      paste0(msg,
             sprintf("  n_eff / nChains below 100 for %s parameters (%s%%).",
                  nb_warning_absolute, 100*nb_warning_absolute/N),
      "\nThis indicates that for these parameters the effective sample size is ",
      "too small to provide reliable Rhat diagnostic")
    
  }
  
  return(list(isProb = (nb_warning_ratio>0 || nb_warning_absolute>0), 
              msg = msg))
}

##' @title Checks potential scale reduction factors
##' 
##' @param object a `stanfit` object
##' 
##' @import rstan
##' 
##' @export
check_rhat <- function(object) {
  stopifnot(inherits(object, "stanfit"))
  
    rhat_thresh <- 1.1
    fit_summary <- summary(object, probs = c(0.5))$summary
    N <- dim(fit_summary)[[1]]
    
    msg <- ""
    no_warning <- TRUE
    for (n in 1:N) {
      rhat <- fit_summary[,6][n]
      if (rhat > rhat_thresh || is.infinite(rhat) || is.nan(rhat)) {
        currentMsg <- sprintf("Rhat for parameter %s is %s!",
                              rownames(fit_summary)[n], rhat)
        no_warning <- FALSE
        if(msg != ""){
          msg <- paste0(msg,"\n")
        }
        msg <- paste0(msg, currentMsg)
      }
    }
    if (no_warning){
      msg <- "Rhat looks reasonable for all parameters"
    } else {
      msg <- paste0(msg, "\n",
                    "  Rhat above ", rhat_thresh,
                    " indicates that the chains very likely have not mixed"
                    )
    }
    
    return(list(isProb = !no_warning,
                msg = msg))
  
}


##' @title Run all diagnostics one at a time
##' 
##' @param object a `stanfit` object
##' 
##' @import rstan
##' 
##' @export
check_all_diagnostics <- function(object) {
  stopifnot(inherits(object, "stanfit"))
  resNeff <- check_n_eff(object)
  resRhat <- check_rhat(object)
  resDiv <- check_div(object)
  resTreeDepth <- check_treedepth(object)
  resEnergy <- check_energy(object)
  
  isProb <- resNeff$isProb ||
    resRhat$isProb ||
    resDiv$isProb ||
    resTreeDepth$isProb ||
    resEnergy$isProb
  
  msg <- paste0(resNeff$msg, "\n",
                resRhat$msg, "\n",
                resDiv$msg, "\n",
                resTreeDepth$msg, "\n",
                resEnergy$msg)
  
  return(list(isProb = isProb,
              msg = msg))
}

##' @title Returns parameter arrays separated into divergent and non-divergent transitions
##'
##' @param object a `stanfit` object
##' 
##' @import rstan
##' 
##' @export
partition_div <- function(object) {
  stopifnot(inherits(object, "stanfit"))
  
  nom_params <- extract(object, permuted=FALSE)
  n_chains <- dim(nom_params)[2]
  params <- as.data.frame(do.call(rbind,
                                  lapply(1:n_chains,
                                         function(n) nom_params[,n,]))
  )

  sampler_params <- get_sampler_params(object,
                                       inc_warmup=FALSE)
  divergent <- do.call(rbind, sampler_params)[,'divergent__']
  params$divergent <- divergent

  div_params <- params[params$divergent == 1,]
  nondiv_params <- params[params$divergent == 0,]

  return(list(div_params, nondiv_params))
}


##' @title Returns total elapsed time that was needed to run stan
##'
##' @param object a `stanfit` object
##' 
##' @import rstan
##' 
##' @export
get_total_elapsed_time <- function(object) {
  stopifnot(inherits(object, "stanfit"))
  
  tmpMat <- get_elapsed_time(object)
  elapsed <- sum(rowSums(tmpMat))
  elapsed
}




