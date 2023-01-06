##' @export
logdNu <- function(x, xlogxcoef, lgammaxcoef, xcoef, const){
  value <- ifelse(x<=0., 0., xlogxcoef * x * log(x) + lgammaxcoef * lgamma(x) + 
                    xcoef * x + const)
}

##' @export
dNu <- function(x, xlogxcoef, lgammaxcoef, xcoef, const){
  value <- ifelse(x<=0, 0., exp( logdNu(x, xlogxcoef, lgammaxcoef, xcoef, const)))
  for(i in 1:length(value)){
    if(!is.finite(value[i])){
      msg <- paste0("not finite value returned in density function (x = ", x[i])
      msg <- paste0(msg, " ; xlogxcoef = ", xlogxcoef)
      msg <- paste0(msg, " ; lgammaxcoef = ", lgammaxcoef)
      msg <- paste0(msg, " ; xcoef = ", xcoef)
      msg <- paste0(msg, " ; const = ", const)
      msg <- paste0(msg, " ; value = ", value[i],")\n")
      stop(msg)
    }
  }
  value
}

##' @export
dANu <- function(x, A, xlogxcoef, lgammaxcoef, xcoef, const){
  1/A * dNu(x/A, xlogxcoef, lgammaxcoef, xcoef, const)
}

##' @export
pANu <- function(q, A, xlogxcoef, lgammaxcoef, xcoef, const){
  p <- rep(0., length(q))
  for(i in seq_along(q)){
    p[i] <- integrate(f = dANu, lower = 0, upper = q,
                      A = A,
                      xlogxcoef = xlogxcoef,
                      lgammaxcoef = lgammaxcoef,
                      xcoef = xcoef,
                      const = const)$value
  }
  p
}

##' @export
qANu <- function(p, A, xlogxcoef, lgammaxcoef, xcoef, const){
  q <- rep(0., length(p))
  for(i in seq_along(p)){
    reverseSignFound <- FALSE
    upper <- 1
    for (j in 1:12){
      if(pANu(q = upper, A, xlogxcoef, lgammaxcoef, xcoef, const) > p[i]){
        reverseSignFound <- TRUE
        break
      }
      upper <- upper * 10
    }
    if(!reverseSignFound){
      stop("error in qANu : no sign reversal found in interval [0:1e10]")
    }
    q[i] <- uniroot(
      f = function(x){
        pANu(x, A, xlogxcoef, lgammaxcoef, xcoef, const) - p[i]
      },
      lower = 0.,
      upper = upper,
      extendInt = "no")$root
  }
  q
}

##' @export
rANu <- function(n, A, xlogxcoef, lgammaxcoef, xcoef, const){
  p <- runif(n)
  sample <- qANu(p, A, xlogxcoef, lgammaxcoef, xcoef, const)
}

##' @export
integrateDANu <- function(A, xlogxcoef, lgammaxcoef, xcoef, const, 
                          lowerBound, upperBound, log = FALSE){
  if(upperBound <= 0.) return (0.)
  if(lowerBound > 0.) return(
    integratePANu(p, A, xlogxcoef, lgammaxcoef, xcoef, const, 
                  lowerBound = 0, upperBound = upperBound) -
      integratePANu(p, A, xlogxcoef, lgammaxcoef, xcoef, const, 
                    lowerBound = 0, upperBound = lowerBound))
  if(lowerBound < 0.) lowerBound = 0.
  
  # first find where the maximum of the function is
  xMax <- 0.
  maxFound <- FALSE
  reverseSignFound <- FALSE
  if(upperBound)
  if(upperBound == Inf){
    upper <- 1
    for (j in 1:13){
      if(xlogxcoef * (log(upper)+1) + lgammaxcoef * digamma(upper) + xcoef < 0){
        reverseSignFound <- TRUE
        break
      }
      upper <- upper * 10
    }
    if(!reverseSignFound){
      stop("error in integrateDANu : no maximum found in interval [0:1e12]")
    }
  } else {
    upper <- upperBound*A
    if(xlogxcoef * (log(upper)+1) + lgammaxcoef * digamma(upper) + xcoef < 0){
      xMax <- upper
      maxFound <- TRUE
    }
  }
  
  if(!maxFound){
    xMax <- uniroot(
      f = function(x){
        xlogxcoef * (log(x)+1) + lgammaxcoef * digamma(x) + xcoef
      },
      lower = 1e-9,
      upper = upper,
      extendInt = "no")$root
  }
  
  yMax <- logdNu(x = xMax, xlogxcoef, lgammaxcoef, xcoef, const)
  
  xMax <- xMax/A
  yMax <- yMax - log(A)
  
  # depending on the value of yMax => use constant to decrease the logvalue to 
  # keep the integrand finite
  addConst <- max(c(0,yMax-100))
  
  value <- integrate(f = dANu, lower = 0, upper = upperBound,
                     A = A, xlogxcoef, lgammaxcoef, xcoef,
                     const = const - addConst)$value
  
  if(log){
    value <- log(value) + addConst
  } else {
    value <- exp(log(value) + addConst)
  }
  value
}

##' @export
meanANu <- function(A, xlogxcoef, lgammaxcoef, xcoef, const){
  
  myDeriv <- function(x, K1, K2, K3, const){
    K1 * x * (2*log(x)+1) + K2 * (lgamma(x) + x * digamma(x)) +
      2 * K3 * x + const
  }
  
  # first find where the maximum of the function is
  xMax <- 0.
  
  positiveDerivFound <- FALSE
  lower <- 1e-9
  for(j in 1:22){
    myDerivLower <- myDeriv(x = lower, 
                            K1 = xlogxcoef,
                            K2 = lgammaxcoef, 
                            K3 = xcoef, 
                            const = const)
    if(myDerivLower > 0.){
      positiveDerivFound <- TRUE
      break
    }
    lower <- lower * 10
  }
  if(!positiveDerivFound){
    stop("error in meanANu : no positive first order derivative in interval [1e-9:1e12]")
  }
  
  reverseSignFound <- FALSE
  upper <- max(c(1, lower))
  for (j in 1:12){
    if(myDeriv(x = upper, 
               K1 = xlogxcoef, 
               K2 = lgammaxcoef, 
               K3 = xcoef, 
               const = const) < 0){
      reverseSignFound <- TRUE
      break
    }
    upper <- upper * 10
  }
  if(!reverseSignFound){
    stop("error in meanANu : no maximum found in interval [0:1e12]")
  }
  
  xMax <- uniroot(
    f = function(x){
      myDeriv(x, 
              K1 = xlogxcoef, 
              K2 = lgammaxcoef, 
              K3 = xcoef, 
              const = const)
    },
    lower = lower,
    upper = upper,
    extendInt = "no")$root
  
  
  yMax <- log(xMax) + logdNu(x = xMax, xlogxcoef, lgammaxcoef, xcoef, const)
  
  # depending on the value of yMax => use constant to decrease the logvalue to 
  # keep the integrand finite
  addConst <- max(c(0,yMax-100))
  divConst <- exp(-addConst)
  
  value <- integrate(
    f = function(x){ divConst * x * dNu(x,
                                        xlogxcoef = xlogxcoef,
                                        lgammaxcoef = lgammaxcoef,
                                        xcoef = xcoef,
                                        const = const)},
    lower = 0, 
    upper = Inf)$value
    
  value <- exp(log(value) + addConst)                   
  
  A * value
}

##' @export
meanANuANu <- function(A, xlogxcoef, lgammaxcoef, xcoef, const){
  myDeriv <- function(x, K1, K2, K3, const){
    K1 * x * x * (3*log(x)+1) + K2 * x * (x * lgamma(x) + 2 * digamma(x)) +
      3 * K3 * x * x + 2 * x * const
  }
  
  # first find where the maximum of the function is
  xMax <- 0.
  
  positiveDerivFound <- FALSE
  lower <- 1e-9
  for(j in 1:22){
    myDerivLower <- myDeriv(x = lower, 
                            K1 = xlogxcoef,
                            K2 = lgammaxcoef, 
                            K3 = xcoef, 
                            const = const)
    if(myDerivLower > 0.){
      positiveDerivFound <- TRUE
      break
    }
    lower <- lower * 10
  }
  if(!positiveDerivFound){
    stop("error in meanANu : no positive first order derivative in interval [1e-9:1e12]")
  }
  
  reverseSignFound <- FALSE
  upper <- max(c(1, lower))
  for (j in 1:12){
    if(myDeriv(x = upper, 
               K1 = xlogxcoef, 
               K2 = lgammaxcoef, 
               K3 = xcoef, 
               const = const) < 0){
      reverseSignFound <- TRUE
      break
    }
    upper <- upper * 10
  }
  if(!reverseSignFound){
    stop("error in meanANu : no maximum found in interval [0:1e12]")
  }
  
  xMax <- uniroot(
    f = function(x){
      myDeriv(x, 
              K1 = xlogxcoef, 
              K2 = lgammaxcoef, 
              K3 = xcoef, 
              const = const)
    },
    lower = lower,
    upper = upper,
    extendInt = "no")$root
  
  
  yMax <- 2 * log(xMax) + logdNu(x = xMax, xlogxcoef, lgammaxcoef, xcoef, const)
  
  # depending on the value of yMax => use constant to decrease the logvalue to 
  # keep the integrand finite
  addConst <- max(c(0,yMax-100))
  divConst <- exp(-addConst)
  
  value <- integrate(
    f = function(x){ divConst * x * x * dNu(x,
                                        xlogxcoef = xlogxcoef,
                                        lgammaxcoef = lgammaxcoef,
                                        xcoef = xcoef,
                                        const = const)},
    lower = 0, 
    upper = Inf)$value
  
  value <- exp(log(value) + addConst)                   
  
  A * A *  value
}

##' @export
modeANu <- function(A, xlogxcoef, lgammaxcoef, xcoef, const){
  upper <- 1
  reverseSignFound <- FALSE
  for (j in 1:13){
    if(xlogxcoef * (log(upper)+1) + lgammaxcoef * digamma(upper) + xcoef < 0){
      reverseSignFound <- TRUE
      break
    }
    upper <- upper * 10
  }
  if(!reverseSignFound){
    stop("error in meanANuApprox : no mode found in interval [0:1e12]")
  }
  
  mode <- uniroot(
    f = function(x){
      xlogxcoef * (log(x)+1) + lgammaxcoef * digamma(x) + xcoef
    },
    lower = 1e-9,
    upper = upper,
    extendInt = "no")$root
  
  return(A*mode)
}

##' @export
varianceANuApproxi <- function(A, xlogxcoef, lgammaxcoef, xcoef, const){
  # variance of Laplace approxi
  mode <- modeANu(A, xlogxcoef, lgammaxcoef, xcoef, const)
  var <- -A*A/(xlogxcoef * (1/mode) + lgammaxcoef * trigamma(mode))
}

nudNu <- function(x, xlogxcoef, lgammaxcoef, xcoef, const){
  x * dNu(x, xlogxcoef, lgammaxcoef, xcoef, const)
}

nu2dNu <- function(x, xlogxcoef, lgammaxcoef, xcoef, const){
  x * x * dNu(x, xlogxcoef, lgammaxcoef, xcoef, const)
}

