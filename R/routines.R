
#--------------------------------------------------------------------------------------------------------------------------
#' Power Generalised Weibull (PGW) hazard function.
#' http://rpubs.com/FJRubio/PGW
#--------------------------------------------------------------------------------------------------------------------------
#' @param eta     : scale parameter
#' @param nu      : shape parameter
#' @param delta   : shape parameter
#' @param t       : positive argument
#' @param log: log scale (TRUE or FALSE)
#' @return the value of the PGW hazard function
#' @export

hpgw <- function(t, eta, nu, delta, log = FALSE){
  val <- log(nu) - log(delta) - nu*log(eta) + (nu-1)*log(t) +
    (1/delta - 1)*log( 1 + (t/eta)^nu )
  if(log) return(val) else return(exp(val))
}

#--------------------------------------------------------------------------------------------------------------------------
#' Power Generalised Weibull (PGW) cumulative hazard function.
#' http://rpubs.com/FJRubio/PGW
#--------------------------------------------------------------------------------------------------------------------------
#' @param eta     : scale parameter
#' @param nu      : shape parameter
#' @param delta   : shape parameter
#' @param t       : positive argument
#' @return the value of the PGW cumulative hazard function
#' @export

chpgw <- function(t, eta, nu, delta){
  val <- -1 + ( 1 + (t/eta)^nu )^(1/delta)
  return(val)
}


#--------------------------------------------------------------------------------------------------------------------------
#' Power Generalised Weibull (PGW) quantile function.
#' http://rpubs.com/FJRubio/PGW
#--------------------------------------------------------------------------------------------------------------------------
#' @param eta     : scale parameter
#' @param nu      : shape parameter
#' @param delta   : shape parameter
#' @param t       : positive argument
#' @return the value of the PGW quantile function
#' @export
qpgw <- function(p, eta, nu, delta){
  out <- eta*(  ( 1 - log(1-p) )^delta - 1 )^(1/nu)
  return(out)
}

#----------------------------------------------------------------------------------------
#' Lognormal (LN) Hazard Function.
#----------------------------------------------------------------------------------------
#' @param mu   : log location parameter
#' @param sigma      : scale parameter
#' @param t       : positive argument
#' #' @param log: log scale (TRUE or FALSE)
#' @return the value of the LN hazard function
#' @export

hlnorm <- function(t,mu,sigma, log = FALSE){
  lpdf0 <-  dlnorm(t,mu,sigma, log = T)
  ls0 <- plnorm(t,mu,sigma, lower.tail = FALSE, log.p = T)
  val <- lpdf0 - ls0
  if(log) return(val) else return(exp(val))
}

#----------------------------------------------------------------------------------------
#' Lognormal (LN) Cumulative Hazard Function.
#----------------------------------------------------------------------------------------
#' @param mu   : log location parameter
#' @param sigma      : scale parameter
#' @param t       : positive argument
#' #' @param log: log scale (TRUE or FALSE)
#' @return the value of the LN cumulative hazard function
#' @export
chlnorm <- function(t,mu,sigma){
  H0 <- -plnorm(t,mu,sigma, lower.tail = FALSE, log.p = TRUE)
  return(H0)
}

#----------------------------------------------------------------------------------------
#' Log-logistic (LL) Hazard Function.
#----------------------------------------------------------------------------------------
#' @param mu   : log location parameter
#' @param sigma      : scale parameter
#' @param t       : positive argument
#' #' @param log: log scale (TRUE or FALSE)
#' @return the value of the LL hazard function
#' @export
hllogis <- function(t,mu,sigma, log = FALSE){
  lpdf0 <-  dlogis(log(t),mu,sigma, log = T) - log(t)
  ls0 <- plogis(log(t),mu,sigma, lower.tail = FALSE, log.p = T)
  val <- lpdf0 - ls0
  if(log) return(val) else return(exp(val))
}

#----------------------------------------------------------------------------------------
#' Lognormal (LL) Cumulative Hazard Function.
#----------------------------------------------------------------------------------------
#' @param mu   : log location parameter
#' @param sigma      : scale parameter
#' @param t       : positive argument
#' #' @param log: log scale (TRUE or FALSE)
#' @return the value of the LL cumulative hazard function
#' @export
chllogis <- function(t,mu,sigma){
  H0 <- -plogis(log(t),mu,sigma, lower.tail = FALSE, log.p = TRUE)
  return(H0)
}

#----------------------------------------------------------------------------------------
#' Gamma (G) Hazard Function.
#----------------------------------------------------------------------------------------
#' @param shape   : shape parameter
#' @param scale      : scale parameter
#' @param t       : positive argument
#' #' @param log: log scale (TRUE or FALSE)
#' @return the value of the G hazard function
#' @export
hgamma <- function(t, shape, scale, log = FALSE){
  lpdf0 <-  dgamma(t, shape = shape, scale = scale, log = T)
  ls0 <- pgamma(t, shape = shape, scale = scale, lower.tail = FALSE, log.p = T)
  val <- lpdf0 - ls0
  if(log) return(val) else return(exp(val))
}

#----------------------------------------------------------------------------------------
#' Gamma (G) Cumulative Hazard Function.
#----------------------------------------------------------------------------------------
#' @param shape   : shape parameter
#' @param scale      : scale parameter
#' @param t       : positive argument
#' #' @param log: log scale (TRUE or FALSE)
#' @return the value of the G cumulative hazard function
#' @export
chgamma <- function(t, shape, scale){
  H0 <- -pgamma(t, shape = shape, scale = scale, lower.tail = FALSE, log.p = TRUE)
  return(H0)
}

########################################################################################################
# Conditional log likelihoods (conditional on u, the random effects).
########################################################################################################

#-------------------------------------------------------------------
#' MEGH I.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param u  : random effects
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @return Conditional log-likelihood with Gamma baseline hazard
#' @export

log_likC1_gamma <- function(par, u, times, status, des, des_t){
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  des <- as.matrix(des)
  des_t <- as.matrix(des_t)
  q <- dim(des_t)[2]
  p <- dim(des)[2]
  times.obs <- times[status]
  ae0 <- par[1]; be0 <- par[2];   alpha <- par[3:(2+q)]; beta <- par[(3+q):(2+p+q)]
  x.alpha <- as.vector(des_t%*%alpha)
  x.beta <- as.vector(des%*%beta)
  x.beta.obs <- x.beta[status]
  exp.x.beta.dif <- as.vector(exp( x.beta + u[ID] - x.alpha ))
  exp.x.alpha <- exp(x.alpha)
  exp.x.alpha.obs <- exp.x.alpha[status]
  lhaz0 <- hgamma(times.obs*exp.x.alpha.obs,ae0,be0, log = TRUE) + x.beta.obs + u[ID][status]
  val <- sum(lhaz0) - sum(chgamma(times*exp.x.alpha,ae0,be0)*exp.x.beta.dif)
  return(sum(val))
}

#-------------------------------------------------------------------
#' MEGH I.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param u  : random effects
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @return Conditional log-likelihood with Log-normal baseline hazard
#' @export

log_likC1_LN <- function(par, u, times, status, des, des_t){
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  des <- as.matrix(des)
  des_t <- as.matrix(des_t)
  q <- dim(des_t)[2]
  p <- dim(des)[2]
  times.obs <- times[status]
  ae0 <- par[1]; be0 <- par[2];   alpha <- par[3:(2+q)]; beta <- par[(3+q):(2+p+q)]
  x.alpha <- as.vector(des_t%*%alpha)
  x.beta <- as.vector(des%*%beta)
  x.beta.obs <- x.beta[status]
  exp.x.beta.dif <- as.vector(exp( x.beta + u[ID] - x.alpha ))
  exp.x.alpha <- exp(x.alpha)
  exp.x.alpha.obs <- exp.x.alpha[status]
  lhaz0 <- hlnorm(times.obs*exp.x.alpha.obs,ae0,be0, log = TRUE) + x.beta.obs + u[ID][status]
  val <- sum(lhaz0) - sum(chlnorm(times*exp.x.alpha,ae0,be0)*exp.x.beta.dif)
  return(sum(val))
}

#-------------------------------------------------------------------
#' MEGH I.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param u  : random effects
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @return Conditional log-likelihood with Log-logistic baseline hazard
#' @export
log_likC1_LL <- function(par, u, times, status, des, des_t){
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  des <- as.matrix(des)
  des_t <- as.matrix(des_t)
  q <- dim(des_t)[2]
  p <- dim(des)[2]
  times.obs <- times[status]
  ae0 <- par[1]; be0 <- par[2];   alpha <- par[3:(2+q)]; beta <- par[(3+q):(2+p+q)]
  x.alpha <- as.vector(des_t%*%alpha)
  x.beta <- as.vector(des%*%beta)
  x.beta.obs <- x.beta[status]
  exp.x.beta.dif <- as.vector(exp( x.beta + u[ID] - x.alpha ))
  exp.x.alpha <- exp(x.alpha)
  exp.x.alpha.obs <- exp.x.alpha[status]
  lhaz0 <- hllogis(times.obs*exp.x.alpha.obs,ae0,be0, log = TRUE) + x.beta.obs + u[ID][status]
  val <- sum(lhaz0) - sum(chllogis(times*exp.x.alpha,ae0,be0)*exp.x.beta.dif)
  return(sum(val))
}

#-------------------------------------------------------------------
#' MEGH I.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param u  : random effects
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @return Conditional log-likelihood with PGW baseline hazard
#' @export
log_likC1_PGW <- function(par, u, times, status, des, des_t){
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  des <- as.matrix(des)
  des_t <- as.matrix(des_t)
  q <- dim(des_t)[2]
  p <- dim(des)[2]
  times.obs <- times[status]
  ae0 <- par[1]; be0 <- par[2];  ce0 <- par[3]; alpha <- par[4:(3+q)]; beta <- par[(4+q):(3+p+q)]
  x.alpha <- as.vector(des_t%*%alpha)
  x.beta <- as.vector(des%*%beta)
  x.beta.obs <- x.beta[status]
  exp.x.beta.dif <- as.vector(exp( x.beta + u[ID] - x.alpha ))
  exp.x.alpha <- exp(x.alpha)
  exp.x.alpha.obs <- exp.x.alpha[status]
  lhaz0 <- hpgw(times.obs*exp.x.alpha.obs,ae0,be0,ce0, log = TRUE) + x.beta.obs + u[ID][status]
  val <- sum(lhaz0) - sum(chpgw(times*exp.x.alpha,ae0,be0,ce0)*exp.x.beta.dif)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param u  : random effects
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @return Conditional log-likelihood with Gamma baseline hazard
#' @export

log_likC2_gamma <- function(par, u, times, status, des, des_t){
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  des <- as.matrix(des)
  des_t <- as.matrix(des_t)
  q <- dim(des_t)[2]
  p <- dim(des)[2]
  times.obs <- times[status]
  ae0 <- par[1]; be0 <- par[2];   alpha <- par[3:(2+q)]; beta <- par[(3+q):(2+p+q)]
  x.alpha <- as.vector(des_t%*%alpha)
  x.beta <- as.vector(des%*%beta)
  x.beta.obs <- x.beta[status]
  exp.x.beta.dif <- as.vector(exp( x.beta - x.alpha ))
  exp.x.alpha <- exp(x.alpha + u[ID])
  exp.x.alpha.obs <- exp.x.alpha[status]
  lhaz0 <- hgamma(times.obs*exp.x.alpha.obs,ae0,be0, log = TRUE) + x.beta.obs + u[ID][status]
  val <- sum(lhaz0) - sum(chgamma(times*exp.x.alpha,ae0,be0)*exp.x.beta.dif)
  return(sum(val))
}

#-------------------------------------------------------------------
#' MEGH II.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param u  : random effects
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @return Conditional log-likelihood with Log-Normal baseline hazard
#' @export
log_likC2_LN <- function(par, u, times, status, des, des_t){
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  des <- as.matrix(des)
  des_t <- as.matrix(des_t)
  q <- dim(des_t)[2]
  p <- dim(des)[2]
  times.obs <- times[status]
  ae0 <- par[1]; be0 <- par[2];   alpha <- par[3:(2+q)]; beta <- par[(3+q):(2+p+q)]
  x.alpha <- as.vector(des_t%*%alpha)
  x.beta <- as.vector(des%*%beta)
  x.beta.obs <- x.beta[status]
  exp.x.beta.dif <- as.vector(exp( x.beta - x.alpha ))
  exp.x.alpha <- exp(x.alpha + u[ID])
  exp.x.alpha.obs <- exp.x.alpha[status]
  lhaz0 <- hlnorm(times.obs*exp.x.alpha.obs,ae0,be0, log = TRUE) + x.beta.obs + u[ID][status]
  val <- sum(lhaz0) - sum(chlnorm(times*exp.x.alpha,ae0,be0)*exp.x.beta.dif)
  return(sum(val))
}

#-------------------------------------------------------------------
#' MEGH II.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param u  : random effects
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @return Conditional log-likelihood with Log-logistic baseline hazard
#' @export
log_likC2_LL <- function(par, u, times, status, des, des_t){
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  des <- as.matrix(des)
  des_t <- as.matrix(des_t)
  q <- dim(des_t)[2]
  p <- dim(des)[2]
  times.obs <- times[status]
  ae0 <- par[1]; be0 <- par[2];   alpha <- par[3:(2+q)]; beta <- par[(3+q):(2+p+q)]
  x.alpha <- as.vector(des_t%*%alpha)
  x.beta <- as.vector(des%*%beta)
  x.beta.obs <- x.beta[status]
  exp.x.beta.dif <- as.vector(exp( x.beta - x.alpha ))
  exp.x.alpha <- exp(x.alpha + u[ID])
  exp.x.alpha.obs <- exp.x.alpha[status]
  lhaz0 <- hllogis(times.obs*exp.x.alpha.obs,ae0,be0, log = TRUE) + x.beta.obs + u[ID][status]
  val <- sum(lhaz0) - sum(chllogis(times*exp.x.alpha,ae0,be0)*exp.x.beta.dif)
  return(sum(val))
}

#-------------------------------------------------------------------
#' MEGH I.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param u  : random effects
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @return Conditional log-likelihood with PGW baseline hazard
#' @export
log_likC2_PGW <- function(par, u, times, status, des, des_t){
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  des <- as.matrix(des)
  des_t <- as.matrix(des_t)
  q <- dim(des_t)[2]
  p <- dim(des)[2]
  times.obs <- times[status]
  ae0 <- par[1]; be0 <- par[2];  ce0 <- par[3]; alpha <- par[4:(3+q)]; beta <- par[(4+q):(3+p+q)]
  x.alpha <- as.vector(des_t%*%alpha)
  x.beta <- as.vector(des%*%beta)
  x.beta.obs <- x.beta[status]
  exp.x.beta.dif <- as.vector(exp( x.beta - x.alpha ))
  exp.x.alpha <- exp(x.alpha + u[ID])
  exp.x.alpha.obs <- exp.x.alpha[status]
  lhaz0 <- hpgw(times.obs*exp.x.alpha.obs,ae0,be0,ce0, log = TRUE) + x.beta.obs + u[ID][status]
  val <- sum(lhaz0) - sum(chpgw(times*exp.x.alpha,ae0,be0,ce0)*exp.x.beta.dif)
  return(sum(val))
}


########################################################################################################
# Conditional individual log likelihoods (conditional on u, the random effects)
########################################################################################################
# par: all model parameters in the original parameterisation
# u  : random effects
# des      : design matrix (p x n)
# des_t     : design matrix for time-dependent effects (q x n)
# status : vital status
# times : survival times
# index : individual ID

#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param  par: all model parameters in the original parameterisation
#' @param u  : random effects
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : individual ID
#' @return Conditional individual log-likelihood with Gamma baseline hazard
#' @export

log_likCInd1_gamma <- function(par, u, times, status, des, des_t, index){
  times <- as.vector(times)
  status <- as.vector(as.logical(status))
  des <- as.matrix(des)
  des_t <- as.matrix(des_t)
  indi <- which(ID==index)
  ID.obs <- ID[status]
  indo <- which(ID.obs==index)
  q <- dim(des_t)[2]
  p <- dim(des)[2]
  ae0 <- par[1]; be0 <- par[2];   alpha <- par[3:(2+q)]; beta <- par[(3+q):(2+p+q)]
  x.alpha <- as.vector(des_t%*%alpha)
  x.beta <- as.vector(des%*%beta)
  x.beta.obs <- x.beta[status]
  exp.x.beta.dif <- as.vector(exp( x.beta + u - x.alpha ))
  exp.x.alpha <- exp(x.alpha)
  exp.x.alpha.obs <- exp.x.alpha[status]
  lhaz0 <- hgamma(times.obs*exp.x.alpha.obs,ae0,be0, log = TRUE) + x.beta.obs + u
  val <- sum(lhaz0) - sum(chgamma(times*exp.x.alpha,ae0,be0)*exp.x.beta.dif)
  return(sum(val))
}

#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param  par: all model parameters in the original parameterisation
#' @param u  : random effects
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : individual ID
#' @return Conditional individual log-likelihood with Log-Normal baseline hazard
#' @export
log_likCInd1_LN <- function(par, u, times, status, des, des_t, index){
  indi <- which(ID==index)
  times <- as.vector(times[indi])
  status <- as.vector(as.logical(status[indi]))
  times.obs <- times[status]
  des <- as.matrix(des[indi,])
  des_t <- as.matrix(des_t[indi,])
  q <- dim(des_t)[2]
  p <- dim(des)[2]
  ae0 <- par[1]; be0 <- par[2];   alpha <- par[3:(2+q)]; beta <- par[(3+q):(2+p+q)]
  x.alpha <- as.vector(des_t%*%alpha)
  x.beta <- as.vector(des%*%beta)
  x.beta.obs <- x.beta[status]
  exp.x.beta.dif <- as.vector(exp( x.beta + u - x.alpha ))
  exp.x.alpha <- exp(x.alpha)
  exp.x.alpha.obs <- exp.x.alpha[status]
  lhaz0 <- hlnorm(times.obs*exp.x.alpha.obs,ae0,be0, log = TRUE) + x.beta.obs + u
  val <- sum(lhaz0) - sum(chlnorm(times*exp.x.alpha,ae0,be0)*exp.x.beta.dif)
  return(sum(val))
}

#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param  par: all model parameters in the original parameterisation
#' @param u  : random effects
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : individual ID
#' @return Conditional individual log-likelihood with Log-logistic baseline hazard
#' @export
log_likCInd1_LL <- function(par, u, times, status, des, des_t, index){
  indi <- which(ID==index)
  times <- as.vector(times[indi])
  status <- as.vector(as.logical(status[indi]))
  times.obs <- times[status]
  des <- as.matrix(des[indi,])
  des_t <- as.matrix(des_t[indi,])
  q <- dim(des_t)[2]
  p <- dim(des)[2]
  ae0 <- par[1]; be0 <- par[2];   alpha <- par[3:(2+q)]; beta <- par[(3+q):(2+p+q)]
  x.alpha <- as.vector(des_t%*%alpha)
  x.beta <- as.vector(des%*%beta)
  x.beta.obs <- x.beta[status]
  exp.x.beta.dif <- as.vector(exp( x.beta + u - x.alpha ))
  exp.x.alpha <- exp(x.alpha)
  exp.x.alpha.obs <- exp.x.alpha[status]
  lhaz0 <- hllogis(times.obs*exp.x.alpha.obs,ae0,be0, log = TRUE) + x.beta.obs + u
  val <- sum(lhaz0) - sum(chllogis(times*exp.x.alpha,ae0,be0)*exp.x.beta.dif)
  return(sum(val))
}

#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param  par: all model parameters in the original parameterisation
#' @param u  : random effects
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : individual ID
#' @return Conditional individual log-likelihood with PGW baseline hazard
#' @export
log_likCInd1_PGW <- function(par, u, times, status, des, des_t, index){
  indi <- which(ID==index)
  times <- as.vector(times[indi])
  status <- as.vector(as.logical(status[indi]))
  times.obs <- times[status]
  des <- as.matrix(des[indi,])
  des_t <- as.matrix(des_t[indi,])
  q <- dim(des_t)[2]
  p <- dim(des)[2]
  ae0 <- par[1]; be0 <- par[2];  ce0 <- par[3]; alpha <- par[4:(3+q)]; beta <- par[(4+q):(3+p+q)]
  x.alpha <- as.vector(des_t%*%alpha)
  x.beta <- as.vector(des%*%beta)
  x.beta.obs <- x.beta[status]
  exp.x.beta.dif <- as.vector(exp( x.beta + u - x.alpha ))
  exp.x.alpha <- exp(x.alpha)
  exp.x.alpha.obs <- exp.x.alpha[status]
  lhaz0 <- hpgw(times.obs*exp.x.alpha.obs,ae0,be0,ce0, log = TRUE) + x.beta.obs + u
  val <- sum(lhaz0) - sum(chpgw(times*exp.x.alpha,ae0,be0,ce0)*exp.x.beta.dif)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II
#-------------------------------------------------------------------
#' @param  par: all model parameters in the original parameterisation
#' @param u  : random effects
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : individual ID
#' @return Conditional individual log-likelihood with Gamma baseline hazard
#' @export
log_likCInd2_gamma <- function(par, u, times, status, des, des_t, index){
  indi <- which(ID==index)
  times <- as.vector(times[indi])
  status <- as.vector(as.logical(status[indi]))
  times.obs <- times[status]
  des <- as.matrix(des[indi,])
  des_t <- as.matrix(des_t[indi,])
  q <- dim(des_t)[2]
  p <- dim(des)[2]
  ae0 <- par[1]; be0 <- par[2];   alpha <- par[3:(2+q)]; beta <- par[(3+q):(2+p+q)]
  x.alpha <- as.vector(des_t%*%alpha)
  x.beta <- as.vector(des%*%beta)
  x.beta.obs <- x.beta[status]
  exp.x.beta.dif <- as.vector(exp( x.beta - x.alpha ))
  exp.x.alpha <- exp(x.alpha + u)
  exp.x.alpha.obs <- exp.x.alpha[status]
  lhaz0 <- hgamma(times.obs*exp.x.alpha.obs,ae0,be0, log = TRUE) + x.beta.obs + u
  val <- sum(lhaz0) - sum(chgamma(times*exp.x.alpha,ae0,be0)*exp.x.beta.dif)
  return(sum(val))
}

#-------------------------------------------------------------------
#' MEGH II
#-------------------------------------------------------------------
#' @param  par: all model parameters in the original parameterisation
#' @param u  : random effects
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : individual ID
#' @return Conditional individual log-likelihood with Log-Normal baseline hazard
#' @export
log_likCInd2_LN <- function(par, u, times, status, des, des_t, index){
  indi <- which(ID==index)
  times <- as.vector(times[indi])
  status <- as.vector(as.logical(status[indi]))
  times.obs <- times[status]
  des <- as.matrix(des[indi,])
  des_t <- as.matrix(des_t[indi,])
  q <- dim(des_t)[2]
  p <- dim(des)[2]
  ae0 <- par[1]; be0 <- par[2];   alpha <- par[3:(2+q)]; beta <- par[(3+q):(2+p+q)]
  x.alpha <- as.vector(des_t%*%alpha)
  x.beta <- as.vector(des%*%beta)
  x.beta.obs <- x.beta[status]
  exp.x.beta.dif <- as.vector(exp( x.beta - x.alpha ))
  exp.x.alpha <- exp(x.alpha + u)
  exp.x.alpha.obs <- exp.x.alpha[status]
  lhaz0 <- hlnorm(times.obs*exp.x.alpha.obs,ae0,be0, log = TRUE) + x.beta.obs + u
  val <- sum(lhaz0) - sum(chlnorm(times*exp.x.alpha,ae0,be0)*exp.x.beta.dif)
  return(sum(val))
}

#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param  par: all model parameters in the original parameterisation
#' @param u  : random effects
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : individual ID
#' @return Conditional individual log-likelihood with Log-logistic baseline hazard
#' @export
log_likCInd2_LL <- function(par, u, times, status, des, des_t, index){
  indi <- which(ID==index)
  times <- as.vector(times[indi])
  status <- as.vector(as.logical(status[indi]))
  times.obs <- times[status]
  des <- as.matrix(des[indi,])
  des_t <- as.matrix(des_t[indi,])
  q <- dim(des_t)[2]
  p <- dim(des)[2]
  ae0 <- par[1]; be0 <- par[2];   alpha <- par[3:(2+q)]; beta <- par[(3+q):(2+p+q)]
  x.alpha <- as.vector(des_t%*%alpha)
  x.beta <- as.vector(des%*%beta)
  x.beta.obs <- x.beta[status]
  exp.x.beta.dif <- as.vector(exp( x.beta - x.alpha ))
  exp.x.alpha <- exp(x.alpha + u)
  exp.x.alpha.obs <- exp.x.alpha[status]
  lhaz0 <- hllogis(times.obs*exp.x.alpha.obs,ae0,be0, log = TRUE) + x.beta.obs + u
  val <- sum(lhaz0) - sum(chllogis(times*exp.x.alpha,ae0,be0)*exp.x.beta.dif)
  return(sum(val))
}

#-------------------------------------------------------------------
#' MEGH II
#-------------------------------------------------------------------
#' @param  par: all model parameters in the original parameterisation
#' @param u  : random effects
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : individual ID
#' @return Conditional individual log-likelihood with PGW baseline hazard
#' @export
log_likCInd2_PGW <- function(par, u, times, status, des, des_t, index){
  indi <- which(ID==index)
  times <- as.vector(times[indi])
  status <- as.vector(as.logical(status[indi]))
  times.obs <- times[status]
  des <- as.matrix(des[indi,])
  des_t <- as.matrix(des_t[indi,])
  q <- dim(des_t)[2]
  p <- dim(des)[2]
  ae0 <- par[1]; be0 <- par[2];  ce0 <- par[3]; alpha <- par[4:(3+q)]; beta <- par[(4+q):(3+p+q)]
  x.alpha <- as.vector(des_t%*%alpha)
  x.beta <- as.vector(des%*%beta)
  x.beta.obs <- x.beta[status]
  exp.x.beta.dif <- as.vector(exp( x.beta - x.alpha ))
  exp.x.alpha <- exp(x.alpha + u)
  exp.x.alpha.obs <- exp.x.alpha[status]
  lhaz0 <- hpgw(times.obs*exp.x.alpha.obs,ae0,be0,ce0, log = TRUE) + x.beta.obs + u
  val <- sum(lhaz0) - sum(chpgw(times*exp.x.alpha,ae0,be0,ce0)*exp.x.beta.dif)
  return(sum(val))
}



###############################################################################################
###############################################################################################
###############################################################################################
# FUNCTIONS FOR NORMAL RANDOM EFFECTS
###############################################################################################
###############################################################################################
###############################################################################################



####################################################################################
#' Function to simulate times to event from a model with MEGH structure.
####################################################################################
#' @param seed : seed for simulation
#' @param theta : parameters of the baseline hazard
#' @param beta : regression parameters multiplying the hazard
#' @param n : sample size
#' @param des : Design matrix for hazard level effects
#' @param des_t : Design matrix for time-dependent effects
#' @param ID : Individual identifiers
#' @param Sigma : Standard deviation of the random effects or covariance matrix
#' @param cens : Censoring times
#' @param restr : random effects structure ("I", "II", "GENERAL")
#' @param distr : baseline hazard distribution ("LN", "LL", "PGW" or "gamma")
#' @return a list containing the survival times, vital status, and ID
#' @export
####################################################################################

simMEGH <- function(seed, des = NULL, des_t = NULL, ID, alpha, beta, theta, Sigma, cens, restr, distr){

  if(!is.null(des)){
    des <- as.matrix(des)
    n <- nrow(des)
    des.beta <- as.vector(des%*%beta)
  } else des.beta <- 0

  if(!is.null(des_t)){
    des_t <- as.matrix(des_t)
    n <- nrow(des_t)
    dest.alpha <- as.vector(des_t%*%alpha)
  } else dest.alpha <- 0

  # Simulation

  # Mixed structure MEGH-I
  if(restr == "I"){
    set.seed(seed)
    U <- rnorm(n = n.clust, mean = 0, sd = Sigma) # random effects
    U.ID <- as.vector(U[ID])
    exp.mt  <- as.vector( exp( -dest.alpha ) )
    exp.mdif <- as.vector( exp( dest.alpha - des.beta - U.ID)  )
  }

  # Mixed structure MEGH-II
  if(restr == "II"){
    set.seed(seed)
    U <- rnorm(n = n.clust, mean = 0, sd = Sigma) # random effects
    U.ID <- as.vector(U[ID])
    exp.mt  <- as.vector( exp( -dest.alpha - U.ID) )
    exp.mdif <- as.vector( exp( dest.alpha - des.beta )  )
  }

  # Mixed structure MEGH
  if(restr == "GENERAL"){
    set.seed(seed)
    U <- rmvnorm(n = n.clust, mean = c(0,0), sigma = Sigma) # random effects
    U.IDt <- as.vector(U[ID,1])
    U.ID <- as.vector(U[ID,2])
    exp.mt  <- as.vector( exp( -dest.alpha - U.IDt) )
    exp.mdif <- as.vector( exp( dest.alpha + U.IDt - des.beta  - U.ID)  )
  }

  set.seed(seed)
  ut1  <- runif(n) # uniform variables

  # PGW baseline hazard
  if(distr == "PGW"){
    times <-  qpgw( 1 - exp( log(1-ut1)*exp.mdif ), theta[1], theta[2], theta[3])*exp.mt
  }
  # LN baseline hazard
  if(distr == "LN"){
    times <-  qlnorm( 1 - exp( log(1-ut1)*exp.mdif ), theta[1], theta[2])*exp.mt
  }
  # LL baseline hazard
  if(distr == "LL"){
    times <-  qllogis( 1 - exp( log(1-ut1)*exp.mdif ), theta[1], theta[2])*exp.mt
  }
  # gamma baseline hazard
  if(distr == "gamma"){
    times <-  qgamma( 1 - exp( log(1-ut1)*exp.mdif ), shape = theta[1], scale = theta[2])*exp.mt
  }

  # Corresponding times and vital status
  status <- as.vector(times < cens)
  times <- as.vector(ifelse(status, times, cens))
  status <- as.numeric(status)
  obj <- list(times = as.vector(times), status = status, ID = ID)
  return(obj)
}


########################################################################################################
# Marginal log likelihoods (original parameterisation)
########################################################################################################
# par: all model parameters in the original parameterisation
# des      : design matrix (p x n)
# des_t     : design matrix for time-dependent effects (q x n)
# status : vital status
# times : survival times
# n.clust : number of clusters

#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return  Marginal log likelihood with Gamma baseline hazard (original parameterisation)
#' @export

ml1_gamma <- function(par, times, status, des, des_t, n.clust){
  sigma <- par[1]; param = par[-1];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = i) + dnorm(t,0,sigma, log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = i) )*dnorm(t,0,sigma) )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return  Marginal log likelihood with Log-Normal baseline hazard (original parameterisation)
#' @export
ml1_LN <- function(par, times, status, des, des_t, n.clust){
  sigma <- par[1]; param = par[-1];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_LN(param, u = t, times, status, des, des_t, index = i) + dnorm(t,0,sigma, log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_LN(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = i) )*dnorm(t,0,sigma) )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return  Marginal log likelihood with Log-logistic baseline hazard (original parameterisation)
#' @export
ml1_LL <- function(par, times, status, des, des_t, n.clust){
  sigma <- par[1]; param = par[-1];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_LL(param, u = t, times, status, des, des_t, index = i) + dnorm(t,0,sigma, log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_LL(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = i) )*dnorm(t,0,sigma) )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return  Marginal log likelihood with PGW baseline hazard (original parameterisation)
#' @export
ml1_PGW <- function(par, times, status, des, des_t, n.clust){
  sigma <- par[1]; param = par[-1];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = i) + dnorm(t,0,sigma, log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = i) )*dnorm(t,0,sigma) )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return  Marginal log likelihood with Gamma baseline hazard (original parameterisation)
#' @export
ml2_gamma <- function(par, times, status, des, des_t, n.clust){
  sigma <- par[1]; param = par[-1];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = i) + dnorm(t,0,sigma, log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = i) )*dnorm(t,0,sigma) )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return  Marginal log likelihood with Log-Normal baseline hazard (original parameterisation)
#' @export
ml2_LN <- function(par, times, status, des, des_t, n.clust){
  sigma <- par[1]; param = par[-1];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_LN(param, u = t, times, status, des, des_t, index = i) + dnorm(t,0,sigma, log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_LN(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = i) )*dnorm(t,0,sigma) )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return  Marginal log likelihood with Log-logistic baseline hazard (original parameterisation)
#' @export
ml2_LL <- function(par, times, status, des, des_t, n.clust){
  sigma <- par[1]; param = par[-1];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_LL(param, u = t, times, status, des, des_t, index = i) + dnorm(t,0,sigma, log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_LL(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = i) )*dnorm(t,0,sigma) )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return  Marginal log likelihood with PGW baseline hazard (original parameterisation)
#' @export
ml2_PGW <- function(par, times, status, des, des_t, n.clust){
  sigma <- par[1]; param = par[-1];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = i) + dnorm(t,0,sigma, log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = i) )*dnorm(t,0,sigma) )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


########################################################################################################
# Marginal log likelihoods (reparameterised)
########################################################################################################
# par: all model parameters reparameterised to the real line
# des      : design matrix (p x n)
# des_t     : design matrix for time-dependent effects (q x n)
# status : vital status
# times : survival times

#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return  Marginal log likelihood with Gamma baseline hazard (reparameterised)
#' @export
ml1_gamma_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-1)
  sigma <- exp(par[1]); param[1:2] = exp(par[2:3]); param[-c(1:2)] = par[-c(1:3)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = i) + dnorm(t,0,sigma, log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = i) )*dnorm(t,0,sigma) )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return  Marginal log likelihood with Log-Normal baseline hazard (reparameterised)
#' @export
ml1_LN_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-1)
  sigma <- exp(par[1]); param[1] = par[2]; param[2] = exp(par[3]); param[-c(1:2)] = par[-c(1:3)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_LN(param, u = t, times, status, des, des_t, index = i) + dnorm(t,0,sigma, log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_LN(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = i) )*dnorm(t,0,sigma) )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return  Marginal log likelihood with Log-logistic baseline hazard (reparameterised)
#' @export
ml1_LL_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-1)
  sigma <- exp(par[1]); param[1] = par[2]; param[2] = exp(par[3]); param[-c(1:2)] = par[-c(1:3)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_LL(param, u = t, times, status, des, des_t, index = i) + dnorm(t,0,sigma, log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_LL(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = i) )*dnorm(t,0,sigma) )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return  Marginal log likelihood with PGW baseline hazard (reparameterised)
#' @export
ml1_PGW_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-1)
  sigma <- exp(par[1]); param[1:3] = exp(par[2:4]); param[-c(1:3)] = par[-c(1:4)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = i) + dnorm(t,0,sigma, log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = i) )*dnorm(t,0,sigma) )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return  Marginal log likelihood with Gamma baseline hazard (reparameterised)
#' @export
ml2_gamma_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-1)
  sigma <- exp(par[1]); param[1:2] = exp(par[2:3]); param[-c(1:2)] = par[-c(1:3)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = i) + dnorm(t,0,sigma, log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = i) )*dnorm(t,0,sigma) )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return  Marginal log likelihood with Log-Normal baseline hazard (reparameterised)
#' @export
ml2_LN_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-1)
  sigma <- exp(par[1]); param[1] = par[2]; param[2] = exp(par[3]); param[-c(1:2)] = par[-c(1:3)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_LN(param, u = t, times, status, des, des_t, index = i) + dnorm(t,0,sigma, log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_LN(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = i) )*dnorm(t,0,sigma) )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return  Marginal log likelihood with Log-logistic baseline hazard (reparameterised)
#' @export
ml2_LL_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-1)
  sigma <- exp(par[1]); param[1] = par[2]; param[2] = exp(par[3]); param[-c(1:2)] = par[-c(1:3)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_LL(param, u = t, times, status, des, des_t, index = i) + dnorm(t,0,sigma, log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_LL(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = i) )*dnorm(t,0,sigma) )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return  Marginal log likelihood with PGW baseline hazard (reparameterised)
#' @export
ml2_PGW_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-1)
  sigma <- exp(par[1]); param[1:3] = exp(par[2:4]); param[-c(1:3)] = par[-c(1:4)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = i) + dnorm(t,0,sigma, log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = i) )*dnorm(t,0,sigma) )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


########################################################################################################
# Individual Marginal log likelihoods (original parameterisation)
########################################################################################################
# par: all model parameters in the original parameterisation
# des      : design matrix (p x n)
# des_t     : design matrix for time-dependent effects (q x n)
# status : vital status
# times : survival times
# index : cluster ID

#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return  Individual Marginal log likelihood with Gamma baseline hazard (original parameterisation)
#' @export
ml1Ind_gamma <- function(par, times, status, des, des_t, index){
  sigma <- par[1]; param = par[-1];
  tempf <- Vectorize(function(t) log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = index) + dnorm(t,0,sigma, log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = index) )*dnorm(t,0,sigma) )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return  Individual Marginal log likelihood with Log-Normal baseline hazard (original parameterisation)
#' @export
ml1Ind_LN <- function(par, times, status, des, des_t, index){
  sigma <- par[1]; param = par[-1];
  tempf <- Vectorize(function(t) log_likCInd1_LN(param, u = t, times, status, des, des_t, index = index) + dnorm(t,0,sigma, log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_LN(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = index) )*dnorm(t,0,sigma) )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return  Individual Marginal log likelihood with Log-logistic baseline hazard (original parameterisation)
#' @export
ml1Ind_LL <- function(par, times, status, des, des_t, index){
  sigma <- par[1]; param = par[-1];
  tempf <- Vectorize(function(t) log_likCInd1_LL(param, u = t, times, status, des, des_t, index = index) + dnorm(t,0,sigma, log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_LL(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = index) )*dnorm(t,0,sigma) )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return  Individual Marginal log likelihood with PGW baseline hazard (original parameterisation)
#' @export
ml1Ind_PGW <- function(par, times, status, des, des_t, index){
  sigma <- par[1]; param = par[-1];
  tempf <- Vectorize(function(t) log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = index) + dnorm(t,0,sigma, log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = index) )*dnorm(t,0,sigma) )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return  Individual Marginal log likelihood with Gamma baseline hazard (original parameterisation)
#' @export
ml2Ind_gamma <- function(par, times, status, des, des_t, index){
  sigma <- par[1]; param = par[-1];
  tempf <- Vectorize(function(t) log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = index) + dnorm(t,0,sigma, log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = index) )*dnorm(t,0,sigma) )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return  Individual Marginal log likelihood with Log-Normal baseline hazard (original parameterisation)
#' @export
ml2Ind_LN <- function(par, times, status, des, des_t, index){
  sigma <- par[1]; param = par[-1];
  tempf <- Vectorize(function(t) log_likCInd2_LN(param, u = t, times, status, des, des_t, index = index) + dnorm(t,0,sigma, log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_LN(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = index) )*dnorm(t,0,sigma) )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return  Individual Marginal log likelihood with Log-logistic baseline hazard (original parameterisation)
#' @export
ml2Ind_LL <- function(par, times, status, des, des_t, index){
  sigma <- par[1]; param = par[-1];
  tempf <- Vectorize(function(t) log_likCInd2_LL(param, u = t, times, status, des, des_t, index = index) + dnorm(t,0,sigma, log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_LL(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = index) )*dnorm(t,0,sigma) )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return  Individual Marginal log likelihood with PGW baseline hazard (original parameterisation)
#' @export
ml2Ind_PGW <- function(par, times, status, des, des_t, index){
  sigma <- par[1]; param = par[-1];
  tempf <- Vectorize(function(t) log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = index) + dnorm(t,0,sigma, log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = index) )*dnorm(t,0,sigma) )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}



########################################################################################################
# Individual Marginal log likelihoods (reparameterised)
########################################################################################################
# par: all model parameters reparameterised to the real line
# des      : design matrix (p x n)
# des_t     : design matrix for time-dependent effects (q x n)
# status : vital status
# times : survival times
# index : cluster ID

#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return  Individual Marginal log likelihood with Gamma baseline hazard (reparameterised)
#' @export
ml1Ind_gamma_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-1)
  sigma <- exp(par[1]); param[1:2] = exp(par[2:3]); param[-c(1:2)] = par[-c(1:3)];
  tempf <- Vectorize(function(t) log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = index) + dnorm(t,0,sigma, log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = index) )*dnorm(t,0,sigma) )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return  Individual Marginal log likelihood with Log-Normal baseline hazard (reparameterised)
#' @export
ml1Ind_LN_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-1)
  sigma <- exp(par[1]); param[1] = par[2]; param[2] = exp(par[3]); param[-c(1:2)] = par[-c(1:3)];
  tempf <- Vectorize(function(t) log_likCInd1_LN(param, u = t, times, status, des, des_t, index = index) + dnorm(t,0,sigma, log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_LN(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = index) )*dnorm(t,0,sigma) )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return  Individual Marginal log likelihood with Log-logistic baseline hazard (reparameterised)
#' @export
ml1Ind_LL_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-1)
  sigma <- exp(par[1]); param[1] = par[2]; param[2] = exp(par[3]); param[-c(1:2)] = par[-c(1:3)];
  tempf <- Vectorize(function(t) log_likCInd1_LL(param, u = t, times, status, des, des_t, index = index) + dnorm(t,0,sigma, log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_LL(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = index) )*dnorm(t,0,sigma) )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return  Individual Marginal log likelihood with PGW baseline hazard (reparameterised)
#' @export
ml1Ind_PGW_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-1)
  sigma <- exp(par[1]); param[1:3] = exp(par[2:4]); param[-c(1:3)] = par[-c(1:4)];
  tempf <- Vectorize(function(t) log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = index) + dnorm(t,0,sigma, log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = index) )*dnorm(t,0,sigma) )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return  Individual Marginal log likelihood with Gamma baseline hazard (reparameterised)
#' @export
ml2Ind_gamma_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-1)
  sigma <- exp(par[1]); param[1:2] = exp(par[2:3]); param[-c(1:2)] = par[-c(1:3)];
  tempf <- Vectorize(function(t) log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = index) + dnorm(t,0,sigma, log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = index) )*dnorm(t,0,sigma) )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return  Individual Marginal log likelihood with Log-Normal baseline hazard (reparameterised)
#' @export
ml2Ind_LN_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-1)
  sigma <- exp(par[1]); param[1] = par[2]; param[2] = exp(par[3]); param[-c(1:2)] = par[-c(1:3)];
  tempf <- Vectorize(function(t) log_likCInd2_LN(param, u = t, times, status, des, des_t, index = index) + dnorm(t,0,sigma, log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_LN(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = index) )*dnorm(t,0,sigma) )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return  Individual Marginal log likelihood with Log-logistic baseline hazard (reparameterised)
#' @export
ml2Ind_LL_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-1)
  sigma <- exp(par[1]); param[1] = par[2]; param[2] = exp(par[3]); param[-c(1:2)] = par[-c(1:3)];
  tempf <- Vectorize(function(t) log_likCInd2_LL(param, u = t, times, status, des, des_t, index = index) + dnorm(t,0,sigma, log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_LL(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = index) )*dnorm(t,0,sigma) )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return  Individual Marginal log likelihood with PGW baseline hazard (reparameterised)
#' @export
ml2Ind_PGW_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-1)
  sigma <- exp(par[1]); param[1:3] = exp(par[2:4]); param[-c(1:3)] = par[-c(1:4)];
  tempf <- Vectorize(function(t) log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = index) + dnorm(t,0,sigma, log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = index) )*dnorm(t,0,sigma) )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#########################################################################################################
# Marginal survival function for cluster "index" at time t
#########################################################################################################
# t: time
# par: Marginal maximum likelihood estimator (original parameterisation)
# des         : Design matrix
# des_t        : Design matrix for time-dependent effects
# index : individual ID
# NMC : number of Monte Carlo samples to be used to produce the marginal

#-------------------------------------------------------------------
#' MEGH I.
#-------------------------------------------------------------------
#' Marginal survival function for cluster "index" at time t
#-------------------------------------------------------------------
#' @param  t: time
#' @param par: Marginal maximum likelihood estimator (original parameterisation)
#' @param des  : Design matrix
#' @param des_t : Design matrix for time-dependent effects
#' @param index : individual ID
#' @param NMC : number of Monte Carlo samples to be used to produce the marginal
#' @return The marginal survival function for cluster "index" at time t with PGW baseline hazard
#' @export
MSurvPGW1 <- function(t, par, des, des_t, index, NMC){
  sigma <- par[1]
  theta <- par[2:4]
  indi <- which(ID==index)
  des <- as.matrix(des); des_t <- as.matrix(des_t)
  des <- as.matrix(des[indi,])
  des_t <- as.matrix(des_t[indi,])
  n.ind <- length(indi)
  q <- dim(des_t)[2]; p <- dim(des)[2]
  alpha <- par[5:(4+q)]; beta <- par[(5+q):(4+p+q)]
  # simulated random effects for Monte Carlo integration
  u <- rnorm(NMC, mean = 0, sd = sigma)
  # Marginal survival
  sfit.PGW <- vector()
  for(j in 1:NMC){
    exp.mlePGW <- as.vector(exp(des_t%*%alpha))
    exp.mlePGW.dif <- as.vector(exp(des%*%beta - des_t%*%alpha + u[j]))
    sfit.PGW[j] <- mean(exp(-chpgw(t*exp.mlePGW,theta[1],theta[2],theta[3])*exp.mlePGW.dif ))
  }

  return(mean(sfit.PGW))

}

#-------------------------------------------------------------------
#' MEGH I.
#-------------------------------------------------------------------
#' Marginal survival function for cluster "index" at time t
#-------------------------------------------------------------------
#' @param  t: time
#' @param par: Marginal maximum likelihood estimator (original parameterisation)
#' @param des  : Design matrix
#' @param des_t : Design matrix for time-dependent effects
#' @param index : individual ID
#' @param NMC : number of Monte Carlo samples to be used to produce the marginal
#' @return The marginal survival function for cluster "index" at time t with Log-logistic baseline hazard
#' @export
MSurvLL1 <- function(t, par, des, des_t, index, NMC){
  sigma <- par[1]
  theta <- par[2:3]
  indi <- which(ID==index)
  des <- as.matrix(des); des_t <- as.matrix(des_t)
  des <- as.matrix(des[indi,])
  des_t <- as.matrix(des_t[indi,])
  n.ind <- length(indi)
  q <- dim(des_t)[2]; p <- dim(des)[2]
  alpha <- par[4:(3+q)]; beta <- par[(4+q):(3+p+q)]
  # simulated random effects for Monte Carlo integration
  u <- rnorm(NMC, mean = 0, sd = sigma)
  # Marginal survival
  sfit.LL <- vector()
  for(j in 1:NMC){
    exp.mleLL <- as.vector(exp(des_t%*%alpha))
    exp.mleLL.dif <- as.vector(exp(des%*%beta - des_t%*%alpha + u[j]))
    sfit.LL[j] <- mean(exp(-chllogis(t*exp.mleLL,theta[1],theta[2])*exp.mleLL.dif ))
  }

  return(mean(sfit.LL))

}


#-------------------------------------------------------------------
#' MEGH II.
#-------------------------------------------------------------------
#' Marginal survival function for cluster "index" at time t
#-------------------------------------------------------------------
#' @param  t: time
#' @param par: Marginal maximum likelihood estimator (original parameterisation)
#' @param des  : Design matrix
#' @param des_t : Design matrix for time-dependent effects
#' @param index : individual ID
#' @param NMC : number of Monte Carlo samples to be used to produce the marginal
#' @return The marginal survival function for cluster "index" at time t with PGW baseline hazard
#' @export
MSurvPGW2 <- function(t, par, des, des_t, index, NMC){
  sigma <- par[1]
  theta <- par[2:4]
  indi <- which(ID==index)
  des <- as.matrix(des); des_t <- as.matrix(des_t)
  des <- as.matrix(des[indi,])
  des_t <- as.matrix(des_t[indi,])
  n.ind <- length(indi)
  q <- dim(des_t)[2]; p <- dim(des)[2]
  alpha <- par[5:(4+q)]; beta <- par[(5+q):(4+p+q)]
  # simulated random effects for Monte Carlo integration
  u <- rnorm(NMC, mean = 0, sd = sigma)
  # Marginal survival
  sfit.PGW <- vector()
  for(j in 1:NMC){
    exp.mlePGW <- as.vector(exp(des_t%*%alpha + u[j]))
    exp.mlePGW.dif <- as.vector(exp(des%*%beta - des_t%*%alpha))
    sfit.PGW[j] <- mean(exp(-chpgw(t*exp.mlePGW,theta[1],theta[2],theta[3])*exp.mlePGW.dif ))
  }

  return(mean(sfit.PGW))

}

#-------------------------------------------------------------------
#' MEGH II.
#-------------------------------------------------------------------
#' Marginal survival function for cluster "index" at time t
#-------------------------------------------------------------------
#' @param  t: time
#' @param par: Marginal maximum likelihood estimator (original parameterisation)
#' @param des  : Design matrix
#' @param des_t : Design matrix for time-dependent effects
#' @param index : individual ID
#' @param NMC : number of Monte Carlo samples to be used to produce the marginal
#' @return The marginal survival function for cluster "index" at time t with Log-logistic baseline hazard
#' @export
MSurvLL2 <- function(t, par, des, des_t, index, NMC){
  sigma <- par[1]
  theta <- par[2:3]
  indi <- which(ID==index)
  des <- as.matrix(des); des_t <- as.matrix(des_t)
  des <- as.matrix(des[indi,])
  des_t <- as.matrix(des_t[indi,])
  n.ind <- length(indi)
  q <- dim(des_t)[2]; p <- dim(des)[2]
  alpha <- par[4:(3+q)]; beta <- par[(4+q):(3+p+q)]
  # simulated random effects for Monte Carlo integration
  u <- rnorm(NMC, mean = 0, sd = sigma)
  # Marginal survival
  sfit.LL <- vector()
  for(j in 1:NMC){
    exp.mleLL <- as.vector(exp(des_t%*%alpha + u[j]))
    exp.mleLL.dif <- as.vector(exp(des%*%beta - des_t%*%alpha))
    sfit.LL[j] <- mean(exp(-chllogis(t*exp.mleLL,theta[1],theta[2])*exp.mleLL.dif ))
  }

  return(mean(sfit.LL))

}




###############################################################################################
###############################################################################################
###############################################################################################
# FUNCTIONS FOR STUDENT-T RANDOM EFFECTS
###############################################################################################
###############################################################################################
###############################################################################################


####################################################################################
# Function to simulate times to event from a model with MEGH structure
####################################################################################
# seed : seed for simulation
# theta: parameters of the baseline hazard
# beta : regression parameters multiplying the hazard
# n : sample size
# des : Design matrix
# des_t : Design matrix for time-dependent effects
# ID : Individual identifiers
# sigma : Standard deviation of the random effects
# nu : degrees of freedom
# cens : Censoring times
# restr : random effects structure ("I", "II" or "III")
# distr : baseline hazard distribution ("LN", "LL", "PGW" or "gamma")
####################################################################################


# FUNCTIONS FOR STUDENT-T RANDOM EFFECTS
###############################################################################################
###############################################################################################
###############################################################################################


############################################################################################################
#' Function to simulate times to event from a model with MEGH structure and Student-t random effects
############################################################################################################
#' @param seed : seed for simulation
#' @param theta: parameters of the baseline hazard
#' @param beta : regression parameters multiplying the hazard
#' @param n : sample size
#' @param des : Design matrix
#' @param des_t : Design matrix for time-dependent effects
#' @param ID : Individual identifiers
#' @param sigma : Standard deviation of the random effects
#' @param nu : degrees of freedom
#' @param cens : Censoring times
#' @param restr : random effects structure ("I", "II" or "III")
#' @param distr : baseline hazard distribution ("LN", "LL", "PGW" or "gamma")
#' @return a list containing the survival times, vital status, and ID
#' @export
####################################################################################

simMEGH_t <- function(seed, des = NULL, des_t = NULL, ID, alpha, beta, theta, sigma, nu, cens, restr, distr){

  if(!is.null(des)){
    des <- as.matrix(des)
    n <- nrow(des)
    des.beta <- as.vector(des%*%beta)
  } else des.beta <- 0

  if(!is.null(des)){
    des_t <- as.matrix(des_t)
    n <- nrow(des_t)
    dest.alpha <- as.vector(des_t%*%alpha)
  }   else dest.alpha <- 0

  # Simulation
  set.seed(seed)
  ut1  <- runif(n) # uniform variables
  U <- sigma*rt(n = n.clust, df = nu) # random effects
  U.ID <- as.vector(U[ID])

  # Mixed structure I
  if(restr == "I"){
    exp.mt  <- as.vector( exp( -dest.alpha ) )
    exp.mdif <- as.vector( exp( dest.alpha - des.beta - U.ID)  )
  }

  # Mixed structure II
  if(restr == "II"){
    exp.mt  <- as.vector( exp( -dest.alpha - U.ID) )
    exp.mdif <- as.vector( exp( dest.alpha - des.beta )  )
  }

  # PGW baseline hazard
  if(distr == "PGW"){
    times <-  qpgw( 1 - exp( log(1-ut1)*exp.mdif ), theta[1], theta[2], theta[3])*exp.mt
  }
  # LN baseline hazard
  if(distr == "LN"){
    times <-  qlnorm( 1 - exp( log(1-ut1)*exp.mdif ), theta[1], theta[2])*exp.mt
  }
  # LL baseline hazard
  if(distr == "LL"){
    times <-  qllogis( 1 - exp( log(1-ut1)*exp.mdif ), theta[1], theta[2])*exp.mt
  }
  # gamma baseline hazard
  if(distr == "gamma"){
    times <-  qgamma( 1 - exp( log(1-ut1)*exp.mdif ), shape = theta[1], scale = theta[2])*exp.mt
  }

  # Corresponding times and vital status
  status <- as.vector(times < cens)
  times <- as.vector(ifelse(status, times, cens))
  status <- as.numeric(status)
  obj <- list(times = as.vector(times), status = status, ID = ID)
  return(obj)
}


########################################################################################################
# Marginal log likelihoods (original parameterisation)
########################################################################################################
# par: all model parameters in the original parameterisation
# des      : design matrix (p x n)
# des_t     : design matrix for time-dependent effects (q x n)
# status : vital status
# times : survival times
# n.clust : number of clusters


#-------------------------------------------------------------------
#' MEGH I. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with Gamma baseline hazard (original parameterisation)
#' @export

tml1_gamma <- function(par, times, status, des, des_t, n.clust){
  sigma <- par[1]; nu = par[2]; param = par[-c(1:2)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = i) +
                         dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = i) )*dt(t/sigma, df = nu)/sigma )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with Log-Normal baseline hazard (original parameterisation)
#' @export
tml1_LN <- function(par, times, status, des, des_t, n.clust){
  sigma <- par[1]; nu = par[2]; param = par[-c(1:2)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_LN(param, u = t, times, status, des, des_t, index = i) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_LN(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = i) )*dt(t/sigma, df = nu)/sigma )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with Log-logistic baseline hazard (original parameterisation)
#' @export
tml1_LL <- function(par, times, status, des, des_t, n.clust){
  sigma <- par[1]; nu = par[2]; param = par[-c(1:2)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_LL(param, u = t, times, status, des, des_t, index = i) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_LL(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = i) )*dt(t/sigma, df = nu)/sigma )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with PGW baseline hazard (original parameterisation)
#' @export
tml1_PGW <- function(par, times, status, des, des_t, n.clust){
  sigma <- par[1]; nu = par[2]; param = par[-c(1:2)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = i) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = i) )*dt(t/sigma, df = nu)/sigma )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with Gamma baseline hazard (original parameterisation)
#' @export
tml2_gamma <- function(par, times, status, des, des_t, n.clust){
  sigma <- par[1]; nu = par[2]; param = par[-c(1:2)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = i) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = i) )*dt(t/sigma, df = nu)/sigma )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with Log-Normal baseline hazard (original parameterisation)
#' @export
tml2_LN <- function(par, times, status, des, des_t, n.clust){
  sigma <- par[1]; nu = par[2]; param = par[-c(1:2)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_LN(param, u = t, times, status, des, des_t, index = i) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_LN(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = i) )*dt(t/sigma, df = nu)/sigma )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with Log-logistic baseline hazard (original parameterisation)
#' @export
tml2_LL <- function(par, times, status, des, des_t, n.clust){
  sigma <- par[1]; nu = par[2]; param = par[-c(1:2)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_LL(param, u = t, times, status, des, des_t, index = i) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_LL(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = i) )*dt(t/sigma, df = nu)/sigma )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with PGW baseline hazard (original parameterisation)
#' @export
tml2_PGW <- function(par, times, status, des, des_t, n.clust){
  sigma <- par[1]; nu = par[2]; param = par[-c(1:2)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = i) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = i) )*dt(t/sigma, df = nu)/sigma )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


########################################################################################################
# Marginal log likelihoods (reparameterised)
########################################################################################################
# par: all model parameters reparameterised to the real line
# des      : design matrix (p x n)
# des_t     : design matrix for time-dependent effects (q x n)
# status : vital status
# times : survival times

#-------------------------------------------------------------------
#' MEGH I. Student-t random effects.
#-------------------------------------------------------------------
#' @param  par: all model parameters reparameterised to the real line
#' @param des : design matrix (p x n)
#' @param des_t  : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @return Marginal log likelihood with Gamma baseline hazard (reparameterised)
#' @export

tml1_gamma_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-2)
  sigma <- exp(par[1]); nu = exp(par[2]); param[1:2] = exp(par[3:4]); param[-c(1:2)] = par[-c(1:4)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = i) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = i) )*dt(t/sigma, df = nu)/sigma )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. Student-t random effects.
#-------------------------------------------------------------------
#' @param  par: all model parameters reparameterised to the real line
#' @param des : design matrix (p x n)
#' @param des_t  : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @return Marginal log likelihood with Log-Normal baseline hazard (reparameterised)
#' @export
tml1_LN_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-2)
  sigma <- exp(par[1]); nu = exp(par[2]); param[1] = par[3]; param[2] = exp(par[4]); param[-c(1:2)] = par[-c(1:4)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_LN(param, u = t, times, status, des, des_t, index = i) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_LN(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = i) )*dt(t/sigma, df = nu)/sigma )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. Student-t random effects.
#-------------------------------------------------------------------
#' @param  par: all model parameters reparameterised to the real line
#' @param des : design matrix (p x n)
#' @param des_t  : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @return Marginal log likelihood with Log-logistic baseline hazard (reparameterised)
#' @export
tml1_LL_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-2)
  sigma <- exp(par[1]); nu = exp(par[2]); param[1] = par[3]; param[2] = exp(par[4]); param[-c(1:2)] = par[-c(1:4)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_LL(param, u = t, times, status, des, des_t, index = i) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_LL(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = i) )*dt(t/sigma, df = nu)/sigma )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. Student-t random effects.
#-------------------------------------------------------------------
#' @param  par: all model parameters reparameterised to the real line
#' @param des : design matrix (p x n)
#' @param des_t  : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @return Marginal log likelihood with PGW baseline hazard (reparameterised)
#' @export
tml1_PGW_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-2)
  sigma <- exp(par[1]); nu = exp(par[2]); param[1:3] = exp(par[3:5]); param[-c(1:3)] = par[-c(1:5)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = i) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = i) )*dt(t/sigma, df = nu)/sigma )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. Student-t random effects.
#-------------------------------------------------------------------
#' @param  par: all model parameters reparameterised to the real line
#' @param des : design matrix (p x n)
#' @param des_t  : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @return Marginal log likelihood with Gamma baseline hazard (reparameterised)
#' @export
tml2_gamma_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-2)
  sigma <- exp(par[1]); nu = exp(par[2]); param[1:2] = exp(par[3:4]); param[-c(1:2)] = par[-c(1:4)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = i) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = i) )*dt(t/sigma, df = nu)/sigma )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. Student-t random effects.
#-------------------------------------------------------------------
#' @param  par: all model parameters reparameterised to the real line
#' @param des : design matrix (p x n)
#' @param des_t  : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @return Marginal log likelihood with Log-Normal baseline hazard (reparameterised)
#' @export
tml2_LN_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-2)
  sigma <- exp(par[1]); nu = exp(par[2]); param[1] = par[3]; param[2] = exp(par[4]); param[-c(1:2)] = par[-c(1:4)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_LN(param, u = t, times, status, des, des_t, index = i) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_LN(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = i) )*dt(t/sigma, df = nu)/sigma )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. Student-t random effects.
#-------------------------------------------------------------------
#' @param  par: all model parameters reparameterised to the real line
#' @param des : design matrix (p x n)
#' @param des_t  : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @return Marginal log likelihood with Log-logistic baseline hazard (reparameterised)
#' @export
tml2_LL_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-2)
  sigma <- exp(par[1]); nu = exp(par[2]); param[1] = par[3]; param[2] = exp(par[4]); param[-c(1:2)] = par[-c(1:4)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_LL(param, u = t, times, status, des, des_t, index = i) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_LL(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = i) )*dt(t/sigma, df = nu)/sigma )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. Student-t random effects.
#-------------------------------------------------------------------
#' @param  par: all model parameters reparameterised to the real line
#' @param des : design matrix (p x n)
#' @param des_t  : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @return Marginal log likelihood with PGW baseline hazard (reparameterised)
#' @export
tml2_PGW_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-2)
  sigma <- exp(par[1]); nu = exp(par[2]); param[1:3] = exp(par[3:5]); param[-c(1:3)] = par[-c(1:5)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = i) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = i) )*dt(t/sigma, df = nu)/sigma )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


########################################################################################################
# Individual Marginal log likelihoods (original parameterisation)
########################################################################################################
# par: all model parameters in the original parameterisation
# des      : design matrix (p x n)
# des_t     : design matrix for time-dependent effects (q x n)
# status : vital status
# times : survival times
# index : cluster ID


#-------------------------------------------------------------------
#' MEGH I. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log likelihood with Gamma baseline hazard (original parameterisation)
#' @export

tml1Ind_gamma <- function(par, times, status, des, des_t, index){
  sigma <- par[1]; nu = par[2]; param = par[-c(1:2)];
  tempf <- Vectorize(function(t) log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = index) +
                       dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = index) )*dt(t/sigma, df = nu)/sigma )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log likelihood with Log-Normal baseline hazard (original parameterisation)
#' @export
tml1Ind_LN <- function(par, times, status, des, des_t, index){
  sigma <- par[1]; nu = par[2]; param = par[-c(1:2)];
  tempf <- Vectorize(function(t) log_likCInd1_LN(param, u = t, times, status, des, des_t, index = index) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_LN(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = index) )*dt(t/sigma, df = nu)/sigma )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log likelihood with Log-logistic baseline hazard (original parameterisation)
#' @export
tml1Ind_LL <- function(par, times, status, des, des_t, index){
  sigma <- par[1]; nu = par[2]; param = par[-c(1:2)];
  tempf <- Vectorize(function(t) log_likCInd1_LL(param, u = t, times, status, des, des_t, index = index) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_LL(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = index) )*dt(t/sigma, df = nu)/sigma )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log likelihood with PGW baseline hazard (original parameterisation)
#' @export
tml1Ind_PGW <- function(par, times, status, des, des_t, index){
  sigma <- par[1]; nu = par[2]; param = par[-c(1:2)];
  tempf <- Vectorize(function(t) log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = index) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = index) )*dt(t/sigma, df = nu)/sigma )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log likelihood with Gamma baseline hazard (original parameterisation)
#' @export
tml2Ind_gamma <- function(par, times, status, des, des_t, index){
  sigma <- par[1]; nu = par[2]; param = par[-c(1:2)];
  tempf <- Vectorize(function(t) log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = index) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = index) )*dt(t/sigma, df = nu)/sigma )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log likelihood with Log-Normal baseline hazard (original parameterisation)
#' @export
tml2Ind_LN <- function(par, times, status, des, des_t, index){
  sigma <- par[1]; nu = par[2]; param = par[-c(1:2)];
  tempf <- Vectorize(function(t) log_likCInd2_LN(param, u = t, times, status, des, des_t, index = index) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_LN(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = index) )*dt(t/sigma, df = nu)/sigma )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log likelihood with Log-logistic baseline hazard (original parameterisation)
#' @export
tml2Ind_LL <- function(par, times, status, des, des_t, index){
  sigma <- par[1]; nu = par[2]; param = par[-c(1:2)];
  tempf <- Vectorize(function(t) log_likCInd2_LL(param, u = t, times, status, des, des_t, index = index) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_LL(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = index) )*dt(t/sigma, df = nu)/sigma )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log likelihood with PGW baseline hazard (original parameterisation)
#' @export
tml2Ind_PGW <- function(par, times, status, des, des_t, index){
  sigma <- par[1]; nu = par[2]; param = par[-c(1:2)];
  tempf <- Vectorize(function(t) log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = index) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = index) )*dt(t/sigma, df = nu)/sigma )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


########################################################################################################
# Individual Marginal log likelihoods (reparameterised)
########################################################################################################
# par: all model parameters reparameterised to the real line
# des      : design matrix (p x n)
# des_t     : design matrix for time-dependent effects (q x n)
# status : vital status
# times : survival times
# index : cluster ID

#-------------------------------------------------------------------
#' MEGH I. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log likelihood with Gamma baseline hazard (reparameterised)
#' @export
tml1Ind_gamma_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-2)
  sigma <- exp(par[1]); nu = exp(par[2]); param[1:2] = exp(par[3:4]); param[-c(1:2)] = par[-c(1:4)];
  tempf <- Vectorize(function(t) log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = index) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = index) )*dt(t/sigma, df = nu)/sigma )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log likelihood with Log-Normal baseline hazard (reparameterised)
#' @export
tml1Ind_LN_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-2)
  sigma <- exp(par[1]); nu = exp(par[2]); param[1] = par[3]; param[2] = exp(par[4]); param[-c(1:2)] = par[-c(1:4)];
  tempf <- Vectorize(function(t) log_likCInd1_LN(param, u = t, times, status, des, des_t, index = index) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_LN(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = index) )*dt(t/sigma, df = nu)/sigma )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log likelihood with Log-logistic baseline hazard (reparameterised)
#' @export
tml1Ind_LL_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-2)
  sigma <- exp(par[1]); nu = exp(par[2]); param[1] = par[3]; param[2] = exp(par[4]); param[-c(1:2)] = par[-c(1:4)];
  tempf <- Vectorize(function(t) log_likCInd1_LL(param, u = t, times, status, des, des_t, index = index) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_LL(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = index) )*dt(t/sigma, df = nu)/sigma )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log likelihood with PGW baseline hazard (reparameterised)
#' @export
tml1Ind_PGW_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-2)
  sigma <- exp(par[1]); nu = exp(par[2]); param[1:3] = exp(par[3:5]); param[-c(1:3)] = par[-c(1:5)];
  tempf <- Vectorize(function(t) log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = index) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = index) )*dt(t/sigma, df = nu)/sigma )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log likelihood with Gamma baseline hazard (reparameterised)
#' @export
tml2Ind_gamma_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-2)
  sigma <- exp(par[1]); nu = exp(par[2]); param[1:2] = exp(par[3:4]); param[-c(1:2)] = par[-c(1:4)];
  tempf <- Vectorize(function(t) log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = index) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = index) )*dt(t/sigma, df = nu)/sigma )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log likelihood with Log-Normal baseline hazard (reparameterised)
#' @export
tml2Ind_LN_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-2)
  sigma <- exp(par[1]); nu = exp(par[2]); param[1] = par[3]; param[2] = exp(par[4]); param[-c(1:2)] = par[-c(1:4)];
  tempf <- Vectorize(function(t) log_likCInd2_LN(param, u = t, times, status, des, des_t, index = index) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_LN(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = index) )*dt(t/sigma, df = nu)/sigma )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log likelihood with Log-Logistic baseline hazard (reparameterised)
#' @export
tml2Ind_LL_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-2)
  sigma <- exp(par[1]); nu = exp(par[2]); param[1] = par[3]; param[2] = exp(par[4]); param[-c(1:2)] = par[-c(1:4)];
  tempf <- Vectorize(function(t) log_likCInd2_LL(param, u = t, times, status, des, des_t, index = index) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_LL(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = index) )*dt(t/sigma, df = nu)/sigma )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. Student-t random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log likelihood with PGW baseline hazard (reparameterised)
#' @export
tml2Ind_PGW_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-2)
  sigma <- exp(par[1]); nu = exp(par[2]); param[1:3] = exp(par[3:5]); param[-c(1:3)] = par[-c(1:5)];
  tempf <- Vectorize(function(t) log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = index) + dt( t/sigma, df =  nu, log = TRUE) - log(sigma) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = index) )*dt(t/sigma, df = nu)/sigma )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}

#########################################################################################################
# Marginal survival function for cluster "index" at time t
#########################################################################################################
# t: time
# par: Marginal maximum likelihood estimator (original parameterisation)
# des         : Design matrix
# des_t        : Design matrix for time-dependent effects
# index : individual ID
# NMC : number of Monte Carlo samples to be used to produce the marginal

#-------------------------------------------------------------------
#' MEGH I. Student-t random effects
#-------------------------------------------------------------------
#' @param t: time
#' @param par: Marginal maximum likelihood estimator (original parameterisation)
#' @param des : Design matrix
#' @param des_t : Design matrix for time-dependent effects
#' @param index : individual ID
#' @param NMC : number of Monte Carlo samples to be used to produce the marginal
#' @return The marginal survival function for cluster "index" at time t with PGW baseline hazard
#' @export
MSurvPGW1t <- function(t, par, des, des_t, index, NMC){
  sigma <- par[1]
  nu <- par[2]
  theta <- par[4:6]
  indi <- which(ID==index)
  des <- as.matrix(des); des_t <- as.matrix(des_t)
  des <- as.matrix(des[indi,])
  des_t <- as.matrix(des_t[indi,])
  n.ind <- length(indi)
  q <- dim(des_t)[2]; p <- dim(des)[2]
  alpha <- par[5:(4+q)]; beta <- par[(5+q):(4+p+q)]
  # simulated random effects for Monte Carlo integration
  u <- sigma*rt(NMC, df = nu)
  # Marginal survival
  sfit.PGW <- vector()
  for(j in 1:NMC){
    exp.mlePGW <- as.vector(exp(des_t%*%alpha))
    exp.mlePGW.dif <- as.vector(exp(des%*%beta - des_t%*%alpha + u[j]))
    sfit.PGW[j] <- mean(exp(-chpgw(t*exp.mlePGW,MLE[1],MLE[2],MLE[3])*exp.mlePGW.dif ))
  }

  return(mean(sfit.PGW))

}



#-------------------------------------------------------------------
#' MEGH II. Student-t random effects
#-------------------------------------------------------------------
#' @param t: time
#' @param par: Marginal maximum likelihood estimator (original parameterisation)
#' @param des : Design matrix
#' @param des_t : Design matrix for time-dependent effects
#' @param index : individual ID
#' @param NMC : number of Monte Carlo samples to be used to produce the marginal
#' @return The marginal survival function for cluster "index" at time t with PGW baseline hazard
#' @export

MSurvPGW2t <- function(t, par, des, des_t, index, NMC){
  sigma <- par[1]
  nu <- par[2]
  theta <- par[4:6]
  indi <- which(ID==index)
  des <- as.matrix(des); des_t <- as.matrix(des_t)
  des <- as.matrix(des[indi,])
  des_t <- as.matrix(des_t[indi,])
  n.ind <- length(indi)
  q <- dim(des_t)[2]; p <- dim(des)[2]
  alpha <- par[5:(4+q)]; beta <- par[(5+q):(4+p+q)]
  # simulated random effects for Monte Carlo integration
  u <- sigma*rt(NMC, df = nu)
  # Marginal survival
  sfit.PGW <- vector()
  for(j in 1:NMC){
    exp.mlePGW <- as.vector(exp(des_t%*%alpha + u[j]))
    exp.mlePGW.dif <- as.vector(exp(des%*%beta - des_t%*%alpha))
    sfit.PGW[j] <- mean(exp(-chpgw(t*exp.mlePGW,MLE[1],MLE[2],MLE[3])*exp.mlePGW.dif ))
  }

  return(mean(sfit.PGW))

}







###############################################################################################
###############################################################################################
###############################################################################################
# FUNCTIONS FOR CENTRED TWO PIECE NORMAL RANDOM EFFECTS
###############################################################################################
###############################################################################################
###############################################################################################



############################################################################################################
#' Function to simulate times to event from a model with MEGH structure and twopiece normal random effects
############################################################################################################
#' @param seed  : seed for simulation
#' @param theta : parameters of the baseline hazard
#' @param beta : regression parameters multiplying the hazard
#' @param n : sample size
#' @param des : Design matrix for hazard-level effects
#' @param des_t : Design matrix for time-dependent effects
#' @param ID : Individual identifiers
#' @param skew : skewness parameter (eps parameterisation) of the random effects
#' @param sigma : Standard deviation of the random effects
#' @param cens : Censoring times
#' @param restr : random effects structure ("I", "II" or "III")
#' @param distr : baseline hazard distribution ("LN", "LL", "PGW" or "gamma")
#' @return list of survival times, vital status, and ID.
#' @export
####################################################################################

simMEGH_TPN <- function(seed, des = NULL, des_t = NULL, ID, alpha, beta, theta, sigma, skew, cens, restr, distr){

  if(!is.null(des)){
    des <- as.matrix(des)
    n <- nrow(des)
    des.beta <- as.vector(des%*%beta)
  } else des.beta <- 0

  if(!is.null(des)){
    des_t <- as.matrix(des_t)
    n <- nrow(des_t)
    dest.alpha <- as.vector(des_t%*%alpha)
  }   else dest.alpha <- 0

  # Simulation
  set.seed(seed)
  ut1  <- runif(n) # uniform variables
  U <- rtp3(n = n.clust, 2*skew*sigma*sqrt(2/pi), sigma, skew, FUN = rnorm, param = "eps")  # random effects
  U.ID <- as.vector(U[ID])

  # Mixed structure I
  if(restr == "I"){
    exp.mt  <- as.vector( exp( -dest.alpha ) )
    exp.mdif <- as.vector( exp( dest.alpha - des.beta - U.ID)  )
  }

  # Mixed structure II
  if(restr == "II"){
    exp.mt  <- as.vector( exp( -dest.alpha - U.ID) )
    exp.mdif <- as.vector( exp( dest.alpha - des.beta )  )
  }

  # PGW baseline hazard
  if(distr == "PGW"){
    times <-  qpgw( 1 - exp( log(1-ut1)*exp.mdif ), theta[1], theta[2], theta[3])*exp.mt
  }
  # LN baseline hazard
  if(distr == "LN"){
    times <-  qlnorm( 1 - exp( log(1-ut1)*exp.mdif ), theta[1], theta[2])*exp.mt
  }
  # LL baseline hazard
  if(distr == "LL"){
    times <-  qllogis( 1 - exp( log(1-ut1)*exp.mdif ), theta[1], theta[2])*exp.mt
  }
  # gamma baseline hazard
  if(distr == "gamma"){
    times <-  qgamma( 1 - exp( log(1-ut1)*exp.mdif ), shape = theta[1], scale = theta[2])*exp.mt
  }

  # Corresponding times and vital status
  status <- as.vector(times < cens)
  times <- as.vector(ifelse(status, times, cens))
  status <- as.numeric(status)
  obj <- list(times = as.vector(times), status = status, ID = ID)
  return(obj)
}


########################################################################################################
# Marginal log likelihoods (original parameterisation)
########################################################################################################
# par: all model parameters in the original parameterisation
# des      : design matrix (p x n)
# des_t     : design matrix for time-dependent effects (q x n)
# status : vital status
# times : survival times
# n.clust : number of clusters


#-------------------------------------------------------------------
#' MEGH I. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix for hazard-level effects (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with Gamma baseline hazard (original parameterisation)
#' @export

tpnml1_gamma <- function(par, times, status, des, des_t, n.clust){
  sigma = par[1]; skew <- par[2]; param = par[-c(1:2)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = i) +
                         dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = i) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix for hazard-level effects (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with Log-Normal baseline hazard (original parameterisation)
#' @export
tpnml1_LN <- function(par, times, status, des, des_t, n.clust){
  sigma = par[1]; skew <- par[2]; param = par[-c(1:2)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_LN(param, u = t, times, status, des, des_t, index = i) +
                         dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_LN(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = i) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix for hazard-level effects (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with Log-logistic baseline hazard (original parameterisation)
#' @export
tpnml1_LL <- function(par, times, status, des, des_t, n.clust){
  sigma = par[1]; skew <- par[2]; param = par[-c(1:2)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_LL(param, u = t, times, status, des, des_t, index = i) +
                         dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_LL(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = i) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix for hazard-level effects (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with PGW baseline hazard (original parameterisation)
#' @export
tpnml1_PGW <- function(par, times, status, des, des_t, n.clust){
  sigma = par[1]; skew <- par[2]; param = par[-c(1:2)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = i) +
                         dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = i) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix for hazard-level effects (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with Gamma baseline hazard (original parameterisation)
#' @export
tpnml2_gamma <- function(par, times, status, des, des_t, n.clust){
  sigma = par[1]; skew <- par[2]; param = par[-c(1:2)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = i) +
                         dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = i) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix for hazard-level effects (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with Log-Normal baseline hazard (original parameterisation)
#' @export
tpnml2_LN <- function(par, times, status, des, des_t, n.clust){
  sigma = par[1]; skew <- par[2]; param = par[-c(1:2)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_LN(param, u = t, times, status, des, des_t, index = i) +
                         dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_LN(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = i) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix for hazard-level effects (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with Log-logistic baseline hazard (original parameterisation)
#' @export
tpnml2_LL <- function(par, times, status, des, des_t, n.clust){
  sigma = par[1]; skew <- par[2]; param = par[-c(1:2)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_LL(param, u = t, times, status, des, des_t, index = i) +
                         dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_LL(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = i) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix for hazard-level effects (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with PGW baseline hazard (original parameterisation)
#' @export
tpnml2_PGW <- function(par, times, status, des, des_t, n.clust){
  sigma = par[1]; skew <- par[2]; param = par[-c(1:2)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = i) +
                         dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = i) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


########################################################################################################
# Marginal log likelihoods (reparameterised)
########################################################################################################
# par: all model parameters reparameterised to the real line
# des      : design matrix (p x n)
# des_t     : design matrix for time-dependent effects (q x n)
# status : vital status
# times : survival times


#-------------------------------------------------------------------
#' MEGH I. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix for hazard-level effects (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with Gamma baseline hazard (reparameterised)
#' @export
tpnml1_gamma_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-2)
  sigma = exp(par[1]); skew <- tanh(par[2]); param[1:2] = exp(par[3:4]); param[-c(1:2)] = par[-c(1:4)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = i) +
                         dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = i) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix for hazard-level effects (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with Log-Normal baseline hazard (reparameterised)
#' @export
tpnml1_LN_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-2)
  sigma = exp(par[1]); skew <- tanh(par[2]); param[1] = par[3]; param[2] = exp(par[4]); param[-c(1:2)] = par[-c(1:4)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_LN(param, u = t, times, status, des, des_t, index = i) +
                         dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_LN(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = i) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix for hazard-level effects (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with Log-logistic baseline hazard (reparameterised)
#' @export
tpnml1_LL_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-2)
  sigma = exp(par[1]); skew <- tanh(par[2]); param[1] = par[3]; param[2] = exp(par[4]); param[-c(1:2)] = par[-c(1:4)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_LL(param, u = t, times, status, des, des_t, index = i) +
                         dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_LL(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = i) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix for hazard-level effects (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with PGW baseline hazard (reparameterised)
#' @export
tpnml1_PGW_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-2)
  sigma = exp(par[1]); skew <- tanh(par[2]); param[1:3] = exp(par[3:5]); param[-c(1:3)] = par[-c(1:5)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = i) +
                         dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = i) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix for hazard-level effects (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with Gamma baseline hazard (reparameterised)
#' @export
tpnml2_gamma_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-2)
  sigma = exp(par[1]); skew <- tanh(par[2]); param[1:2] = exp(par[3:4]); param[-c(1:2)] = par[-c(1:4)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = i) +
                         dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = i) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix for hazard-level effects (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with Log-Normal baseline hazard (reparameterised)
#' @export
tpnml2_LN_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-2)
  sigma = exp(par[1]); skew <- tanh(par[2]); param[1] = par[3]; param[2] = exp(par[4]); param[-c(1:2)] = par[-c(1:4)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_LN(param, u = t, times, status, des, des_t, index = i) +
                         dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_LN(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = i) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix for hazard-level effects (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with Log-logistic baseline hazard (reparameterised)
#' @export
tpnml2_LL_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-2)
  sigma = exp(par[1]); skew <- tanh(par[2]); param[1] = par[3]; param[2] = exp(par[4]); param[-c(1:2)] = par[-c(1:4)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_LL(param, u = t, times, status, des, des_t, index = i) +
                         dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_LL(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = i) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix for hazard-level effects (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param n.clust : number of clusters
#' @return Marginal log likelihood with PGW baseline hazard (reparameterised)
#' @export
tpnml2_PGW_rep <- function(par, times, status, des, des_t, n.clust){
  param = rep(0,length(par)-2)
  sigma = exp(par[1]); skew <- tanh(par[2]); param[1:3] = exp(par[3:5]); param[-c(1:3)] = par[-c(1:5)];
  val <- vector()
  for(i in 1:n.clust){
    tempf <- Vectorize(function(t) log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = i) +
                         dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
    u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
    lc.fun <- Vectorize(function(t) exp( log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = i) -
                                           log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = i) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
    val[i] <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = i)
  }
  return(sum(val))
}


########################################################################################################
# Individual Marginal log likelihoods (original parameterisation)
########################################################################################################
# par: all model parameters in the original parameterisation
# des      : design matrix (p x n)
# des_t     : design matrix for time-dependent effects (q x n)
# status : vital status
# times : survival times
# index : cluster ID


#-------------------------------------------------------------------
#' MEGH I. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log-likelihood with Gamma baseline hazard (original parameterisation)
#' @export

tpnml1Ind_gamma <- function(par, times, status, des, des_t, index){
  sigma = par[1]; skew <- par[2]; param = par[-c(1:2)];
  tempf <- Vectorize(function(t) log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = index) +
                       dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = index) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log-likelihood with Log-Normal baseline hazard (original parameterisation)
#' @export
tpnml1Ind_LN <- function(par, times, status, des, des_t, index){
  sigma = par[1]; skew <- par[2]; param = par[-c(1:2)];
  tempf <- Vectorize(function(t) log_likCInd1_LN(param, u = t, times, status, des, des_t, index = index) +
                       dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_LN(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = index) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log-likelihood with Log-logistic baseline hazard (original parameterisation)
#' @export
tpnml1Ind_LL <- function(par, times, status, des, des_t, index){
  sigma = par[1]; skew <- par[2]; param = par[-c(1:2)];
  tempf <- Vectorize(function(t) log_likCInd1_LL(param, u = t, times, status, des, des_t, index = index) +
                       dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_LL(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = index) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log-likelihood with PGW baseline hazard (original parameterisation)
#' @export
tpnml1Ind_PGW <- function(par, times, status, des, des_t, index){
  sigma = par[1]; skew <- par[2]; param = par[-c(1:2)];
  tempf <- Vectorize(function(t) log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = index) +
                       dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = index) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log-likelihood with Gamma baseline hazard (original parameterisation)
#' @export
tpnml2Ind_gamma <- function(par, times, status, des, des_t, index){
  sigma = par[1]; skew <- par[2]; param = par[-c(1:2)];
  tempf <- Vectorize(function(t) log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = index) +
                       dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = index) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log-likelihood with Log-Normal baseline hazard (original parameterisation)
#' @export
tpnml2Ind_LN <- function(par, times, status, des, des_t, index){
  sigma = par[1]; skew <- par[2]; param = par[-c(1:2)];
  tempf <- Vectorize(function(t) log_likCInd2_LN(param, u = t, times, status, des, des_t, index = index) +
                       dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_LN(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = index) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log-likelihood with Log-logistic baseline hazard (original parameterisation)
#' @export
tpnml2Ind_LL <- function(par, times, status, des, des_t, index){
  sigma = par[1]; skew <- par[2]; param = par[-c(1:2)];
  tempf <- Vectorize(function(t) log_likCInd2_LL(param, u = t, times, status, des, des_t, index = index) +
                       dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_LL(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = index) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters in the original parameterisation
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log-likelihood with PGW baseline hazard (original parameterisation)
#' @export
tpnml2Ind_PGW <- function(par, times, status, des, des_t, index){
  sigma = par[1]; skew <- par[2]; param = par[-c(1:2)];
  tempf <- Vectorize(function(t) log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = index) +
                       dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = index) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


########################################################################################################
# Individual Marginal log likelihoods (reparameterised)
########################################################################################################
# par: all model parameters reparameterised to the real line
# des      : design matrix (p x n)
# des_t     : design matrix for time-dependent effects (q x n)
# status : vital status
# times : survival times
# index : cluster ID

#-------------------------------------------------------------------
#' MEGH I. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log-likelihood with Gamma baseline hazard (reparameterised)
#' @export
tpnml1Ind_gamma_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-2)
  sigma = exp(par[1]); skew <- tanh(par[2]); param[1:2] = exp(par[3:4]); param[-c(1:2)] = par[-c(1:4)];
  tempf <- Vectorize(function(t) log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = index) +
                       dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_gamma(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = index) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_gamma(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log-likelihood with Log-Normal baseline hazard (reparameterised)
#' @export
tpnml1Ind_LN_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-2)
  sigma = exp(par[1]); skew <- tanh(par[2]); param[1] = par[3]; param[2] = exp(par[4]); param[-c(1:2)] = par[-c(1:4)];
  tempf <- Vectorize(function(t) log_likCInd1_LN(param, u = t, times, status, des, des_t, index = index) +
                       dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_LN(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = index) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LN(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log-likelihood with Log-logistic baseline hazard (reparameterised)
#' @export
tpnml1Ind_LL_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-2)
  sigma = exp(par[1]); skew <- tanh(par[2]); param[1] = par[3]; param[2] = exp(par[4]); param[-c(1:2)] = par[-c(1:4)];
  tempf <- Vectorize(function(t) log_likCInd1_LL(param, u = t, times, status, des, des_t, index = index) +
                       dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_LL(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = index) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_LL(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH I. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log-likelihood with PGW baseline hazard (reparameterised)
#' @export
tpnml1Ind_PGW_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-2)
  sigma = exp(par[1]); skew <- tanh(par[2]); param[1:3] = exp(par[3:5]); param[-c(1:3)] = par[-c(1:5)];
  tempf <- Vectorize(function(t) log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = index) +
                       dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd1_PGW(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = index) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd1_PGW(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log-likelihood with Gamma baseline hazard (reparameterised)
#' @export
tpnml2Ind_gamma_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-2)
  sigma = exp(par[1]); skew <- tanh(par[2]); param[1:2] = exp(par[3:4]); param[-c(1:2)] = par[-c(1:4)];
  tempf <- Vectorize(function(t) log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = index) +
                       dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_gamma(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = index) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_gamma(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log-likelihood with Log-Normal baseline hazard (reparameterised)
#' @export
tpnml2Ind_LN_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-2)
  sigma = exp(par[1]); skew <- tanh(par[2]); param[1] = par[3]; param[2] = exp(par[4]); param[-c(1:2)] = par[-c(1:4)];
  tempf <- Vectorize(function(t) log_likCInd2_LN(param, u = t, times, status, des, des_t, index = index) +
                       dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_LN(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = index) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LN(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log-likelihood with Log-logistic baseline hazard (reparameterised)
#' @export
tpnml2Ind_LL_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-2)
  sigma = exp(par[1]); skew <- tanh(par[2]); param[1] = par[3]; param[2] = exp(par[4]); param[-c(1:2)] = par[-c(1:4)];
  tempf <- Vectorize(function(t) log_likCInd2_LL(param, u = t, times, status, des, des_t, index = index) +
                       dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_LL(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = index) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_LL(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#-------------------------------------------------------------------
#' MEGH II. twopiece normal random effects.
#-------------------------------------------------------------------
#' @param par: all model parameters (reparameterised)
#' @param des : design matrix (p x n)
#' @param des_t : design matrix for time-dependent effects (q x n)
#' @param status : vital status
#' @param times : survival times
#' @param index : cluster ID
#' @return Individual Marginal log-likelihood with PGW baseline hazard (reparameterised)
#' @export
tpnml2Ind_PGW_rep <- function(par, times, status, des, des_t, index){
  param = rep(0,length(par)-2)
  sigma = exp(par[1]); skew <- tanh(par[2]); param[1:3] = exp(par[3:5]); param[-c(1:3)] = par[-c(1:5)];
  tempf <- Vectorize(function(t) log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = index) +
                       dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps", log = TRUE) )
  u.max <- optimize(tempf, lower = -10, upper = 10, maximum = TRUE)$maximum
  lc.fun <- Vectorize(function(t) exp( log_likCInd2_PGW(param, u = t, times, status, des, des_t, index = index) -
                                         log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = index) )*dtp3(t,2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = dnorm, param ="eps") )
  val <- log(integrate(lc.fun,-10,10)$value) + log_likCInd2_PGW(param, u = u.max, times, status, des, des_t, index = index)
  return(sum(val))
}


#########################################################################################################
# Marginal survival function for cluster "index" at time t
#########################################################################################################
# t: time
# par: Marginal maximum likelihood estimator (original parameterisation)
# des         : Design matrix
# des_t        : Design matrix for time-dependent effects
# index : individual ID
# NMC : number of Monte Carlo samples to be used to produce the marginal

#-------------------------------------------------------------------
#' MEGH I. twopiece normal random effects
#-------------------------------------------------------------------
#' @param t: time
#' @param par: Marginal maximum likelihood estimator (original parameterisation)
#' @param des  : Design matrix
#' @param des_t : Design matrix for time-dependent effects
#' @param index : individual ID
#' @param NMC : number of Monte Carlo samples to be used to produce the marginal
#' @return Marginal survival function with PGW baseline hazard for cluster "index" at time t
#' @export

MSurvPGW1tpn <- function(t, par, des, des_t, index, NMC){
  sigma <- par[1]
  skew <- par[2]
  theta <- par[3:5]
  indi <- which(ID==index)
  des <- as.matrix(des); des_t <- as.matrix(des_t)
  des <- as.matrix(des[indi,])
  des_t <- as.matrix(des_t[indi,])
  n.ind <- length(indi)
  q <- dim(des_t)[2]; p <- dim(des)[2]
  alpha <- par[6:(5+q)]; beta <- par[(6+q):(5+p+q)]
  # simulated random effects for Monte Carlo integration
  u <- rtp3(NMC, 2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = rnorm, param ="eps")
  # Marginal survival
  sfit.PGW <- vector()
  for(j in 1:NMC){
    exp.mlePGW <- as.vector(exp(des_t%*%alpha))
    exp.mlePGW.dif <- as.vector(exp(des%*%beta - des_t%*%alpha + u[j]))
    sfit.PGW[j] <- mean(exp(-chpgw(t*exp.mlePGW,theta[1],theta[2],theta[3])*exp.mlePGW.dif ))
  }

  return(mean(sfit.PGW))

}



#-------------------------------------------------------------------
#' MEGH II. twopiece normal random effects
#-------------------------------------------------------------------
#' @param t: time
#' @param par: Marginal maximum likelihood estimator (original parameterisation)
#' @param des  : Design matrix
#' @param des_t : Design matrix for time-dependent effects
#' @param index : individual ID
#' @param NMC : number of Monte Carlo samples to be used to produce the marginal
#' @return Marginal survival function with PGW baseline hazard for cluster "index" at time t
#' @export

MSurvPGW2tpn <- function(t, par, des, des_t, index, NMC){
  sigma <- par[1]
  skew <- par[2]
  theta <- par[3:5]
  indi <- which(ID==index)
  des <- as.matrix(des); des_t <- as.matrix(des_t)
  des <- as.matrix(des[indi,])
  des_t <- as.matrix(des_t[indi,])
  n.ind <- length(indi)
  q <- dim(des_t)[2]; p <- dim(des)[2]
  alpha <- par[6:(5+q)]; beta <- par[(6+q):(5+p+q)]
  # simulated random effects for Monte Carlo integration
  u <- rtp3(NMC, 2*skew*sigma*sqrt(2/pi),sigma,skew, FUN = rnorm, param ="eps")
  # Marginal survival
  sfit.PGW <- vector()
  for(j in 1:NMC){
    exp.mlePGW <- as.vector(exp(des_t%*%alpha + u[j]))
    exp.mlePGW.dif <- as.vector(exp(des%*%beta - des_t%*%alpha))
    sfit.PGW[j] <- mean(exp(-chpgw(t*exp.mlePGW,theta[1],theta[2],theta[3])*exp.mlePGW.dif ))
  }

  return(mean(sfit.PGW))

}

