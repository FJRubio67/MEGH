##############################################################################################
# MEGH
##############################################################################################

# Required packages
library(mvtnorm)
library(MEGH)
library(spBayesSurv)

#********************************************************
# Design matrices and ID indicators based on real data
#********************************************************

data(LeukSurv)
#View(LeukSurv)
#dim(LeukSurv)

X <- as.matrix(cbind(scale(LeukSurv$age), LeukSurv$sex, scale(LeukSurv$wbc), scale(LeukSurv$tpi) ))
colnames(X) <- cbind("std age", "sex", "wbc", "tpi")
Xt <- as.matrix(cbind(scale(LeukSurv$age)))
n.ind <- dim(X)[1]  # number of individuals
n.clust <- length(unique(LeukSurv$district)) # number of clusters
q <- dim(Xt)[2]
p <- dim(X)[2]
r <- n.clust
ID <- as.vector(LeukSurv$district)


#True values of the parameters
Sigma <- cbind(c(1,0.5),c(0.5,1))
alpha <- 0.96 
beta <- c(1.00,0.08,0.22,0.10) 
theta <- c(0.2,1.50,3.0)
true_parameters <- c(sigma,theta,alpha,beta)
npar <- length(true_parameters)
cens <- 4 # censoring times


#********************************************************
# Simulation
#********************************************************
#*
simulated_survival <- simMEGH(seed=1234, des=X, des_t=Xt, ID=ID, 
                              alpha=alpha, beta=beta, theta=theta, Sigma=Sigma, 
                              cens=rep(cens,n.ind), restr="GENERAL", distr="PGW")  

summary(simulated_survival)

# Histogram of survival times
hist(simulated_survival$times, xlab = "times", ylab ="density", probability = TRUE, main = "")
box()

# Proportion of observed survival times
mean(simulated_survival$status)





##############################################################################################
# MEGH-I
##############################################################################################

# Required packages
library(mvtnorm)
library(MEGH)
library(spBayesSurv)

#********************************************************
# Design matrices and ID indicators based on real data
#********************************************************

data(LeukSurv)
#View(LeukSurv)
#dim(LeukSurv)

X <- as.matrix(cbind(scale(LeukSurv$age), LeukSurv$sex, scale(LeukSurv$wbc), scale(LeukSurv$tpi) ))
colnames(X) <- cbind("std age", "sex", "wbc", "tpi")
Xt <- as.matrix(cbind(scale(LeukSurv$age)))
n.ind <- dim(X)[1]  # number of individuals
n.clust <- length(unique(LeukSurv$district)) # number of clusters
q <- dim(Xt)[2]
p <- dim(X)[2]
r <- n.clust
ID <- as.vector(LeukSurv$district)


#True values of the parameters
Sigma <- 1
alpha <- 0.96 
beta <- c(1.00,0.08,0.22,0.10) 
theta <- c(0.2,1.50,3.0)
true_parameters <- c(sigma,theta,alpha,beta)
npar <- length(true_parameters)
cens <- 4 # censoring times


#********************************************************
# Simulation
#********************************************************
#*
simulated_survival <- simMEGH(seed=1234, des=X, des_t=Xt, ID=ID, 
                              alpha=alpha, beta=beta, theta=theta, Sigma=Sigma, 
                              cens=rep(cens,n.ind), restr="I", distr="PGW")  

summary(simulated_survival)

# Histogram of survival times
hist(simulated_survival$times, xlab = "times", ylab ="density", probability = TRUE, main = "")
box()

# Proportion of observed survival times
mean(simulated_survival$status)







##############################################################################################
# MEGH-II
##############################################################################################

# Required packages
library(mvtnorm)
library(MEGH)
library(spBayesSurv)

#********************************************************
# Design matrices and ID indicators based on real data
#********************************************************

data(LeukSurv)
#View(LeukSurv)
#dim(LeukSurv)

X <- as.matrix(cbind(scale(LeukSurv$age), LeukSurv$sex, scale(LeukSurv$wbc), scale(LeukSurv$tpi) ))
colnames(X) <- cbind("std age", "sex", "wbc", "tpi")
Xt <- as.matrix(cbind(scale(LeukSurv$age)))
n.ind <- dim(X)[1]  # number of individuals
n.clust <- length(unique(LeukSurv$district)) # number of clusters
q <- dim(Xt)[2]
p <- dim(X)[2]
r <- n.clust
ID <- as.vector(LeukSurv$district)


#True values of the parameters
Sigma <- 1
alpha <- 0.96 
beta <- c(1.00,0.08,0.22,0.10) 
theta <- c(0.2,1.50,3.0)
true_parameters <- c(sigma,theta,alpha,beta)
npar <- length(true_parameters)
cens <- 4 # censoring times


#********************************************************
# Simulation
#********************************************************
#*
simulated_survival <- simMEGH(seed=1234, des=X, des_t=Xt, ID=ID, 
                              alpha=alpha, beta=beta, theta=theta, Sigma=Sigma, 
                              cens=rep(cens,n.ind), restr="II", distr="PGW")  

summary(simulated_survival)

# Histogram of survival times
hist(simulated_survival$times, xlab = "times", ylab ="density", probability = TRUE, main = "")
box()

# Proportion of observed survival times
mean(simulated_survival$status)

