library(mvtnorm)

### The data simulation set-up ###
p = 2000; #number of covariates
n_test = 2.5e4 #test data, will be generated only once

load("correlationmatrix.Rdata") #loading CorrX
load("betas_for_simulation_0.75.Rdata")
load("betas_for_simulation_0.85.Rdata")
#generate test covariate matrix
set.seed(25)
X_test <- rmvnorm(n_test, mean = rep(0,p), sigma = CorrX, method = "chol")
#generate test response by a logistic regression simulation (Y sampled from bernoulli)
probs_test <- 1/(1+exp(-as.numeric(X_test %*% betas)))
probs_test1 <- 1/(1+exp(-as.numeric(X_test %*% betas1)))

set.seed(1333)
Y_test <- sapply(probs_test,function(x) rbinom(n=1,size=1,prob=x))
Y_test <- sapply(probs_test1,function(x) rbinom(n=1,size=1,prob=x))




