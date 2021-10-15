###########################################################
### Computing the true AUC, using a very large test set ###
###########################################################


library(pROC)
library(cvAUC)
library(glmnet)
library(mvtnorm)
library(randomForest)

### Description how to use seeds ###
load("seeds___.Rdata") #load seeds into environment to simulate the same data. Choose for ___ the desired simulation setting.
# Seeds are structures as follows:
# seeds_x_y_z:  1. x denotes high or low signal. 0.75 corresponds to lambda=0,001 (low signal) and 0.85 corresponds to lambda=0.01 (high signal)
#               2. y denotes the learner
#               3. z denotes the samples size, either N=100, or N=200.
# Confidence bounds and corresponding true aucs should be computed using the same seeds.




load("X_test.Rdata") #load test data of covariates
load("Y_test_0.75.Rdata") #load Y_test, check whether the high or low signal situation is required
#load("Y_test_0.85.Rdata")
load("betas_for_simulation_0.75.Rdata") #load betas, check whether the high or low signal situation is required
#load("betas_for_simulation_0.85.Rdata")
n_train = 100; #number of training samples, will be generated Nsim times  
p = 2000
learner <- "lasso"
##### for lasso!!! #####

real_aucs <-c()
for (i in 1:length(seeds)){
  #generating same training data by specifying the seed
  set.seed(seeds[i])
  X_train <- rmvnorm(n_train, mean = rep(0,p), sigma = CorrX, method = "chol")
  probs <- 1/(1+exp(-as.numeric(X_train %*% betas1)))
  set.seed(seeds[i])
  Y_train <- sapply(probs,function(x) rbinom(n=1,size=1,prob=x))
  # compute real auc
  cvreg <- cv.glmnet(X_train,Y_train,alpha=1,family="binomial", nfolds=10, type.measure="deviance")
  lam <- cvreg$lambda.min
  regfit <- glmnet(X_train,Y_train,alpha=1,family="binomial")
  pred <- as.numeric(predict(regfit, newx=X_test,s=lam, type="response"))
  auc_i <-pROC::auc(Y_test,pred)[1]
  real_aucs[i]<-auc_i
}

### for ridge!!! #####
real_aucs <-c()
for (i in 1:length(seeds)){
  #generating same training data by specifying the seed
  set.seed(seeds[i])
  X_train <- rmvnorm(n_train, mean = rep(0,p), sigma = CorrX, method = "chol")
  probs <- 1/(1+exp(-as.numeric(X_train %*% betas)))
  set.seed(seeds[i])
  Y_train <- sapply(probs,function(x) rbinom(n=1,size=1,prob=x))
  # compute real auc
  cvreg <- cv.glmnet(X_train,Y_train,alpha=0,family="binomial", nfolds=10, type.measure="deviance")
  lam <- cvreg$lambda.min
  regfit <- glmnet(X_train,Y_train,alpha=0,family="binomial")
  pred <- as.numeric(predict(regfit, newx=X_test,s=lam, type="response"))
  auc_i <-pROC::auc(Y_test,pred)[1]
  real_aucs[i]<-auc_i
}

### for random forest!!! ###
real_aucs<-c()
for (i in 1:length(seeds)){
  #generating same training data by specifying the seed
  set.seed(seeds[i])
  X_train <- rmvnorm(n_train, mean = rep(0,p), sigma = CorrX, method = "chol")
  probs <- 1/(1+exp(-as.numeric(X_train %*% betas1)))
  set.seed(seeds[i])
  Y_train <- sapply(probs,function(x) rbinom(n=1,size=1,prob=x))
  Y_train <- as.factor(Y_train)
  # compute real auc
  databoth <- data.frame(Y_train,X_train)
  if(p>=1000) mtryp <- sqrt(p) else mtryp <- p/3
  rrfit <- rfsrc(Y_train ~ ., mtry=mtryp, var.used="all.trees",ntree=100, data = databoth, importance="none")
  predrf <- predict(rrfit, newdata = data.frame(X_test), importance = F)$predicted[,2]
  auc_i <-pROC::auc(Y_test,predrf)[1]
  real_aucs_0.85_rf_n200[i]<-auc_i
}