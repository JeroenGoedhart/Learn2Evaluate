##########################################################################
### Computing the predictive performance estimates using the bootstrap ###
##########################################################################


library(glmnet)
library(randomForestSRC)
library(pROC)
library(mvtnorm)
library(multiridge)
n_train=200; p=2000; 
B_replicates = 500
learner = "lasso"

### function to compute AUCs ###
aucf <- function(resppred, sm=F){
  #resppred <- rfres[[1]][[1]];sm=F
  resp <- resppred[,1]; pred <-resppred[,2]
  ci <- as.numeric(pROC::auc(resp,pred,smooth=sm)) 
  return(ci) 
}

### function to compute bootstrapped bound ###
bootstrapped_auc <- function(B_replicates, learner,rndnumber,conf.level){
  set.seed(rndnumber)
  X <- rmvnorm(n_train, mean = rep(0,p), sigma = CorrX, method = "chol")
  if (learner=="ridge"){
    XXblocks <- X %*% t(X)
    }
  #simulate response by bernoulli distribution (logistic regression setting)
  probs <- 1/(1+exp(-as.numeric(X %*% betas1)))
  set.seed(rndnumber)
  Y <- sapply(probs,function(x) rbinom(n=1,size=1,prob=x))
  aucs <- c() #empty vector to store all bootstrapped aucs
  for (i in 1:B_replicates) {
    #defining a bootstrap replicate
    set.seed(round(runif(1,0,1)*i^2*B_replicates+5*runif(1,0,1)*sqrt(i)+rexp(1,3)*i-100))
    boot <- sample(n_train, n_train, replace = TRUE)
    Xtrain <- X[boot,]
    Ytrain <- Y[boot]
    Xtest <- X[-boot,]
    Ytest <- Y[-boot]
    if (learner=="ridge"){
      XXTtr <-  list(XXblocks[boot,boot])
      }
    if (learner=="ridge"){
      XXTtest <- list(XXblocks[-boot,boot])
      }
    if (learner=="ridge"){
      #find optimal lambda
      ridge <- fastCV2(XXTtr,Ytrain, model="logistic",traceCV=FALSE, kfold = 10)
      optlambdas <- ridge$lambdas
      #augment covariates with optimal lambda
      XXTpen <- SigmaFromBlocks(XXTtr,penalties=optlambdas)
      #Fit model
      fit <- IWLSridge(XXTpen,Y=Ytrain, model="logistic")
      #predictions
      XXTtestpen <- SigmaFromBlocks(XXTtest,penalties=optlambdas)
      predlinpred <- predictIWLS(fit,Sigmanew=XXTtestpen)
      #model performance evaluation
      aucs[i] <- Scoring(predlinpred, Ytest, model = "logistic", score = ,"auc", print = F)
    }
    if (learner=="lasso"){
      cvreg <- cv.glmnet(Xtrain,Ytrain,alpha=1,family="binomial", nfolds=10, type.measure="deviance")
      lam <- cvreg$lambda.min
      regfit <- glmnet(Xtrain,Ytrain,alpha=1,family="binomial")
      pred <- as.numeric(predict(regfit, newx=Xtest,s=lam, type="response"))
      resppred <-cbind(Ytest,pred)
      aucs[i] <- aucf(resppred = resppred, sm=F)
    }
    if (learner=="rf"){
      Ytrain <- factor(Ytrain); Ytest <- factor(Ytest)
      databoth <- data.frame(Ytrain,Xtrain)
      if(p>=1000) mtryp <- sqrt(p) else mtryp <- p/3
      rrfit <- rfsrc(Ytrain ~ ., mtry=mtryp, var.used="all.trees",ntree=100, data = databoth, importance="none")
      pred <- predict(rrfit, newdata = data.frame(Xtest), importance = F)$predicted[,2]
      resppred <-cbind(Ytest,pred)
      aucs[i] <- aucf(resppred = resppred, sm=F)
    }
  }
  auc_avg <- mean(aucs)
  conf.bound <- quantile(aucs,conf.level)
  estimates <- c(auc_avg,conf.bound)
  names(estimates) <- c("point estimate","confidence bound")
  return(list(aucs,estimates))
}

library(foreach)
library(doParallel)

### Description how to use seeds ###
load("seeds___.Rdata") #load seeds into environment to simulate the same data. Choose for ___ the desired simulation setting.
# Seeds are structures as follows:
# seeds_x_y_z:  1. x denotes high or low signal. 0.75 corresponds to lambda=0,001 (low signal) and 0.85 corresponds to lambda=0.01 (high signal)
#               2. y denotes the learner
#               3. z denotes the samples size, either N=100, or N=200.
# Confidence bounds and corresponding true aucs should be computed using the same seeds.



expo = c('aucf','bootstrapped_auc')
nc = detectCores()
cl = makeCluster(nc-1)
registerDoParallel(cl)

ci_boot <- foreach(i = 1:length(seeds), .packages = c('mvtnorm','glmnet','multiridge','pROC','randomForestSRC'), 
                   .export = ls(globalenv()), .inorder = TRUE)  %dopar% {
                     conf_bootstrap <- bootstrapped_auc(B_replicates = B_replicates, rndnumber = seeds2[i],
                                                        learner = learner, conf.level = 0.05)[[2]]
                   }
stopCluster(cl)

###########################################################################################################################################3

