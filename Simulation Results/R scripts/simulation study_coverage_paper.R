############################################################################################################
### Simulation study to assess coverage of PLC methodology to construct lower bound confidence intervals ###
############################################################################################################

library(mvtnorm)   # to simulate from a multivariate normal
library(glmnet)   # to fit lasso/ridge regression
library(randomForest) #to perform random forest
library(cobs)      # to fit constrained splines 
library(pROC)      # to calculate the AUC and corresponding confidence interval
library(foreach)
library(parallel)
library(doParallel)
library(cvAUC)  #to apply the method of Le Dell et al.
library(randomForestSRC) # to fit a random forest


## The data simulation set-up: defining parameters ##
Nsim = 1000; #number of simulations, for each simulation a confidence bound is returned
n_train = 100; #number of training samples, will be generated Nsim times  
p = 2000; #number of covariates
nmin = 20; nmax=n_train; nev=10;nrep=50;nrepcv=5; ## parameters to generate PLC
learner = "rf"

source('side functions for coverage simulation.R') ##load side functions into environment


## make sure working directory is set to map containing files below ##
load("correlationmatrix.Rdata")
load("betas_for_simulation_0.85.Rdata")# higher signal beta's (lambda = 0,01 in paper)
load("betas_for_simulation_0.75.Rdata")# lower signal beta's (lambda = 0,001 in paper)

### function to compute predictive performance estimates by Learn2Evaluate and 10-fold CV usng an asymptotic variance estimates (Le Dell et al.)
plcs <- function(n_train, p, nmin, nev, nrep, nrepcv, learner, rndnumber){
  

## Step 1: simulate the data ##

# simulate design matrix
set.seed(rndnumber)
X_train <- rmvnorm(n_train, mean = rep(0,p), sigma = CorrX, method = "chol")


#simulate response by bernoulli distribution (logistic regression setting)
probs <- 1/(1+exp(-as.numeric(X_train %*% betas1))) #betas 1: higher signal
                                                    #betas 2: lower signal
set.seed(rndnumber)
Y_train <- sapply(probs,function(x) rbinom(n=1,size=1,prob=x))


## step 2: assess cv (10 fold) auc value and compute confidence interval with de laan's method ##

ci_auc_delaan <- auc_delaan(X = X_train, Y = Y_train, V = 10, learner = learner) #computes estimate by (cross-validation) and confidence interval


## step 3: create predictive learning curve ##
#Fit a model to subset of training data and make predictions on left out subset of training data

if (learner=="lasso"){
  plcres0 <- PLCregr(Y=Y_train,X=as.matrix(X_train),nmin = nmin, nev=nev, nrep=nrep, nrepcv=nrepcv, alpha=1) ##learning and prediction for lasso
}

if (learner=="ridge"){
  plcres0 <- PLCregr(Y=Y_train,X=as.matrix(X_train),nmin = nmin,nev=nev,nrep=nrep, nrepcv=nrepcv, alpha=0) ##learning and prediction for ridge
}
if (learner=="rf"){
  plcres0 <- PLCrf(Y=as.factor(Y_train),X=as.data.frame(X_train),nmin = nmin,nev=nev,nrep=nrep) ##learning and prediction for rf
}


plcres <- plcres0[[1]]
nseqtr <- plcres0[[2]]

#compute AUCs for all training sample sizes (outer lapply), and for all repeats (inner lapply)
aucs <- lapply(plcres,function(ex) unlist(lapply(ex,aucf, sm=F)))
aucsmn <- unlist(lapply(aucs,mean)) #compute mean for all training sizes

#store matrix with training sizes and corresponding aucs and confidence bounds 
allaucsmn <- cbind(Training set Size = nseqtr, AUC = aucsmn) #contains the training sizes (column 1) and corresponding aucs (column 2)
colnames(allaucsmn) <- c("Training Size", "AUC")
remove(plcres,plcres0)


## step 4: fit power law and constrained regression spline to data (= training size versus AUC) ##

#fit spline to aucs and training sizes
splines <- Splinefitting(aucs = allaucsmn[,2], ntrain = allaucsmn[,1], constraint = c("increase", "concave"))

#fit power law to aucs and training sizes
power_law <- power_law_fit(aucs = allaucsmn[,2], ntrain = allaucsmn[,1])
spline_prediction <- max(splines[,2]) # spline prediction of auc at maximal training size
powerlaw_prediction <- max(power_law[,2]) #power law prediction of auc at maximal training size
predictions <-c(spline_prediction,powerlaw_prediction)
names(predictions) <- c("spline","powerlaw")

## step 5: compute optimal training size based on: 1. bias (2%); 2. minimize MSE ##

#compute optimal training size for both the power law and spline fit
opt_n_bias <- ntrain_bias(powerlaw_results = power_law, spline_results = splines, bias = 0.02) #method 1: bias
opt_n_MSE <- ntrain_MSE(powerlaw_results = power_law, spline_results = splines, 
                        n_proportion = length(which(Y_train==1))/length(Y_train)) #method 2: minimize MSE

remove(splines, power_law, aucs,aucsmn)

## step 6: retrain the model at derived training sizes (4 in total) ##

#defining the optimal training sizes
ntrain_spline_bias <-ceiling(opt_n_bias[[1]][1]); AUC_spline_bias <-opt_n_bias[[1]][2]
ntrain_spline_mse <-ceiling(opt_n_MSE[[1]][1]); AUC_spline_MSE <-opt_n_MSE[[1]][2]
ntrain_powerl_bias <-ceiling(opt_n_bias[[2]][1]); AUC_powerl_bias <-opt_n_bias[[2]][2]
ntrain_powerl_mse <-ceiling(opt_n_MSE[[2]][1]); AUC_powerl_MSE <-opt_n_MSE[[2]][2]
optimal_ntrains <- c(ntrain_spline_bias, ntrain_spline_mse, ntrain_powerl_bias, ntrain_powerl_mse)
AUC_preds <- c(AUC_spline_bias,AUC_spline_MSE,AUC_powerl_bias,AUC_powerl_MSE)
names(optimal_ntrains) <-c("ntrain_spline_bias", "ntrain_spline_mse", 
                           "ntrain_powerl_bias", "ntrain_powerl_mse")

#retrain the models and derive lower confidence bound
if (learner=="ridge"){
  auc_and_ci <-list()
  for (i in 1:length(optimal_ntrains)){
    subsamp <-Subs(Y=Y_train, model = "logistic", balance = T, 
                   ntrain =optimal_ntrains[i], fixedsubs = T,nrepeat = 30) #creating new subsamples
    ci_aucs_i <-c() #empty vector to store all confidence bounds 
    for (k in 1:30) {
      samout<-subsamp[[k]]
      Xtr <- X_train[-samout,]
      Xte <- X_train[samout,]
      respin <- Y_train[-samout]
      respout <- Y_train[samout]
      cvreg <- cv.glmnet(Xtr,respin,alpha=0,family="binomial", nfolds=10, type.measure="deviance")
      lam <- cvreg$lambda.min
      regfit <- glmnet(Xtr,respin,alpha=0,family="binomial")
      pred <- as.numeric(predict(regfit, newx=Xte,s=lam, type="response"))
      resppred <-cbind(respout,pred)
      ci_aucs_i[k] <- ci_auc_delong(resppred = resppred,sm=F,alpha = 0.05)
    }
    ci_auc_i<-mean(ci_aucs_i) # mean of all computed ci's
    auc_and_ci[[i]]<-c(AUC_preds[i],ci_auc_i, optimal_ntrains[i])
    names(auc_and_ci[[i]]) <- c("auc","ci", "ntrain")
  }
}
  

if (learner=="lasso"){
  auc_and_ci <-list()
  for (i in 1:length(optimal_ntrains)){
    subsamp <-Subs(Y=Y_train, model = "logistic", balance = T, 
                   ntrain =optimal_ntrains[i], fixedsubs = T,nrepeat = 30) #creating new subsamples
    ci_aucs_i <-c() #empty vector to store all confidence bounds 
    for (k in 1:30) {
      samout<-subsamp[[k]]
      Xtr <- X_train[-samout,]
      Xte <- X_train[samout,]
      respin <- Y_train[-samout]
      respout <- Y_train[samout]
      cvreg <- cv.glmnet(Xtr,respin,alpha=1,family="binomial", nfolds=10, type.measure="deviance")
      lam <- cvreg$lambda.min
      regfit <- glmnet(Xtr,respin,alpha=1,family="binomial")
      pred <- as.numeric(predict(regfit, newx=Xte,s=lam, type="response"))
      resppred <-cbind(respout,pred)
      ci_aucs_i[k] <- ci_auc_delong(resppred = resppred,sm=F,alpha = 0.05)
    }
    ci_auc_i<-mean(ci_aucs_i) # mean of all computed ci's
    auc_and_ci[[i]]<-c(AUC_preds[i],ci_auc_i, optimal_ntrains[i])
    names(auc_and_ci[[i]]) <- c("auc","ci", "ntrain")
  }
}
if (learner=="rf"){
  auc_and_ci <-list()
  for (i in 1:length(optimal_ntrains)){
    subsamp <-Subs(Y=Y_train, model = "logistic", balance = T, 
                   ntrain =optimal_ntrains[i], fixedsubs = T,nrepeat = 30) #creating new subsamples
    ci_aucs_i <-c() #empty vector to store all confidence bounds 
    for (k in 1:30) {
      samout<-subsamp[[k]]
      Xtr <- X_train[-samout,]
      Xte <- X_train[samout,]
      respin <- Y_train[-samout]
      respout <- Y_train[samout]
      rffit <- randomForest(Xtr,as.factor(respin))
      pred <- predict(rffit, newdata=Xte,type="prob")[,2]
      resppred <-cbind(respout,pred)
      ci_aucs_i[k] <- ci_auc_delong(resppred = resppred,sm=F,alpha = 0.05)
    }
    ci_auc_i<-mean(ci_aucs_i) # mean of all computed ci's
    auc_and_ci[[i]]<-c(AUC_preds[i],ci_auc_i, optimal_ntrains[i])
    names(auc_and_ci[[i]]) <- c("auc","ci", "ntrain")
  }
}
remove(X_train)

names(auc_and_ci) <- c("spline bias", "spline MSE", "powerlaw bias", "powerlaw MSE")
allresults <- list(auc_and_ci, predictions,ci_auc_delaan, allaucsmn)
names(allresults) <- c("derived ci's","auc estimates at max ntrain","auc and ci De Laan's method", "learning curve data")
return(allresults) # end result is a list with 2 elements, [[1]]: predicted auc, ci and ntrain of powerlaw and spline for both bias and MSE
} 

### single computation of predictive performance estimates
learning_curve <- plcs(n_train = n_train, p=2000, nmin = 10, nev = nev, nrep = nrep, nrepcv = nrepcv, learner = learner, rndnumber = 1)
learning_curve


### parallel computation ###
############################

### Description how to use seeds ###
load("seeds___.Rdata") #load seeds into environment to simulate the same data. Choose for ___ the desired simulation setting.
# Seeds are structures as follows:
# seeds_x_y_z:  1. x denotes high or low signal. 0.75 corresponds to lambda=0,001 (low signal) and 0.85 corresponds to lambda=0.01 (high signal)
#               2. y denotes the learner
#               3. z denotes the samples size, either N=100, or N=200.
# Confidence bounds and corresponding true aucs should be computed using the same seeds.



library(foreach)
library(doParallel)
expo = c('aucf','ci_auc_delong','ntrain_bias','ntrain_MSE', 
         'PLCregr','plcs','power_law_fit','Splinefitting','Subs', 'var_AUC')
nc = detectCores()
cl = makeCluster(nc-5)
registerDoParallel(cl)

ci_plcs <- foreach(i = 1:1000, .packages = c('mvtnorm','glmnet','cobs','pROC','cvAUC','randomForest'), 
                   .export = ls(globalenv()), .inorder = TRUE)  %dopar% {
                      ci_plc <- plcs(n_train = n_train, p = p, nmin = 10, nev = nev, nrep = nrep, nrepcv = nrepcv, 
                                     learner = learner, rndnumber = seeds[i])
                   }
stopCluster(cl)


################################################################################################################################