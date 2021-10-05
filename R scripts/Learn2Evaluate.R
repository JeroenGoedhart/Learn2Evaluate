############################################################################
####################   Learn2Evaluate   ####################################
############################################################################

# This script contains the main function "Learn2Evaluate" to estimate the predictive performance
# by employing learning curves. The function is for now limited to the AUC as performance measure. 
# Only lasso, ridge, and random forest may be used as a predictive model 
# For lasso and ridge the penalty parameter is estimated using repeated CV and for random forest the standard hyperparameters are used.
# AUC point estimates are given by the learning curve estimate at the full sample size
# For the confidence bound we first determine a good split in a training set and a test set by employing the learning curve
# This good split is determined by MSE minimization or Bias Control (the user can specify which one)

# Learn2Evaluate consists of 4  global steps:
# step 1: compute repeated hold-out estimates of AUC for different training set sizes
# step 2: fit a learning curve that explains the relationship between estimated AUC and training set size
# step 3: obtain desired details from learning curve: 
                  #1. AUC estimate at full sample size 
                  #2. training set size to construct a confidence bound
# step 4: construct confidence bound by applying repeated hold-out estimation at determined training set size



### main function Learn2Evaluate requires several side functions ###
source("side functions Learn2Evaluate.R")
## Those side functions are:
# 1.  PLCregr
# 2.  PLCrf
# 3.  aucf
# 4.  PowerLawFit
# 5.  SplineFit
# 6.  ntrain_MSE
# 7.  ntrain_bias
# 8.  ci_auc_delong
# 9.  Subs
# 10. var_AUC


### Load/install required dependend packages ###
#1. pROC: to compute the AUC 
#2. cobs: to fit a constrained regression spline learning curve
#3. glmnet: to fit a lasso or a ridge
#4. randomforestSRC: to fit a random forest

packages = c("pROC", "cobs",
             "glmnet", "randomForestSRC")
package.check <- lapply(packages,FUN = function(x) {
                        if (!require(x, character.only = TRUE)) {
                        install.packages(x, dependencies = TRUE)}
                        else {
                        library(x, character.only = TRUE)}}
)

# The function Learn2Evaluate requires the following input
# 1. X: information about covariates, should be defined as a dataframe
# 2. y: response variable, should be binary
# 3. nmin: the minimum training set size, defaults to 20
# 4. nev: the number of training set sizes for which repeated hold-out is performed, defaults to 10
# 5. nrep: the number of repetitions of each repeated hold-out estimation, defaults to 50
# 6. nrepcv: the number of repetitions of cross-validation. This is only used for lasso and ridge (to estimate penalty parameter).
# 7. curve: fit of the learning curve. Can be either "IPL" (=inverse power law) or "CRS" (=constrained regression spline)
# 8. learner: The learner you want to use. For now, the only possibilities are: 1. "ridge", 2. "lasso", 3. "rf" (=random forest).
# 9. method: the method to determine the optimal training set size to construct a confidence bound. Can be either "bias" (=allow 2% margin on empirical bias) or "MSE" (= minimize MSE w.r.t. training set size)
# 10. Plot: do you want to plot the learning curve (yes = TRUE, no = FALSE)


# Running the function yields the following output in a list:
# [[1]]: AUC point estimate
# [[2]]: AUC confidence bound [1]: without bias correction, [2]: with bias correction
# [[3]]: the training set size at which the confidence bounds are constructed
# [[4]]: the learning trajectory, i.e. the repeated hold-out AUC point estimates at different training set sizes

#####################################################################################################################################
#####################################################################################################################################
Learn2Evaluate <- function(X=NULL, Y=NULL, nmin = 20, nev = 10, nrep = 50, nrepcv = 5, 
                           curve = "IPL", learner = "rf", method = "MSE", Plot = TRUE) {
  
  # control statements
  if(length(X)==0) {stop("Covariates not defined")}
  if(length(Y)==0) {stop("Response not defined")}
  if (!(is.data.frame(X))) {stop("Data should be a data frame")}
  if(is.factor(Y)) {Y <-as.numeric(Y)-1}
  if(!all(Y==1 | Y==0)) {stop("Response should be numeric and a 0 - 1 variable")}
  if(!nrow(X)==length(Y)) {stop("Sample size response and covariates should match")}
  #if(nrow(X) < 70) {stop("Sample size should be at least N = 70")}
  if(nev < 6) {stop("Number of different training set sizes should be more than 5")}
  if(nrep < 20) {stop("Number of repeats should be at least 20")}
  if(!(curve == "IPL" || curve =="CRS")) {stop("Curve fit not correctly specified, should be IPL or CRS")}
  if(!(learner == "lasso" || learner =="ridge" || learner =="rf")) {stop("Learner not correctly specified, should be ridge, lasso, or rf")}
  if(!(method == "MSE" || method =="Bias")) {stop("Method not correctly specified, should be MSE or Bias")}
  
  
  #####  STEP 1  #####
  
  print("Step 1: Computing repeated hold-out estimates at different training set sizes")
  ## computing the predicted outcomes on test set (=hold-out set) for different training set sizes and for each training set size nrep times
  if(learner == "lasso") {
    plcres0 <-PLCregr(Y=Y,X=as.matrix(X), nmin = nmin, nev=nev, nmax=nrow(X)-10, nrep=nrep, nrepcv=nrepcv, alpha=1)}
  
  if(learner == "ridge") {
    plcres0 <- PLCregr(Y=Y,X=as.matrix(X), nmin=nmin, nev=nev, nmax=nrow(X)-10, nrep=nrep, nrepcv=nrepcv, alpha=0)}
  
  if (learner == "rf") {
    plcres0 <- PLCrf(Y = Y, X = X, nmin = nmin, nev = nev, nmax = nrow(X) - 10, nrep=nrep)
  }

  plcres <- plcres0[[1]]   # a list of lists containing for each nev (list 1) and nrep (list 2) the predictions and responses of the hold-out set
  nseqtr <- plcres0[[2]]   # contains the different training set sizes
  aucs <- lapply(plcres,function(ex) unlist(lapply(ex,aucf, sm=F))) #computes for each response, prediction combination the AUC
  
  ## computing the average AUC (repeated hold-out estimate based on nrep AUCs) for each training set size
  aucsmn <- unlist(lapply(aucs,mean)) 
  allaucsmn <- cbind(ntrain = nseqtr, AUC = aucsmn) #data points to fit the learning curve: column1 consists of the training set sizes and column 2 consists of the corresponding AUC estimates
  
  
  #####  STEP 2  #####
  
  print("Step 2: Fitting the learning curve")
  if (curve == "IPL") {
  LearningCurve <- PowerlawFit(ntrain = allaucsmn[,1], aucs = allaucsmn[,2], nmin = nmin, nmax = nrow(X))
  }
  if (curve == "CRS") {
    LearningCurve <- SplineFit(ntrain = allaucsmn[,1], aucs = allaucsmn[,2], constraint = c("increase", "concave"), nmin = nmin, nmax = nrow(X))
  }
  
  
  #####  STEP 3  #####
  
  print("Step 1: Obtaining details from learning curve")
  ## 1. Obtaining AUC point estimate from learning curve  
  point_estimate <- max(LearningCurve[,2]) 
  
  
  ## 2. Determining a good training set size to construct a confidence bound
  # done by either MSE minimization (method = MSE) or bias control (method = Bias)

  if (method == "MSE"){
    opt_n <- ntrain_MSE(learning_curve = LearningCurve ,n_proportion = length(which(Y==1))/length(Y))
  }
  if (method == "Bias"){
    opt_n <- ntrain_bias(learning_curve = LearningCurve, bias = 0.02)
  }
  # opt_n contains the determined training set size with the corresponding learning curve AUC estimate at this training set size
  
  
  #####  STEP 4  #####
  
  print("Step 4: Constructing the confidence bound")
  # definining training sets of size n_opt 
  print("Determining Confidence Bound")
  subsamp <-Subs(Y=Y, model = "logistic", balance = T, 
                 ntrain =opt_n[1], fixedsubs = T,nrepeat = nrep)
  ci_aucs <-c() #empty vector to store all confidence bounds
  if (learner=="ridge"){
    X <- as.matrix(X)
    # determining penalty parameter by repeated cross-validation
    lams <- c()
    for(k in 1:nrepcv){
      samout <- subsamp[[k]]
      Xtr <- X[-samout,]
      respin <- Y[-samout]
      cvreg <- cv.glmnet(Xtr,respin,alpha=0,family="binomial", nfolds=10, type.measure="deviance")
      lams <- c(lams,cvreg$lambda.min)
    }
    lam <- median(lams) #estimated penalty parameter
    
    # determining final confidence bound by training and testing model with determined penalty parameter
    for(k in 1:nrep){
      samout <- subsamp[[k]]
      Xtr <- X[-samout,] #training
      Xte <- X[samout,]  #testing
      respin <- Y[-samout] #training
      respout <- Y[samout] #testing
      regfit <- glmnet(Xtr,respin,alpha=0,family="binomial") #fitting the model
      pred <- as.numeric(predict(regfit, newx=Xte,s=lam, type="response")) #obtaining predictions on test set
      resppred <-cbind(respout,pred)
      ci_aucs[k] <- ci_auc_delong(resppred = resppred,sm=F,alpha = 0.05) #computing lower confidence bound for auc
    }
  }
  if (learner=="lasso"){
    X <- as.matrix(X)
    # determining penalty parameter
    lams <- c() #empty vector to store penalty parameters (nrepcv in total)
    for(k in 1:nrepcv){
      samout <- subsamp[[k]]
      Xtr <- X[-samout,]
      respin <- Y[-samout]
      cvreg <- cv.glmnet(Xtr,respin,alpha=1,family="binomial", nfolds=10, type.measure="deviance")
      lams <- c(lams,cvreg$lambda.min)
    }
    lam <- median(lams) #estimated penalty parameter
    
    # determining final confidence bound by training and testing model with determined penalty parameter
    for(k in 1:nrep){
      samout <- subsamp[[k]]
      Xtr <- X[-samout,] #training
      Xte <- X[samout,]  #testing
      respin <- Y[-samout] #training
      respout <- Y[samout] #testing
      regfit <- glmnet(Xtr,respin,alpha=1,family="binomial") #fitting the model
      pred <- as.numeric(predict(regfit, newx=Xte,s=lam, type="response")) #obtaining predictions on test set
      resppred <-cbind(respout,pred)
      ci_aucs[k] <- ci_auc_delong(resppred = resppred,sm=F,alpha = 0.05) #computing lower confidence bound for AUC
    }
  }
    
  if (learner=="rf"){
    Y <- factor(Y) #used random forest package requires binary variable to be a factor
    for (k in 1:nrep) {
      samout <- subsamp[[k]]
      Xtr <- X[-samout,] #training
      Xte <- X[samout,]  #testing
      respin <- Y[-samout] #training
      respout <- Y[samout] #testing
      databoth <- data.frame(respin,Xtr) #required for random forest package
      if(ncol(X)>=1000) mtryp <- sqrt(ncol(X)) else mtryp <- ncol(X)/3
      rrfit <- rfsrc(respin ~ ., mtry=mtryp, var.used="all.trees",ntree=100, data = databoth, importance="none") #fitting the model
      predrf <- predict(rrfit, newdata = Xte, importance = F)$predicted[,2] #obtaining predictions on test set 
      resppred <-cbind(respout,predrf)
      ci_aucs[k] <- ci_auc_delong(resppred = resppred,sm=F,alpha = 0.05) #computing lower confidence bound for AUC
    }
  }
  ci_auc <- mean(ci_aucs); names(ci_auc) <- "ConfBound"
  ci_auc_bc <- ci_auc +(point_estimate - opt_n[2]) ; names(ci_auc_bc) <- "BC_ConfBound"
  conf_bound <- c(ci_auc,ci_auc_bc)
  
  ###  Plotting The Learning Curve and AUC estimates  ###
  if (Plot == T) {
  p <- plot(allaucsmn[,1],allaucsmn[,2], xlim=c(nmin,nrow(X)), ylim=c(0.5,1), pch = 16,
            xlab="Training Size", ylab="AUC", family = "serif", font.lab=2, font=2, cex.lab=1.1) #plotting the learning trajectory
  p <- lines(LearningCurve[,1],LearningCurve[,2], col="red", lwd=2.5) #plotting the fitted learning curve
  p <- points(x=nrow(X),y=point_estimate, col = "black", pch =8, lwd=2) #plotting the point estimate of Learn2Evaluate
  p <- points(x=nrow(X),y=ci_auc,col = "black", pch =24, bg="black") #plotting confidence bound of Learn2Evaluate
  p <- points(x=nrow(X),y=ci_auc_bc ,col = "black", pch =25, bg="black") #plotting bias corrected confidence bound
  }
  
  
  
  results <- list(AUC = point_estimate, ConfidenceBounds = conf_bound, TrainingSetSize = opt_n[1], LearningTrajectory = allaucsmn)      
  return(results)
}
