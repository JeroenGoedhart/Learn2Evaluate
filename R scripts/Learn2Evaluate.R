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
source("side functions Learn2Evaluate.R") #set working directory to map containing source scripts
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
#5. R.utils: progress bar

packages = c("pROC","cobs","glmnet", "randomForestSRC","R.utils")
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
# 10. measure: which performance measure is of interest. Can be either "AUC" or "PMSE" (mean square error of prediction)


# Running the function yields the following output in a list:
# [[1]]: AUC point estimate
# [[2]]: AUC confidence bound [1]: without bias correction, [2]: with bias correction
# [[3]]: the training set size at which the confidence bounds are constructed
# [[4]]: the learning trajectory, i.e. the repeated hold-out AUC point estimates at different training set sizes
# [[5]]: large matrix of data points that can be used to plot the learning curve
# [[6]]: the used details of Learn2Evaluate, it specifies the number of repeats of each subsample size, the method to fit the learning curve, the method to determine the optimal training size, the learner, and the performance metric.

#####################################################################################################################################
#####################################################################################################################################
Learn2Evaluate <- function(X=NULL, Y=NULL, nmin = 20, nev = 10, nrep = 50, nrepcv = 5, 
                           curve = "IPL", learner = "lasso", method = "MSE", measure = "AUC") {
  
  # control statements
  if(length(X)==0) {stop("Covariates not defined")}
  if(length(Y)==0) {stop("Response not defined")}
  if (!(is.data.frame(X))) {stop("Data should be a data frame")}
  if(is.factor(Y)) {Y <-as.numeric(Y)-1}
  #if(!all(Y==1 | Y==0)) {stop("Response should be numeric and a 0 - 1 variable")}
  if(!nrow(X)==length(Y)) {stop("Sample size response and covariates should match")}
  if(nrow(X) < 70) {stop("Sample size should be at least N = 70")}
  if(nev < 6) {stop("Number of different training set sizes should be more than 5")}
  if(nrep < 20) {stop("Number of repeats should be at least 20")}
  if(!(curve == "IPL" || curve =="CRS")) {stop("Curve fit not correctly specified, should be IPL or CRS")}
  if(!(learner == "lasso" || learner =="ridge" || learner =="rf")) {stop("Learner not correctly specified, should be ridge, lasso, or rf")}
  if(!(method == "MSE" || method =="Bias")) {stop("Method not correctly specified, should be MSE or Bias")}
  if(!(measure == "PMSE" || measure =="AUC")) {stop("Performance measure not correctly specified, should be AUC or PMSE")}
  
  
  #####  STEP 1  #####
  
  nseq <- round(seq(nmin,nrow(X)-10,length.out=nev)) #subsamples to evaluate, homogeneously spread over sample size regime
  resall <- list() #empty list to store all predictions and corresponding responses
  ## computing the predicted outcomes on test set (=hold-out set) for different training set sizes and for each training set size nrep times
  pb = utils::txtProgressBar(min = 1, max = nev+1,
                             style = 3, width = 50,
                             title = 'Running Learn2Evaluate.....')
  if (measure == "AUC") {
    model = "logistic"
    distr = "binomial"
    balance = T
    metric = "increasing"
    constraints = c("increase", "concave")
  }
  if (measure == "PMSE") {
    model = "linear"
    distr = "gaussian"
    balance = F
    metric = "decreasing"
    constraints = c("decrease","convex")
  }
  
  if(learner == "lasso") {
    print("lasso")
    X <- as.matrix(X)
    for (i in 1:nev) {
      utils::setTxtProgressBar(pb, i)
      Subsamp <- Subs(Y=Y,model=model,balance=balance,ntrain=nseq[i],fixedsubs=TRUE,nrepeat=nrep)
      plcres0 <-PLCregr(Y=Y,X=X, nfolds=5,nrepcv = nrepcv, alfa=1, nrep = nrep, subs = Subsamp,
                        family = distr)
      resall[[i]] <- plcres0
      names(resall)[i] <- as.character(nseq[i])
    }
  }

  if(learner == "ridge") {
    print("ridge")
    X <- as.matrix(X)
    resall <- list()
    for (i in 1:nev) {
      utils::setTxtProgressBar(pb, i)
      Subsamp <- Subs(Y=Y,model=model,balance=balance,ntrain=nseq[i],fixedsubs=TRUE,nrepeat=nrep)
      plcres0 <-PLCregr(Y=Y,X= X, nfolds=5,nrepcv = nrepcv, alfa=0,nrep = nrep, subs = Subsamp,
                        family = distr)
      resall[[i]] <- plcres0
      names(resall)[i] <- as.character(nseq[i])
    }
  }
    
  
  if (learner == "rf") {
    print("rf")
    if(measure =="AUC"){Y <- factor(Y)}
    resall <- list()
    for (i in 1:nev) {
      utils::setTxtProgressBar(pb, i)
      Subsamp <- Subs(Y=Y,model=model,balance=balance,ntrain=nseq[i],fixedsubs=TRUE,nrepeat=nrep)
      plcres0 <-PLCrf(Y=Y,X= X, nrep = nrep, subs = Subsamp, measure = measure)
      resall[[i]] <- plcres0
      names(resall)[i] <- as.character(nseq[i])
    }
  }
  if (measure == "AUC"){
  estimates <- lapply(resall,function(ex) unlist(lapply(ex,aucf, sm=F))) #computes for each response, prediction combination the AUC
  
  ## computing the average AUC (repeated hold-out estimate based on nrep AUCs) for each training set size
  estimatesmn <- unlist(lapply(estimates,mean)) 
  LearnCurve <- cbind(ntrain = nseq, AUC = estimatesmn) #data points to fit the learning curve: column1 consists of the training set sizes and column 2 consists of the corresponding AUC estimates
  }
  
  if (measure == "PMSE") {
  estimates <- lapply(resall,function(ex) unlist(lapply(ex,pMSE)))
  estimatesmn <- unlist(lapply(estimates,mean))
  LearnCurve <- cbind(training_size = nseq, pMSE = estimatesmn)
  }
  
  #####  STEP 2  #####
  
  if (curve == "IPL") {
  LearningCurve <- PowerlawFit(ntrain = LearnCurve[,1], estimates = LearnCurve[,2], nmin = nmin, nmax = nrow(X), metric = metric)
  }
  if (curve == "CRS") {
    LearningCurve <- SplineFit(ntrain = LearnCurve[,1], estimates = LearnCurve[,2], constraint = constraints, nmin = nmin, nmax = nrow(X))
  }
  
  
  #####  STEP 3  #####
  
  ## 1. Obtaining AUC point estimate from learning curve
  if (measure == "AUC") {
  point_estimate <- max(LearningCurve[,2])}
  if (measure == "PMSE") {
  point_estimate <- min(LearningCurve[,2])}
  
  ## 2. Determining a good training set size to construct a confidence bound
  # done by either MSE minimization (method = MSE) or bias control (method = Bias)
  if (method == "Bias"){
    opt_n <- ntrain_bias(learning_curve = LearningCurve, bias = 0.02, metric = metric)
  }
  if (method == "MSE"){
    if (measure == "AUC"){
    opt_n <- ntrain_MSE_AUC(learning_curve = LearningCurve ,n_proportion = length(which(Y==1))/length(Y))
    }
    if (measure == "PMSE"){
    opt_n <- ntrain_MSE_pMSE(learning_curve = LearningCurve)
    }
  }
  
  # opt_n contains the determined training set size with the corresponding learning curve AUC estimate at this training set size
  
  
  #####  STEP 4  #####
  Subsamp <- Subs(Y=Y,model=model,balance=balance,ntrain=opt_n[1],fixedsubs=TRUE,nrepeat=nrep)
  
  cis <-c() #empty vector to store all confidence bounds
  if (learner=="ridge"){
    X <- as.matrix(X)
    # determining penalty parameter by repeated cross-validation
    lams <- c()
    for(k in 1:nrepcv){
      samout <- Subsamp[[k]]
      Xtr <- X[-samout,]
      respin <- Y[-samout]
      cvreg <- cv.glmnet(Xtr,respin,alpha=0,family=distr, nfolds=10, type.measure="deviance")
      lams <- c(lams,cvreg$lambda.min)
    }
    lam <- median(lams) #estimated penalty parameter
    
    # determining final confidence bound by training and testing model with determined penalty parameter
    for(k in 1:nrep){
      samout <- Subsamp[[k]]
      Xtr <- X[-samout,] #training
      Xte <- X[samout,]  #testing
      respin <- Y[-samout] #training
      respout <- Y[samout] #testing
      regfit <- glmnet(Xtr,respin,alpha=0,family=distr) #fitting the model
      pred <- as.numeric(predict(regfit, newx=Xte,s=lam, type="response")) #obtaining predictions on test set
      resppred <-cbind(respout,pred)
      if (measure == "AUC"){
      cis[k] <- suppressMessages(ci_auc_delong(resppred = resppred,sm=F,alpha = 0.05)) #computing lower confidence bound for auc
      }
      if (measure == "PMSE"){
        pMSE <- pMSE(resppred)
        n <- length(resppred[,1])
        cis[k] <- ci_pMSE(pMSE = pMSE, deg_free = n)[2]
      }
    }
  }
  if (learner =="lasso"){
    X <- as.matrix(X)
    # determining penalty parameter
    lams <- c() #empty vector to store penalty parameters (nrepcv in total)
    for(k in 1:nrepcv){
      samout <- Subsamp[[k]]
      Xtr <- X[-samout,]
      respin <- Y[-samout]
      cvreg <- cv.glmnet(Xtr,respin,alpha=1,family=distr, nfolds=10, type.measure="deviance")
      lams <- c(lams,cvreg$lambda.min)
    }
    lam <- median(lams) #estimated penalty parameter
    
    # determining final confidence bound by training and testing model with determined penalty parameter
    for(k in 1:nrep){
      samout <- Subsamp[[k]]
      Xtr <- X[-samout,] #training
      Xte <- X[samout,]  #testing
      respin <- Y[-samout] #training
      respout <- Y[samout] #testing
      regfit <- glmnet(Xtr,respin,alpha=1,family=distr) #fitting the model
      pred <- as.numeric(predict(regfit, newx=Xte,s=lam, type="response")) #obtaining predictions on test set
      resppred <-cbind(respout,pred)
      #computing lower confidence bound for AUC
      if (measure == "AUC"){
        cis[k] <- suppressMessages(ci_auc_delong(resppred = resppred,sm=F,alpha = 0.05)) #computing lower confidence bound for auc
      }
      if (measure == "PMSE"){
        pMSE <- pMSE(resppred)
        n <- length(resppred[,1])
        cis[k] <- ci_pMSE(pMSE = pMSE, deg_free = n)[2]
      }
    }
  }
    
  if (learner=="rf"){
    if(measure =="AUC"){Y <- factor(Y)} #used random forest package requires binary variable to be a factor
    for (k in 1:nrep) {
      samout <- Subsamp[[k]]
      Xtr <- X[-samout,] #training
      Xte <- X[samout,]  #testing
      respin <- Y[-samout] #training
      respout <- Y[samout] #testing
      databoth <- data.frame(respin,Xtr) #required for random forest package
      if(ncol(X)>=1000) {mtryp <- sqrt(ncol(X))} else {mtryp <- ncol(X)/3}
      rrfit <- rfsrc(respin ~ ., mtry=mtryp, var.used="all.trees",ntree=100, data = databoth, importance="none") #fitting the model
      if (measure == "AUC"){
      predrf <- predict(rrfit, newdata = Xte, importance = F)$predicted[,2]} #obtaining predictions on test set 
      if (measure=="PMSE"){
      predrf <- predict(rrfit, newdata = Xte, importance = F)$predicted} 
      resppred <-cbind(respout,predrf)
      #computing lower confidence bound for AUC or PMSE
      if (measure == "AUC"){
        cis[k] <- suppressMessages(ci_auc_delong(resppred = resppred,sm=F,alpha = 0.05)) #computing lower confidence bound for auc
      }
      if (measure == "PMSE"){
        pMSE <- pMSE(resppred)
        n <- length(resppred[,1])
        cis[k] <- ci_pMSE(pMSE = pMSE, deg_free = n)[2]
      }
    }
  }
  ci <- median(cis); names(ci) <- "ConfBound"
  ci_bc <- ci +(point_estimate - opt_n[2])
  names(ci_bc) <- "BC_ConfBound"
  conf_bound <- c(ci,ci_bc)
  utils::setTxtProgressBar(pb, nev+1)
  
  
  results <- list(PointEstimate = point_estimate, 
                  ConfidenceBounds = conf_bound, 
                  TrainingSetSize = opt_n[1], 
                  LearningTrajectory = LearnCurve, 
                  LearningCurve = LearningCurve, 
                  details = c(nrep = nrep, curve = curve, method = method, learner = learner, measure = measure))      
  return(results)
}
###################################################################################################################################################################################################

######  Plotting The Learning Curve and AUC repeated hold-out estimates  #####
##############################################################################

# Function to plot the curve. As input the list object of the main function Learn2Evaluate is required.
# Additionally, the input Add (=F, =T) is required.
# This specifies whether you want to add a curve to an already existing plot, so that you can compare different learners.
# Only Add=T if you already have a learning curve plotted in your R studio environment.

PlotCurve <- function(Learn2Evaluate, Add = F){
  if(Add == TRUE & is.null(dev.list()['RStudioGD'])) {stop("Cannot add to empty plot")}
  ##obtaining details from Learn2Evaluate object
  aucs <- Learn2Evaluate[[4]]
  learner <- Learn2Evaluate[[6]][4]
  LearningCurve <- Learn2Evaluate[[5]]
  point_estimate <- Learn2Evaluate[[1]]
  ConfBound <- Learn2Evaluate[[2]][1]
  BC_ConfBound <- Learn2Evaluate[[2]][2]
  nmax <- max(aucs[,1])+10
  nmin <- min(aucs[,1])
  measure <- Learn2Evaluate[[6]][5]
  if (measure == "AUC"){ylim = c(0.5,1)}
  if (measure == "PMSE"){ylim =c(0,ConfBound+0.1)}
  if (measure == "AUC"){ylab = "AUC"}
  if (measure == "PMSE"){ylab = "PMSE"}
  if (learner == "lasso") {col = "#395D9CFF"; pos = nmax+2}
  if (learner == "ridge") {col = "#0B0405FF";pos = nmax}
  if (learner == "rf") {col = "#60CEACFF";pos = nmax-2}
  
  if (Add == F) {
    p <- plot(aucs[,1],aucs[,2], xlim=c(nmin,nmax+1), ylim=ylim, pch = 16, col = col,
              xlab="Training Size", ylab=ylab, family = "serif", font.lab=2, font=2, cex.lab=1.1) #plotting the learning trajectory
    p <- lines(LearningCurve[,1],LearningCurve[,2], col=col, lwd=2.5) #plotting the fitted learning curve
    p <- points(x=nmax,y=point_estimate, col = col, pch =8, lwd=2) #plotting the point estimate of Learn2Evaluate
    p <- points(x=pos,y=ConfBound,col = col, pch =24, bg=col) #plotting confidence bound of Learn2Evaluate
    p <- points(x=pos,y=BC_ConfBound ,col = col, pch =25, bg=col) #plotting bias corrected confidence bound
  }
    
  if (Add == T) {
    p <- points(aucs[,1],aucs[,2], pch = 16, col=col) #plotting the learning trajectory
    p <- lines(LearningCurve[,1],LearningCurve[,2], col=col, lwd=2.5) #plotting the fitted learning curve
    p <- points(x=nmax,y=point_estimate, col = col, pch =8, lwd=2) #plotting the point estimate of Learn2Evaluate
    p <- points(x=pos,y=ConfBound,col = col, pch =24, bg=col) #plotting confidence bound of Learn2Evaluate
    p <- points(x=pos,y=BC_ConfBound ,col = col, pch =25, bg=col) #plotting bias corrected confidence bound
  }
}


