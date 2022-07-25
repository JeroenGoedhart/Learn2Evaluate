### functions used in the main function Learn2Evaluate ###

# It contains the following functions:
# subs: used to obtain balanced training subsamples (equal amount of 1's and 0's as response values )
# PLCregr: this is the function that computes all fits and predictions (so for all different training sizes and repeats) of the learner (lasso or ridge in this case)
# PLCrf: this is the function that computes all fits and predictions (so for all different training sizes and repeats) of the learner (random forest in this case)
# aucf: function that computes the area under the curve for the predicted values obtained by PLCregr
# ci_auc_delong: the function that calculates the confidence bound by using delongs nonparametric method
# PowerlawFit: when the points for the plc are obtained (auc versus ntrain), this function fits a power law model to the data
# SplineFit: additionally to the power_law_fit we also fit a smoothing B-spline with additional constraints (monotonic and concave) to the data
# ntrain_bias: function that finds the minimal training size that is within 2% of the auc predicted at the maximal training size
# var_AUC: function that calculates asymptotic variance of AUC (based on gaussian assumption) according to the result of Bamber
# ntrain_MSE: function that minimizes emperical MSE of auc, minimized with respect to training size

#Function for creating (balanced) folds
Subs <- function(Y,model="logistic",balance=TRUE,ntrain,fixedsubs=TRUE,nrepeat=1){ #response is required for balanced CV
  #response: response vector, length n
  #model: "logistic", "cox", etc
  #balance: should the splits balance levels of the response?
  #kfold: scalar, the number of folds
  #fixedsubs: should the folds be fixed? (for reproducibility)
  #nrepeat: number of repeats of the CV
  #Output: list object with kfold elements containing the sample indices of the left-out samples per fold
  
  #Y <- rbinom(100,1,.95);model="logistic";balance=FALSE;ntrain=90;fixedsubs=TRUE;nrepeat=1
  
  response <- Y
  if(model=="linear") balance <- FALSE
  subrep <- function(rep){
    nsam <- length(response)
    ntest <- nsam-ntrain
    if (fixedsubs) set.seed(3534+rep-1) else set.seed(NULL)
    if (!balance) {
      subs <- sort(sample(1:nsam,ntest))
    }
    else {
      if (model == "logistic"){
        if (class(response) == "factor")
          nev <- which((as.numeric(response) - 1) == 1)
        else nev <- which(response == 1)
      }
      if (model == "cox") nev <- which(response[, 2] == 1)
      nonev <- (1:nsam)[-nev]
      nsamev <- length(nev)
      ntestev <- min(max(1,round(ntest/nsam * nsamev)),(nsamev-1))
      ntestnonev <- ntest - ntestev
      subs <- sort(c(sample(nev,ntestev), sample(nonev,ntestnonev)))
    }
    return(subs)
  }
  return(lapply(1:nrepeat,subrep))
}


#PLCs for regularized regression (ridge, lasso, Elastic Net)
PLCregr <- function(Y,X, nfolds=5, alfa=1, nrepcv=5,nrep=50, subs = Subsamp,
                    family){ #alpha =1: lasso; alpha = 0: ridge
  lams <- c()
  for(k in 1:nrepcv){
    samout <- subs[[k]]
    Xtr <- X[-samout,]
    Xte <- X[samout,]
    respin <- Y[-samout]
    respout <- Y[samout]
    cvreg <- suppressWarnings(cv.glmnet(Xtr,respin,alpha=alfa,family=family, nfolds=nfolds, type.measure="deviance"))
    lams <- c(lams,cvreg$lambda.min)
  }
  
  lam <- median(lams)
  res <-list()
  for(k in 1:nrep){
    #k <-1
    samout <- subs[[k]]
    Xtr <- X[-samout,]
    Xte <- X[samout,]
    respin <- Y[-samout]
    respout <- Y[samout]
    regfit <- glmnet(Xtr,respin,alpha=alfa,family=family)
    pred <- as.numeric(predict(regfit, newx=Xte,s=lam, type="response"))
    res[[k]]<- cbind(respout=respout,pred=pred)
  }
  return(res)
}

#PLCs for Random Forest
PLCrf <- function(Y,X,nrep=50, subs = Subsamp, measure){
  res <-list()
  for(k in 1:nrep){
    samout <- subs[[k]]
    Xtr <- X[-samout,]
    Xte <- X[samout,]
    respin <- Y[-samout]
    respout <- Y[samout]
    databoth <- data.frame(respin,Xtr)
    if(ncol(X)>=1000) {mtryp <- sqrt(ncol(X))} else {mtryp <- ncol(X)/3}
    rrfit <- rfsrc(respin ~ ., mtry=mtryp, var.used="all.trees",ntree=100, data = databoth, importance="none")
    if (measure == "AUC"){
      pred <- predict(rrfit, newdata = Xte, importance = F)$predicted[,2] #obtaining predictions on test set 
      res[[k]]<- cbind(respout=as.numeric(respout)-1,pred=pred)}
    
    if (measure=="PMSE"){
      pred <- predict(rrfit, newdata = Xte, importance = F)$predicted
      res[[k]]<- cbind(respout=respout,pred=pred)}
    
  }
  return(res)
}

## function to compute Area under the receiver operating curve (AUC)
aucf <- function(resppred, sm=F){
  #resppred <- rfres[[1]][[1]];sm=F
  resp <- resppred[,1]; pred <-resppred[,2]
  ci <- suppressMessages(as.numeric(pROC::auc(resp,pred,smooth=sm)))  
  return(ci) 
}

## function to compute confidence bound of AUC using de Long method
ci_auc_delong <- function(resppred,sm=F, alpha=0.05){
  #resp <- plc1[[1]][[1]]$respout;pred <- plc1[[1]][[1]]$pred
  resp <- resppred[,1]; pred <-resppred[,2]
  ci <- as.numeric(pROC::ci.auc(resp,pred,smooth=sm, conf.level=1-2*alpha, method="delong")) 
  return(ci[1]) 
}

pMSE <- function(resppred){
  resp <- resppred[,1]; pred <- resppred[,2]
  pMSE <- mean((resp-pred)^2)
  return(pMSE)
}

ci_pMSE <- function(pMSE, deg_free, p_lower = 0.05, p_upper = 0.95){
  c(pred_lower = deg_free / qchisq(p_upper, df = deg_free)*pMSE,
    pred_upper = deg_free / qchisq(p_lower, df = deg_free)*pMSE)
}


## function to fit a power law
PowerlawFit <- function(ntrain, estimates, nmin, nmax, metric = "increasing"){
  X<-ntrain
  y<-estimates
  if (metric =="increasing"){
    compCost<-function(X, y, parameters){  
      a <- parameters[1]
      b <- parameters[2]
      c <- parameters[3]
      m <- length(y)
      J <- sum(((a-b*X^c)- y)^2)/(2*m)
      return(J) 
    }
    parameters <- optim(par = c(0.75,0.2,-0.5), fn = compCost, 
                        X = X, y = y, method = "L-BFGS-B", 
                        lower = c(0.5,0.0001,-5), 
                        upper = c(1,10,0))[[1]]
    ntrain <- seq(nmin, nmax ,0.01)
    performance_pred <-parameters[1]-parameters[2]*ntrain^parameters[3]
  }
  
  if (metric == "decreasing"){
    compCost<-function(X, y, parameters){  
      a <- parameters[1]
      b <- parameters[2]
      c <- parameters[3]
      m <- length(y)
      J <- sum(((a+b*X^c)- y)^2)/(2*m)
      return(J) 
    }
    parameters <- optim(par = c(min(y),2,-1), fn = compCost, 
                        X = X, y = y, control=list(maxit=100000))[[1]]
    ntrain <- seq(nmin, nmax ,0.01)
    performance_pred <-parameters[1]+parameters[2]*ntrain^parameters[3]
  }
  results <- cbind(ntrain,performance_pred)
  colnames(results) <- c("Training Size", metric) 
  return(results)
}

## function to fit a spline: increasing and concave for classification metrics, and decreasing and convex for regression metrics
SplineFit <- function (ntrain, estimates, constraint, nmin, nmax) {
  training_size <- seq(nmin, nmax ,0.01)
  splinefit <- cobs(x=ntrain, y=estimates, constraint = constraint, lambda = 0, pointwise = NULL, method = "quantile") 
  performance_pred <- predict(splinefit, training_size)[,2]
  return(cbind(training_size,performance_pred))
}

S_CurveFit <- function(ntrain, estimates, nmin, nmax){
  X<-ntrain
  y<-estimates
  compCost<-function(X, y, parameters){
    A <- parameters[1] #left asymptote  
    K <- parameters[2] #right asymptote
    M <- parameters[3]
    B <- parameters[4]
    nu <- parameters[5]
    m <- length(y)
    S <- A + (K-A)/((1+exp(-B*(X-M)))^(1/nu))
    LeastSquares <- sum((S-y)^2)/2*m
    return(LeastSquares)
  }
  parameters <- optim(par = c(max(y),min(y),50,0.1,3),
                      fn = compCost, X = X, y = y,
                      control = list(maxit =10000))[[1]]
  ntrain <- seq(nmin, nmax ,0.01)
  A <- parameters[1] #left asymptote  
  K <- parameters[2] #right asymptote
  M <- parameters[3]
  B <- parameters[4]
  nu <- parameters[5]
  performance_pred <- A + (K-A)/((1+exp(-B*(ntrain-M)))^(1/nu))
  IP <- M - (1/B)*log(nu)
  results <- cbind(ntrain,performance_pred)
  return(list(results,parameters,IP))
}
  
## function to compute optimal training size using bias methodology
ntrain_bias <- function (learning_curve, bias = 0.02, metric = "increasing") {
  # find maximal auc value
  if (metric == "increasing"){
    PointEstimate <- max(learning_curve[,2])
    # selecting bias value
    nopt <- learning_curve[which(learning_curve[,2] > (1-bias)*PointEstimate)[1],] #determines the optimal training set size and corresponding AUC value
  }
  
  if (metric == "decreasing"){
    PointEstimate <- min(learning_curve[,2])
    # selecting bias value
    nopt <- learning_curve[learning_curve[,2] < (1-bias)*PointEstimate,][1,] #determines the optimal training set size and corresponding PMSE value
  }
  return(nopt)
}


# defining function to compute var(AUC)
var_AUC <- function(AUC, n1, n2) {
  q1 = AUC/(2-AUC)
  q2 = 2*AUC^2/(1+AUC)
  var = (AUC*(1-AUC) +(n1-1)*(q1-AUC^2) +(n2-1)*(q2-AUC^2))/(n1*n2)
  return(var)
}

# defining function to compute var(MSE)
var_pMSE <- function(pMSE,n) {
  var = pMSE*sqrt((2/n))
  return(var)
}


# function to find ntrain that minimizes MSE
ntrain_MSE_AUC <- function (learning_curve,  n_proportion) {
  training_sizes <- learning_curve[,1]
  nmax <- max(learning_curve[,1])
  AUC <- max(learning_curve[,2])
  MSE <- c()
  for (i in 1:length(training_sizes)) {
    MSE[i] <- (learning_curve[i,2]-AUC)^2 + var_AUC(learning_curve[i,2], (nmax-learning_curve[i,1])*n_proportion, (nmax-learning_curve[i,1])*(1-n_proportion))
  }
  min_MSE <- which(MSE == min(MSE))
  
  LearningCurve_optim <- c(learning_curve[min_MSE,1], learning_curve[min_MSE,2])
  return(LearningCurve_optim)
}

ntrain_MSE_pMSE <- function (learning_curve) {
  training_sizes <- learning_curve[,1]
  nmax <- max(learning_curve[,1])
  pMSE <- min(learning_curve[,2])
  MSE <- c()
  for (i in 1:length(training_sizes)) {
    MSE[i] <- (learning_curve[i,2]-pMSE)^2 + var_pMSE(learning_curve[i,2], nmax-learning_curve[i,1])
  }
                                                                           
  min_MSE <- which(MSE == min(MSE))
  LearningCurve_optim <- c(learning_curve[min_MSE,1], learning_curve[min_MSE,2])
  return(LearningCurve_optim)
}



