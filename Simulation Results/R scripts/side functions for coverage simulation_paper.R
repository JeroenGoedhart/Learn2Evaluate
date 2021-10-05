### functions used for coverage simulation ###

#This file contains all side functions used to compute the coverage of our plc methodology.
# It contains the following functions:
                # subs: used to obtain balanced training subsamples (equal amount of 1's and 0's as response values )
                # PLCregr: this is the function that computes all fits and predictions (so for all different training sizes and repeats) of the learner (lasso or ridge in this case)
                # auc: function that computes the area under the curve for the predicted values obtained by PLCregr
                # ci_auc_delong: the function that calculates the confidence bound by using delongs nonparametric method
                # ci_auc_boot: calculates confidence bound (5%) of auc estimate by using bootstrapping of test set
                # power_law_fit: when the points for the plc are obtained (auc versus ntrain), this function fits a power law model to the data
                # Splinefitting: additionally to the power_law_fit we also fit a smoothing B-spline with additional constraints (monotonic and concave) to the data
                # ntrain_bias: function that finds the minimal training size that is within 2% of the auc predicted at the maximal training size
                # ntrain_MSE: function that minimizes emperical MSE of auc, minimized with respect to training size
                # var_AUC: function that calculates asymptotic variance of AUC (based on gaussian assumption) according to the result of Bamber
                # auc_delaan: calculates auc and its confidence interval by cross validation and 


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
      sub <- sort(c(sample(nev,ntestev), sample(nonev,ntestnonev)))
    }
    return(sub)
  }
  return(lapply(1:nrepeat,subrep))
}

#PLCs for regularized regression (ridge, lasso, Elastic Net)
PLCregr <- function(Y,X,nmin=20,nmax=n_train-10,nev = nev, nrepcv=nrepcv,nfolds=10, nrep=nrep, nested=FALSE,alpha = alpha, fixedsubs=TRUE){ #alpha =1: lasso; alpha = 0: ridge
  #Y<- resp;X<- t(datStd);nmin=20;nmax=nrow(X)-10;nev = 10; nrepcv=3;nrep=5;fixedsubs=TRUE;nested <- FALSE
  nseq <- round(seq(nmin,nmax,length.out=nev))
  print("Sub sample sizes to be evaluated:")
  print(nseq)
  resall <- list()
  #perfridge_mn <-   perflasso_mn <- perfrf_mn <- perfridge2_mn <-   perflasso2_mn <- perfrf2_mn <- c()
  
  #First: create the subsamples; for all learners the same subsamples are used
  for(nseqi in nseq){
    # nseqi <- nseq[2]
    print(paste("Sample size",nseqi))
    if(nested){
      if(nseqi == nseq[1]) subsamp <- Subs(Y=Y,model="logistic",balance=TRUE,ntrain=nseqi,fixedsubs=fixedsubs,nrepeat=nrep) else {
        nremove <- nseqi - nseqprev
        subsamp <- lapply(subsamp,function(ind) {set.seed(24234); el <- length(ind); return(ind[-sample(1:el,nremove)])})
      }
      nseqprev <- nseqi 
    } else {
      subsamp <- Subs(Y=Y,model="logistic",balance=TRUE,ntrain=nseqi,fixedsubs=TRUE,nrepeat=nrep)
    }
    
    print("Cross-validating tuning parameters")
    lams <- c()
    for(k in 1:nrepcv){
      samout <- subsamp[[k]]
      Xtr <- X[-samout,]
      Xte <- X[samout,]
      respin <- Y[-samout]
      respout <- Y[samout]
      cvreg <- cv.glmnet(Xtr,respin,alpha=alpha,family="binomial", nfolds=nfolds, type.measure="deviance")
      lams <- c(lams,cvreg$lambda.min)
      
    }
    lam <- median(lams)
    
    
    res <-list()
    for(k in 1:nrep){
      #k <-1
      print(paste("REPEAT",k))
      samout <- subsamp[[k]]
      Xtr <- X[-samout,]
      Xte <- X[samout,]
      respin <- Y[-samout]
      respout <- Y[samout]
      
      regfit <- glmnet(Xtr,respin,alpha=alpha,family="binomial")
      pred <- as.numeric(predict(regfit, newx=Xte,s=lam, type="response"))
      res<- c(res,list(cbind(respout=respout,pred=pred)))
    }
    resall <- c(resall,list(res))
  }
  return(list(resall, nseq=nseq, ntot= n_train))
}

#PLCs for Random Forest
PLCrf <- function(Y,X,nmin=20,nmax=nrow(X)-10,nev = 10,nrep=10, nested=FALSE, fixedsubs=TRUE){
  #X<- t(datStd);nmin=20;nmax=nrow(X)-10;nev = 10; nrepcv=3;nrep=5;learners = c("lasso","ridge","rf");nested <- FALSE
  resp<-as.factor(Y)
  nseq <- round(seq(nmin,nmax,length.out=nev))
  print("Sub sample sizes to be evaluated:")
  print(nseq)
  resall <- list()
  #perfridge_mn <-   perflasso_mn <- perfrf_mn <- perfridge2_mn <-   perflasso2_mn <- perfrf2_mn <- c()
  
  #First: create the subsamples; for all learners the same subsamples are used
  for(nseqi in nseq){
    # nseqi <- nseq[2]
    print(paste("Sample size",nseqi))
    if(nested){
      if(nseqi == nseq[1]) subsamp <- Subs(Y=resp,model="logistic",balance=TRUE,ntrain=nseqi,fixedsubs=fixedsubs,nrepeat=nrep) else {
        nremove <- nseqi - nseqprev
        subsamp <- lapply(subsamp,function(ind) {set.seed(24234); el <- length(ind); return(ind[-sample(1:el,nremove)])})
      }
      
      nseqprev <- nseqi 
    } else {
      subsamp <- Subs(Y=resp,model="logistic",balance=TRUE,ntrain=nseqi,fixedsubs=TRUE,nrepeat=nrep)
    }
    
    res <-list()
    for(k in 1:nrep){
      #k <-1
      print(paste("REPEAT",k))
      samout <- subsamp[[k]]
      Xtr <- X[-samout,]
      Xte <- X[samout,]
      respin <- resp[-samout]
      respout <- resp[samout]
      
      # rid <- cv.glmnet(Xtr,respin,alpha=0,family="cox", type.measure="deviance")
      # lam <- rid$lambda.min
      databoth <- data.frame(respin,Xtr)
      if(p>=1000) mtryp <- sqrt(p) else mtryp <- p/3
      rrfit <- rfsrc(respin ~ ., mtry=mtryp, var.used="all.trees",ntree=100, data = databoth, importance="none")
      predrf <- predict(rrfit, newdata = data.frame(Xte), importance = F)$predicted[,2]
      res <- c(res,list(cbind(respout=respout,pred=predrf)))
    }
    resall <- c(resall,list(res))
  }
  return(list(resall, nseq=nseq, ntot= nrow(X)))
}




## function to compute area under the curve
aucf <- function(resppred, sm=F){
  #resppred <- rfres[[1]][[1]];sm=F
  resp <- resppred[,1]; pred <-resppred[,2]
  ci <- as.numeric(pROC::auc(resp,pred,smooth=sm)) 
  return(ci) 
}

ci_auc_delong <- function(resppred,sm=F, alpha=0.05){
  #resp <- plc1[[1]][[1]]$respout;pred <- plc1[[1]][[1]]$pred
  resp <- resppred[,1]; pred <-resppred[,2]
  ci <- as.numeric(pROC::ci.auc(resp,pred,smooth=sm, conf.level=1-2*alpha, method="delong")) 
  return(ci[1]) 
}
## function to fit a power law
power_law_fit <- function(aucs,ntrain){
  X<-ntrain
  y<-aucs
  
  compCost<-function(X, y, parameters){  
    a <- parameters[1]
    b <- parameters[2]
    c <- parameters[3]
    m <- length(y)
    J <- sum(((a-b*X^c)- y)^2)/(2*m)
    return(J) 
  }
  parameters <- optim(par = c(0.75,0.2,-0.5), fn = compCost, X = X, y = y, method = "L-BFGS-B", lower = c(0.5,0.0001,-0.99999), upper = c(1,10,0))[[1]]
  ntrain <- seq(nmin, nmax ,0.01)
  auc_pred <-parameters[1]-parameters[2]*ntrain^parameters[3]
  return(cbind(ntrain,auc_pred))
}

## function to fit a monotonic, concave spline
Splinefitting <- function (aucs, ntrain, constraint) {
  training_size <- seq(nmin, nmax ,0.01)
  con <- rbind(c( 1,min(training_size),0.5), # f(min(x)) >= 0.5
               c(-1,max(training_size),1))   # f(max(x)) <= 1
                         
               
  # ridge
  splinefit <- cobs(x=ntrain, y=aucs, constraint = constraint, lambda = 0, pointwise = con, method = "quantile") 
  spline_predictions <- predict(splinefit, training_size)[,2]
  return(cbind(training_size,spline_predictions))
}

## function to compute optimal training size using bias methodology
ntrain_bias <- function (powerlaw_results, spline_results, bias) {
  # find maximal auc value
  max_AUC_powerlaw <- max(powerlaw_results[,2])
  max_AUC_spline <-max(spline_results[,2])
  # selecting bias value
  bias_powerlaw <- powerlaw_results[which(powerlaw_results[,2] > (1-bias)*max_AUC_powerlaw)[1],]
  bias_spline <- spline_results[which(spline_results[,2] > (1-bias)*max_AUC_spline)[1],]
  return(list(bias_spline, bias_powerlaw))
}

# defining function to compute var(AUC)
var_AUC <- function(AUC, n1, n2) {
  q1 = AUC/(2-AUC)
  q2 = 2*AUC^2/(1+AUC)
  var = (AUC*(1-AUC) +(n1-1)*(q1-AUC^2) +(n2-1)*(q2-AUC^2))/(n1*n2)
}

# function to find ntrain that minimizes MSE
ntrain_MSE <- function (powerlaw_results, spline_results,  n_proportion) {
  training_size <- seq(nmin, nmax ,0.01)
  MSE_powerlaw <- c()
  MSE_spline <- c()
  for (i in 1:length(training_size)) {
    MSE_powerlaw[i] <- (powerlaw_results[i,2]-max(powerlaw_results[,2]))^2 + var_AUC(powerlaw_results[i,2], (nmax-powerlaw_results[i,1])*n_proportion, (nmax-powerlaw_results[i,1])*(1-n_proportion))
    MSE_spline[i] <- (spline_results[i,2]-max(spline_results[,2]))^2 + var_AUC(spline_results[i,2], (nmax-spline_results[i,1])*n_proportion, (nmax-spline_results[i,1])*(1-n_proportion))
  }
  min_powerlaw <- which(MSE_powerlaw == min(MSE_powerlaw))
  min_spline <-which(MSE_spline == min(MSE_spline))
  
  powerlaw_optim <- c(powerlaw_results[min_powerlaw,1], powerlaw_results[min_powerlaw,2])
  spline_optim <- c(spline_results[min_spline,1], spline_results[min_spline,2])
  return(list(spline_optim, powerlaw_optim))
}

auc_delaan <- function(X,Y, V, learner){
  folds <- cvFolds(Y=Y, V = V)
  predictions <- list()
  response <- list()
  for (i in 1:V){
    if (learner=="ridge"){
      cvreg <- cv.glmnet(x = X[-folds[[i]],], y = Y[-folds[[i]]], alpha=0,family="binomial", nfolds=10, type.measure="deviance")
      lam_i <- cvreg$lambda.min
      regfit_i <- glmnet(X[-folds[[i]],],Y[-folds[[i]]],alpha=0,family="binomial")
      pred_i <- as.numeric(predict(regfit_i, newx=X[folds[[i]],],s=lam_i, type="response"))
      predictions[[i]] <- pred_i
      response[[i]]<-Y[folds[[i]]]
    }
    if (learner=="lasso"){
      cvreg <- cv.glmnet(x = X[-folds[[i]],], y = Y[-folds[[i]]], alpha=0,family="binomial", nfolds=10, type.measure="deviance")
      lam_i <- cvreg$lambda.min
      regfit_i <- glmnet(X[-folds[[i]],],Y[-folds[[i]]],alpha=0,family="binomial")
      pred_i <- as.numeric(predict(regfit_i, newx=X[folds[[i]],],s=lam_i, type="response"))
      predictions[[i]] <- pred_i
      response[[i]]<-Y[folds[[i]]]
    }
    if (learner =="rf"){
      rffit_i <- randomForest(x=X[-folds[[i]],], y=as.factor(Y[-folds[[i]]]))
      pred_i <- predict(rffit_i, newdata=X[folds[[i]],], type = "prob")[,2]
      predictions[[i]] <- pred_i
      response[[i]]<-Y[folds[[i]]]
    }
    
  }
  out <- ci.cvAUC(predictions=predictions, labels=response, confidence=0.90)
  return(c(out[[1]], out[[3]][1]))
}

## creating folds (cross validation)
cvFolds <- function(Y, V){ #Create CV folds (stratify by outcome)
  Y0 <- split(sample(which(Y==0)), rep(1:V, length=length(which(Y==0))))
  Y1 <- split(sample(which(Y==1)), rep(1:V, length=length(which(Y==1))))
  folds <- vector("list", length=V)
  for (v in seq(V)) {folds[[v]] <- c(Y0[[v]], Y1[[v]])}
  return(folds)
}
