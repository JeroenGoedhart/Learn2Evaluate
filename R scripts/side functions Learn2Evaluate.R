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
    nsam <- length(response)  #total sample size
    ntest <- nsam-ntrain      #test sample size
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
PLCregr <- function(Y,X,nmin=20,nmax=N-10,nev = nev, nrepcv=nrepcv,nfolds=5, nrep=nrep, nested=FALSE,alpha = alpha, fixedsubs=TRUE){ #alpha =1: lasso; alpha = 0: ridge
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
  return(list(resall, nseq=nseq, ntot= nmax))
}

#PLCs for Random Forest
PLCrf <- function(Y,X,nmin=20,nmax=nrow(X)-10,nev = 10,nrep=10, nested=FALSE, fixedsubs=TRUE){
  #X<- t(datStd);nmin=20;nmax=nrow(X)-10;nev = 10; nrepcv=3;nrep=5;learners = c("lasso","ridge","rf");nested <- FALSE
  resp<-factor(Y)
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
      rrfit <- rfsrc(respin ~ ., mtry=sqrt(ncol(X)), var.used="all.trees",ntree=100, data = databoth, importance="none")
      predrf <- predict(rrfit, newdata = data.frame(Xte), importance = F)$predicted[,2]
      res <- c(res,list(cbind(respout=as.numeric(respout)-1,pred=predrf)))
    }
    resall <- c(resall,list(res))
  }
  return(list(resall, nseq=nseq, ntot= nrow(X)))
}

## function to compute Area under the receiver operating curve (AUC)
aucf <- function(resppred, sm=F){
  #resppred <- rfres[[1]][[1]];sm=F
  resp <- resppred[,1]; pred <-resppred[,2]
  ci <- as.numeric(pROC::auc(resp,pred,smooth=sm)) 
  return(ci) 
}

## function to compute confidence bound of AUC using de Long method
ci_auc_delong <- function(resppred,sm=F, alpha=0.05){
  #resp <- plc1[[1]][[1]]$respout;pred <- plc1[[1]][[1]]$pred
  resp <- resppred[,1]; pred <-resppred[,2]
  ci <- as.numeric(pROC::ci.auc(resp,pred,smooth=sm, conf.level=1-2*alpha, method="delong")) 
  return(ci[1]) 
}
## function to fit a power law
PowerlawFit <- function(ntrain, aucs, nmin, nmax){
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
SplineFit <- function (ntrain, aucs, constraint, nmin, nmax) {
  training_size <- seq(nmin, nmax ,0.01)
  con <- rbind(c( 1,min(training_size),0.5), # f(min(x)) >= 0.5
               c(-1,max(training_size),1))   # f(max(x)) <= 1
  
  
  splinefit <- cobs(x=ntrain, y=aucs, constraint = constraint, lambda = 0, pointwise = con, method = "quantile") 
  auc_pred <- predict(splinefit, training_size)[,2]
  return(cbind(training_size,auc_pred))
}

## function to compute optimal training size using bias methodology
ntrain_bias <- function (learning_curve, bias) {
  # find maximal auc value
  AUC_point_estimate <- max(learning_curve[,2])
    # selecting bias value
  nopt <- learning_curve[which(learning_curve[,2] > (1-bias)*AUC_point_estimate)[1],] #determines the optimal training set size and corresponding AUC value
  return(nopt)
}

# defining function to compute var(AUC)
var_AUC <- function(AUC, n1, n2) {
  q1 = AUC/(2-AUC)
  q2 = 2*AUC^2/(1+AUC)
  var = (AUC*(1-AUC) +(n1-1)*(q1-AUC^2) +(n2-1)*(q2-AUC^2))/(n1*n2)
}

# function to find ntrain that minimizes MSE
ntrain_MSE <- function (learning_curve,  n_proportion) {
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



