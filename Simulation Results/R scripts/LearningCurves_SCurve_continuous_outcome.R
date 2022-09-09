gc()

### Simulation Continuous Outcome -S shaped curve ###
setwd("C:/Users/VNOB-0732/Desktop/R files/Learn2Evaluate/Learn2Evaluate")
library(glmnet) 
library(mvtnorm)
library(randomForest)

### All side functions used ###
###############################

PLCregr <- function(Y,X, nfolds=5, alfa=1, nrepcv=10,nrep=50, subs = Subsamp,
                    family){ #alpha =1: lasso; alpha = 0: ridge
  lams <- c() # cv lambda values
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
  res_oracle <- list()
  
  for(k in 1:nrep){
    #k <-1
    samout <- subs[[k]]
    Xtr <- X[-samout,]
    Xte <- X[samout,]
    respin <- Y[-samout]
    respout <- Y[samout]
    regfit <- glmnet(Xtr,respin,alpha=alfa,family=family)
    pred <- as.numeric(predict(regfit, newx=Xte,s=lam, type="response"))
    pred_oracle <- as.numeric(predict(regfit, newx=X_test,s=lam, type="response"))
    res[[k]] <- cbind(respout=respout,pred=pred)
    res_oracle[[k]] <- cbind(Y_test=Y_test,pred=pred_oracle)
  }
  return(list(res,res_oracle))
}
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

pMSE <- function(resppred){
  resp <- resppred[,1]; pred <- resppred[,2]
  pMSE <- mean((resp-pred)^2)
  return(pMSE)
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
  return(results)
}

SplineFit <- function (ntrain, estimates, constraint, nmin, nmax) {
  training_size <- seq(nmin, nmax ,0.01)
  splinefit <- cobs(x=ntrain, y=estimates, constraint = constraint, lambda = 0, pointwise = NULL, method = "quantile") 
  performance_pred <- predict(splinefit, training_size)[,2]
  #return(cbind(training_size,performance_pred))
  return(performance_pred)
}

####################################################

# setting up the simulation
n_train = 100; #number of training samples, will be generated Nsim times  
p = 500; #number of covariates
nmin = 10; nmax=n_train; nev=20;nrep=50;nrepcv=5; ## parameters to generate PLC

# simulating X and Y
betas <- c(rep(0.6,5),rep(0,p-5))
sigma <- 0.5
n_test <- 5000
X_test <- rmvnorm(n_test, mean = rep(0,p), sigma = diag(1,p,p), method = "chol")
Y_test <- as.numeric(X_test %*% betas) + rnorm(n_test, mean = 0, sd = sigma)

nseq <- round(seq(nmin,n_train-10,length.out=nev)) #subsamples to evaluate, homogeneously spread over sample size regime
measure <- "PMSE"
if (measure == "PMSE") {
  model = "linear"
  distr = "gaussian"
  balance = F
  metric = "decreasing"
  constraints = c("decrease","convex")
}
load("seeds.Rdata") # used seeds to simulate data

results <- list() # storage of results
for (i in 1:length(seeds)) {
  print(paste("iteration",i, sep = " "))
  set.seed(seeds[i])
  X_train <- rmvnorm(n_train, mean = rep(0,p), sigma = diag(1,p,p), method = "chol")
  Y_train <- as.numeric(X_train %*% betas) + rnorm(n_train, mean = 0, sd = sigma)
  resall <- list() #empty list to store all predictions and corresponding responses
  resall_oracle <- list()
  for (j in 1:nev) {
    Subsamp <- Subs(Y=Y_train,model=model,balance=balance,ntrain=nseq[j],fixedsubs=TRUE,nrepeat=nrep)
    plcres0 <-PLCregr(Y=Y_train,X=X_train, nfolds=5,nrepcv = nrepcv, alfa=1, nrep = nrep, subs = Subsamp,
                      family = distr)
    resall[[j]] <- plcres0[[1]]
    resall_oracle[[j]]<-plcres0[[2]]
    names(resall)[j] <- as.character(nseq[j])
    names(resall_oracle)[j] <- as.character(nseq[j])
  }
  estimates <- lapply(resall,function(ex) unlist(lapply(ex,pMSE)))
  estimatesmn <- unlist(lapply(estimates,mean))
  estimates_oracle <- lapply(resall_oracle,function(ex) unlist(lapply(ex,pMSE)))
  estimatesmn_oracle <- unlist(lapply(estimates_oracle,mean))
  LearnCurve <- cbind(training_size = nseq, pMSE = estimatesmn, pMSE_oracle = estimatesmn_oracle)
  results[[i]] <- LearnCurve
}

filenm <- "LearningTrajectories_estimate_vs_oracle_0.6_0.5_100_500.Rdata"
save(results, file = filenm)
load("LearningTrajectories_estimate_vs_oracle.Rdata")


#### extracting results
n_train = 100
nmin = 10; nmax=n_train; nev=20;
nseq <- round(seq(nmin,n_train-10,length.out=nev))

# computes for each considered sample size the difference between
# estimates and oracle predictive performance

results_combined <- do.call(rbind.data.frame, results)
diff <- list()
for (i in 1:length(nseq)) {
  ids<-which(results_combined$training_size==nseq[i])
  samples <- results_combined[ids,]
  diff[[i]] <- samples[,2]-samples[,3]
}

pdf("Boxplot.pdf", width=7,height=6)
boxplot(diff, outline = F, ylim= c(-2,2), ylab = "Error", xlab="Training Size", names = as.character(nseq), col="#009E73", show.names=T)
ticks<-c(-2,0,2)
axis(2,at=ticks,labels=ticks)
dev.off()

### Smoothing the results ###
#############################

### 1. Scurves ###
ntrain <- seq(nmin, nmax ,0.01)
length(ntrain)
results[[1]]
Scurves_estimates <- lapply(results, function (x) S_CurveFit(ntrain = x[,1], estimates = x[,2], nmin = nmin, nmax = nmax))
Scurves_oracle <- lapply(results, function (x) S_CurveFit(ntrain = x[,1], estimates = x[,3], nmin = nmin, nmax = nmax))

Scurves_estimates <- do.call(cbind, Scurves_estimates)
Scurves_oracle <- do.call(cbind, Scurves_oracle)

mean <- rowMeans(Scurves_estimates)
quant1 <- apply(Scurves_estimates, 1, quantile, probs = c(0.95),  na.rm = TRUE)
quant2 <- apply(Scurves_estimates, 1, quantile, probs = c(0.05),  na.rm = TRUE) 

mean_oracle <- rowMeans(Scurves_oracle)
quant1_oracle <- apply(Scurves_oracle, 1, quantile, probs = c(0.95),  na.rm = TRUE)
quant2_oracle <- apply(Scurves_oracle, 1, quantile, probs = c(0.05),  na.rm = TRUE) 

# plotting
pdf("Scurves_estimates_oracle_0.6_0.5.pdf", width=7,height=6)
plot(ntrain, mean, type = "l", lwd =2, xlab = "Training Size", ylab = "PMSE",font.lab=2, font=2, ylim = c(0,4))
lines(ntrain, quant1, lty = 2)
lines(ntrain, quant2, lty = 2)
lines(ntrain, mean_oracle, col = "red")
#lines(ntrain, quant1_oracle, lty = 2, col = "red")
#lines(ntrain, quant2_oracle, lty = 2, col = "red")
legend("topright", legend = c("Estimates", "Oracle"), col = c("black","red"), lty =1, bty = "n", cex = 0.8)
dev.off()

### plotting individual curves ###
set.seed(11)
ids <-sample(1000, 3, replace = F)
curve_estimate_1 <- Scurves_estimates[[ids[3]]]
curve_oracle_1 <- Scurves_oracle[[ids[3]]]
pdf("Ind3_curves_estimates_oracle_0.6_0.5.pdf", width=7,height=6)
plot(ntrain, curve_estimate_1[,2], type = "l", lwd =1, xlab = "Training Size", ylab = "PMSE",font.lab=2, font=2, ylim = c(0,4))

lines(ntrain, curve_oracle_1[,2], col = "red", lwd =2)
#lines(ntrain, quant1_oracle, lty = 2, col = "red")
#lines(ntrain, quant2_oracle, lty = 2, col = "red")
legend("topright", legend = c("Estimates", "Oracle"), col = c("black","red"), lty =1, bty = "n", cex = 0.8)
dev.off()









### 2. Constrained smoothing splines ###
library(cobs)
Spline_estimates <- lapply(results, function (x) SplineFit(ntrain = x[,1], estimates = x[,2], nmin = nmin, nmax = nmax, constraint = "decrease"))
Spline_oracle <- lapply(results, function (x) SplineFit(ntrain = x[,1], estimates = x[,3], nmin = nmin, nmax = nmax,constraint = "decrease"))
Spline_estimates <- do.call(cbind, Spline_estimates)
Spline_oracle <- do.call(cbind, Spline_oracle)
ntrain <- seq(nmin, nmax ,0.01)

mean <- rowMeans(Spline_estimates)
quant1 <- apply(Spline_estimates, 1, quantile, probs = c(0.975),  na.rm = TRUE)
quant2 <- apply(Spline_estimates, 1, quantile, probs = c(0.025),  na.rm = TRUE) 

mean_oracle <- rowMeans(Spline_oracle)
quant1_oracle <- apply(Spline_oracle, 1, quantile, probs = c(0.975),  na.rm = TRUE)
quant2_oracle <- apply(Spline_oracle, 1, quantile, probs = c(0.025),  na.rm = TRUE) 

# plotting
pdf("Splines_oracle_vs_estimates_0,6_0.5.pdf", width=7,height=6)
plot(ntrain, mean, type = "l", lwd =2, xlab = "Training Size", ylab = "PMSE",font.lab=2, font=2, ylim = c(0,4))
lines(ntrain, quant1, lty = 2)
lines(ntrain, quant2, lty = 2)
lines(ntrain, mean_oracle, col = "red")
#lines(ntrain, quant1_oracle, lty = 2, col = "red")
#lines(ntrain, quant2_oracle, lty = 2, col = "red")
legend("topright", legend = c("Estimates", "Oracle"), col = c("black","red"), lty =1, bty = "n", cex = 0.8)
dev.off()


legend()
LearnCurve_oracle <-cbind(training_size = nseq, pMSE = estimatesmn_oracle)
LearnCurve <- results[[1]]
plot(LearnCurve[,1],LearnCurve[,2], xlim=c(nmin,n_train), ylim = c(0,var(Y_train)+0.5), xlab="Sample Size", ylab="PMSE", pch = 1, cex =0.7)
points(LearnCurve[,1],LearnCurve[,3],pch = 2, cex =0.7, col = "red")
points(LearnCurve_oracle[,1],LearnCurve_oracle[,2], pch = 2, col="#395D9CFF", cex = 0.7,xlab="Sample Size", ylab="PMSE")

a <- S_CurveFit(ntrain = LearnCurve[,1],
                estimates = LearnCurve[,2],
                nmin = nmin,nmax = nmax)
lines(a[[1]][,1],a[[1]][,2],lwd=1)
points(LearnCurve_oracle[,1],LearnCurve_oracle[,2], pch = 3, col="red", cex = 0.7)




PLCregr_oracle <- function(Y,X, nfolds=5, alfa=1, nrepcv=10,nrep=50, subs = Subsamp,
                           family){ #alpha =1: lasso; alpha = 0: ridge
  lams <- c() # cv lambda values
  for(k in 1:nrepcv){
    samout <- subs[[k]]
    Xtr <- X[-samout,]
    respin <- Y[-samout]
    cvreg <- suppressWarnings(cv.glmnet(Xtr,respin,alpha=alfa,family=family, nfolds=nfolds, type.measure="deviance"))
    lams <- c(lams,cvreg$lambda.min)
  }
  
  lam <- median(lams)
  res_oracle <- list()
  
  for(k in 1:nrep){
    #k <-1
    samout <- subs[[k]]
    Xtr <- X[-samout,]
    respin <- Y[-samout]
    regfit <- glmnet(Xtr,respin,alpha=alfa,family=family)
    pred_oracle <- as.numeric(predict(regfit, newx=X_test,s=lam, type="response"))
    res_oracle[[k]] <- cbind(Y_test=Y_test,pred=pred_oracle)
  }
  return(res_oracle)
}

