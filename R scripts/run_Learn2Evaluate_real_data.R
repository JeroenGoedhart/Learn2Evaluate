############################################################################
###############  Application of Learn2Evaluate   ###########################
############################################################################


load("Bloodplatelet_RNAseq.Rdata") #loads object dataSqrt into environment

### Short Data Description ###
#The data consists of RNAseq data obtained from bloodplatelets of in total 274 patients.
# 55 are control (=HC group) and the remaining patients have 1 of the in total six tumor types:
      
            ## Breast          CRC (colecteral)          GBM (Glioblastoma)          HC (control)       
            ## n=40            n=41                      n=40                        n=55 

            ## Liver          NSCLC (Non-small-cell lung)          Pancreas (Glioblastoma)          
            ## n=14           n=60                                 n=35                        

# We exclude liver, because the sample size is too small (n=14)
# Here we apply Learn2Evaluate to this dataset. We will only consider binary classification.
# Thus, we may perform 15 binary classifications in total:
  #compareGroups = c("Breast,"CRC")
  #compareGroups = c("Breast,"HC")
  #compareGroups = c("Breast,"GBM")
  #compareGroups = c("Breast,"Pancreas")
  compareGroups = c("Breast","NSCLC")
  #compareGroups = c("NSCLC","Pancreas")
  #compareGroups = c("NSCLC","HC")
  #compareGroups = c("NSCLC","GBM")
  #compareGroups = c("NSCLC","CRC")
  #compareGroups = c("CRC,"Pancreas")
  #compareGroups = c("CRC,"HC")
  #compareGroups = c("CRC,"GBM")
  #compareGroups = c("GBM,"HC")
  #compareGroups = c("GBM,"Pancreas")
  #compareGroups = c("HC,"Pancreas")

#### Data Preprocessing ####
############################

#Prepare subsetting data columns, so that resulting data set contains samples of groups to be compared
g1 = compareGroups[1]; g2 = compareGroups[2]

id1 = which(group==g1) #index of group 1
id2 = which(group==g2) #index of group 2

#subset response vector as well [group indicator]
resp = factor(group[c(id1,id2)])

#subset the data
datSqrt2 = dataSqrt[,c(id1,id2)]

#check
colnames(datSqrt2)

#filtering of features based on sd (removingfeatures whose sd=0, no variablity across samples)
sds = apply(datSqrt2,1,sd)
id.del = which(sds==0)
datSqrt2 = datSqrt2[-id.del,]

#standardize data; important for regularized regression approaches
datStd = t(apply(datSqrt2,1,function(x){(x-mean(x))/sd(x)}))    
dim(datStd)


####  Application ##########
############################

#load required functions into environment
source("side functions Learn2Evaluate.R")
source("Learn2Evaluate.R")
X <-as.data.frame(t(datStd)) #Learn2Evaluate requires X to be a dataframe. (transpose is taken to make sure that rows represent samples)

## function to compute predictive performance estimates with Learn2Evaluate
L2E <- Learn2Evaluate(X = X, Y=resp, nmin = 20, nev = 10, nrep = 20, nrepcv = 5, 
                      curve = "IPL", method = "MSE", learner = "rf", Plot = T)
