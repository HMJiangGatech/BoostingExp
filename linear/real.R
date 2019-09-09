## Notice that the precision have to be fine tunned for getting comparable results
# loading and installing required packages
library(Rcpp)
library(glmnet)
library(picasso)
library(oldpicasso)
library(gcdnet)
library(ncvreg)
library(glmnet)
library(spams)

source("../scripts.R")
sourceCpp("../utils.cpp")

# Experiment parameters
# skip some comparison
useRealData = TRUE
useSimData = FALSE
trialN = 10
# for simulated data set
n = 100
d = 1000
ratio = 0.01
nlambda = 30

# Linear Regression
set.seed(111)

# Real Data
if(useRealData)
{
  cat("===================DrivFace\n")
  load("../data/DrivFace.RData")
  x=as.matrix(x)
  y=as.matrix(y)
  x=scale(x)
  y=scale(y)
  test_gausnet(list(X=x,Y=y),alg="asp-newton",trialN = trialN,prec=3*1e-8,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="greedy",    trialN = trialN,prec=2.5*1e-4,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="cyclic",    trialN = trialN,prec=2*1e-4,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="glmnet",    trialN = trialN,prec=1*1e-7,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="gcdnet",    trialN = trialN,prec=1*1e-8,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="ncvreg",    trialN = trialN,prec=1*1e-4,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="fista",     trialN = 1,     prec=1*1e-3,ratio=ratio, nlambda = nlambda)

  cat("===================GHG\n")
  load("../data/GHG.RData")
  x=scale(x)
  y=scale(y)
  test_gausnet(list(X=x,Y=y),alg="asp-newton",trialN = trialN,prec=1*1e-8,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="greedy",    trialN = trialN,prec=9*1e-5,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="cyclic",    trialN = trialN,prec=8*1e-5,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="glmnet",    trialN = trialN,prec=2*1e-8,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="gcdnet",    trialN = trialN,prec=1*1e-8,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="ncvreg",    trialN = trialN,prec=5*1e-5,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="fista",     trialN = 1,     prec=1*1e-3,ratio=ratio, nlambda = nlambda)

  
  cat("===================riboflavin\n")
  load("../data/riboflavin.RData")
  x=matrix(as.numeric(x),nrow=71,byrow = F)
  y=as.matrix(y)
  x=scale(x)
  y=scale(y)
  test_gausnet(list(X=x,Y=y),alg="asp-newton",trialN = trialN,prec=1*1e-10,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="greedy",    trialN = trialN,prec=2*1e-5,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="cyclic",    trialN = trialN,prec=2*1e-5,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="glmnet",    trialN = trialN,prec=5*1e-10,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="gcdnet",    trialN = trialN,prec=3*1e-11,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="ncvreg",    trialN = trialN,prec=1*1e-11,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="fista",     trialN = 1,     prec=1*1e-3,ratio=ratio, nlambda = nlambda)

}
