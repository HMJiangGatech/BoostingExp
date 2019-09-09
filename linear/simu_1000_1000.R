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
useRealData = FALSE
useSimData = TRUE
trialN = 10
# for simulated data set
n = 1000
d = 1000
ratio = 0.01
nlambda = 30

# Linear Regression
set.seed(111)

if( useSimData)
{
  # Simulated data
  print('c=0.25')
  sim_data <- generate_sim(n=n, d=d, c=0.25, seed=112)
  
  test_gausnet(sim_data,alg="asp-newton",trialN = trialN,prec=1*1e-8,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="greedy",    trialN = trialN,prec=3*1e-5,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="cyclic",    trialN = trialN,prec=5*1e-5,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="glmnet",    trialN = trialN,prec=2*1e-8,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="gcdnet",    trialN = trialN,prec=3*1e-6,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="ncvreg",    trialN = trialN,prec=1*1e-5,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="fista",     trialN = 1,prec=1*1e-5,ratio=ratio, nlambda = nlambda,fista_it = 100)
  rm(sim_data)
  
  print('c=0.75')
  sim_data <- generate_sim(n=n, d=d, c=0.75, seed=112)
  test_gausnet(sim_data,alg="asp-newton",trialN = trialN,prec=1*1e-9,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="greedy",    trialN = trialN,prec=1.5*1e-5,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="cyclic",    trialN = trialN,prec=5*1e-6,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="glmnet",    trialN = trialN,prec=3*1e-9,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="gcdnet",    trialN = trialN,prec=4*1e-7,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="ncvreg",    trialN = trialN,prec=1*1e-5,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="fista",     trialN = 1,prec=1*1e-5,ratio=ratio, nlambda = nlambda,fista_it = 100)
  rm(sim_data)
}

# Real Data
if(useRealData)
{
  load("DrivFace.RData")
  x=as.matrix(x)
  y=as.matrix(y)
  x=scale(x)
  y=scale(y)
  test_gausnet(list(X=x,Y=y),alg="asp-newton",trialN = trialN,prec=2*1e-6,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="greedy",    trialN = trialN,prec=5*1e-4,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="cyclic",    trialN = trialN,prec=1.5*1e-4,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="glmnet",    trialN = trialN,prec=2*1e-6,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="gcdnet",    trialN = trialN,prec=4*1e-7,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="ncvreg",    trialN = trialN,prec=1*1e-3,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="fista",     trialN = 1,prec=1*1e-3,ratio=ratio, nlambda = nlambda)


  load("GHG.RData")
  x=scale(x)
  y=scale(y)
  test_gausnet(list(X=x,Y=y),alg="asp-newton",trialN = trialN,prec=1*1e-6,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="greedy",    trialN = trialN,prec=5*1e-4,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="cyclic",    trialN = trialN,prec=5*1e-5,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="glmnet",    trialN = trialN,prec=1*1e-6,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="gcdnet",    trialN = trialN,prec=3*1e-7,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="ncvreg",    trialN = trialN,prec=6*1e-4,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="fista",     trialN = 1,prec=1*1e-3,ratio=ratio, nlambda = nlambda)


  load("riboflavin.RData")
  x=matrix(as.numeric(x),nrow=71,byrow = F)
  y=as.matrix(y)
  x=scale(x)
  y=scale(y)
  test_gausnet(list(X=x,Y=y),alg="asp-newton",trialN = trialN,prec=3*1e-8,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="greedy",    trialN = trialN,prec=2*1e-5,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="cyclic",    trialN = trialN,prec=1*1e-6,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="glmnet",    trialN = trialN,prec=3*1e-8,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="gcdnet",    trialN = trialN,prec=5*1e-9,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="ncvreg",    trialN = trialN,prec=1*1e-5,ratio=ratio, nlambda = nlambda)
  test_gausnet(list(X=x,Y=y),alg="fista",     trialN = 1,prec=1*1e-1,ratio=ratio, nlambda = nlambda)

  test_gausnet(DrivFace,skip=skip)
}
