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
n = 10000
d = 10000
ratio = 0.01
nlambda = 30

# Linear Regression
set.seed(111)

if( useSimData)
{
  # Simulated data
  print('c=0.25')
  sim_data <- generate_sim(n=n, d=d, c=0.25, seed=112, sigma=0.1)
  
  test_gausnet(sim_data,alg="asp-newton",trialN = trialN,prec=1*1e-9,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="greedy",    trialN = trialN,prec=4*1e-5,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="cyclic",    trialN = trialN,prec=5*1e-5,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="glmnet",    trialN = trialN,prec=7*1e-8,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="gcdnet",    trialN = trialN,prec=3*1e-6,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="ncvreg",    trialN = trialN,prec=2*1e-5,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="fista",     trialN = 1,prec=1*1e-5,ratio=ratio, nlambda = nlambda,fista_it = 100)
  rm(sim_data)
  
  print('c=0.75')
  sim_data <- generate_sim(n=n, d=d, c=0.75, seed=112, sigma=0.1)
  test_gausnet(sim_data,alg="asp-newton",trialN = trialN,prec=3*1e-9,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="greedy",    trialN = trialN,prec=2*1e-5,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="cyclic",    trialN = trialN,prec=2*1e-5,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="glmnet",    trialN = trialN,prec=3*1e-9,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="gcdnet",    trialN = trialN,prec=8*1e-7,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="ncvreg",    trialN = trialN,prec=1*1e-5,ratio=ratio, nlambda = nlambda)
  test_gausnet(sim_data,alg="fista",     trialN = 1,prec=1*1e-5,ratio=ratio, nlambda = nlambda,fista_it = 100)
  rm(sim_data)
}
