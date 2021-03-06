#----------------------------------------------------------------------------------#
# Package: picasso                                                                 #
# picasso(): The user interface for picasso()                                      #
# Author: Xingguo Li                                                               #
# Email: <xingguo.leo@gmail.com>                                                   #
# Date: Sep 1st, 2015                                                              #
# Version: 0.4.5                                                                   #
#----------------------------------------------------------------------------------#

oldpicasso <- function(X, 
                    Y, 
                    lambda = NULL,
                    nlambda = NULL,
                    lambda.min.ratio = NULL,
                    lambda.min = NULL,
                    family = "gaussian",
                    method = "l1",
                    alg = "greedy",
                    opt = "naive",
                    gamma = 3,
                    df = NULL,
                    sym = "or",
                    standardize = TRUE,
                    perturb = TRUE,
                    max.act.in = 3,
                    truncation = 1e-2, 
                    prec = 1e-4,
                    max.ite = 1e3,
                    verbose = TRUE)
{
  if(family=="gaussian"){
    if(is.matrix(Y)==FALSE) {
      Y = as.matrix(Y)
    }
    p = ncol(Y)
    if(p==1){
      out = picasso.lasso(X = X, Y = Y, lambda = lambda, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                          lambda.min = lambda.min, method = method, alg = alg, opt = opt, gamma = gamma, df = df, 
                          standardize = standardize, max.act.in = max.act.in, truncation = truncation, prec = prec, 
                          max.ite = max.ite, verbose = verbose)
    }
  }
  if(family=="binomial"){
    if(is.matrix(Y)==FALSE) {
      Y = as.matrix(Y)
    }
    out = picasso.logit(X = X, Y = Y, lambda = lambda, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                        lambda.min = lambda.min, method = method, alg = alg, gamma = gamma, standardize = standardize, 
                        max.act.in = max.act.in, truncation = truncation, prec = prec, max.ite = max.ite, 
                        verbose = verbose)
  }
  if(family=="graph"){
    out = picasso.scio(X = X, lambda = lambda, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                       lambda.min = lambda.min, method = method, alg = alg, opt = opt, gamma = gamma, sym = sym, 
                       max.act.in = max.act.in, truncation = truncation, prec = prec, max.ite = max.ite, 
                       standardize = standardize, perturb = perturb, verbose = verbose)
  }
  out$family = family
  return(out)
}
