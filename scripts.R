lognet_KKT <- function(data, beta, intercept, lambda){
  return(compute_lognet_KKT(data$X, c(data$Y), c(as.matrix(data$X)%*%beta), c(beta), intercept, lambda))
}

lognet_loss <- function(data, beta, intercept, lambda){
  return(compute_lognet_loss(data$X, c(data$Y), c(data$X%*%beta), c(beta), intercept, lambda))
}

gausnet_loss <- function(data, beta, intercept, lambda){
  return(compute_gausnet_loss(data$X, c(data$Y), c(data$X%*%beta), c(beta), intercept, lambda))
}

sqrtlasso_loss <- function(data, beta, intercept, lambda){
  return(compute_sqrtlasso_loss(data$X, c(data$Y), c(data$X%*%beta), c(beta), intercept, lambda))
}

gausnet_KKT <- function(data, beta, intercept, lambda){
  return(compute_gausnet_KKT(data$X, c(data$Y), c(as.matrix(data$X)%*%beta), c(beta), intercept, lambda))
}


poi_KKT <- function(data, beta, intercept, lambda){
  return(compute_poi_KKT(data$X, c(data$Y), c(as.matrix(data$X)%*%beta), c(beta), intercept, lambda))
}

generate_sim_poi <- function(n, d, c, seed=111) {
  set.seed(seed)
  
  X <- matrix(rnorm(n*d),nrow=n)+rnorm(n,1)*0.5*(c+(c^2+4*c)^0.5)
  s <- 20
  true_beta <- c(runif(s), rep(0, d-s)) 
  lambda <- exp(X%*%true_beta)
  Y <- rpois(n, lambda)
  X <- X[Y<1e4,] # prevent too large Y  test_gausnet(list(X=x,Y=y),skip=skip,trialN = trialN,prec=1
  Y <- Y[Y<1e4] # prevent too large Y
  return(list(X=X, Y=c(Y), true_beta=c(true_beta)))
}

generate_sim_lognet <- function(n, d, c, seed=111) {
  set.seed(seed)
  X <- matrix(rnorm(n*d),nrow=n)+rnorm(n,1)*0.5*(c+(c^2+4*c)^0.5)
  s <- 20
  true_beta <- c(runif(s), rep(0, d-s)) 
  linpred <- X%*%true_beta+rnorm(n)/10
  prob <- exp(linpred)/(1 + exp(linpred))
  runis <- runif(n,0,1)
  Y <- ifelse(runis < prob,1,0)
  return(list(X=X, Y=c(Y), true_beta=c(true_beta)))
}

generate_sim<- function(n, d, c, seed=111, sigma=1) {
  set.seed(seed)
  
  X <- scale(matrix(rnorm(n*d),nrow=n)+rnorm(n,1)*0.5*(c+(c^2+4*c)^0.5))
  
  s <- 18
  true_beta <- rep(0,d)
  true_beta[sample(d,18)] = rep(c(3,2,1.5),6)
  Y <- X%*%true_beta+rnorm(n)*sigma
  Y <- Y
  Y = Y - mean(Y)
  return(list(X=X, Y=c(Y), true_beta=c(true_beta)))
}

pathfista <- function(data, lambdas, tol=1e-6, max_it=100, lostfamily='logistic'){
  nlambda <- length(lambdas)

  
  n <- dim(data$X)[1]
  d <- dim(data$X)[2]
  
  out <- list(beta = matrix(rep(0,d*nlambda), ncol=nlambda),
             b0 = rep(0,n),
             lambda = lambdas, t=c())
  
  bprev <- matrix(rep(0, d+1), ncol=1)
  tt = c()
  ct = 0;
  y <- matrix(data$Y, ncol=1)
  x <- cbind(data$X, matrix(rep(1,n), ncol=1))
  for (i in 1:nlambda){
     t<-system.time(fit <- spams.fistaFlat(y, x, bprev,
                            loss=lostfamily, regul = 'l1', tol=tol, max_it = max_it,
                            intercept = FALSE, # no not regularize last row of beta
                            lambda1 = lambdas[i]))
     out$beta[,i] <- fit[1:d,1]
     out$b0[i] <- fit[d+1,1]
     bprev <- fit
    
     # lambda = lambdas[i]
     # f <- function(beta){ 0.5*norm(x%*%beta - y, "F")^2 }
     # gradf <- function(beta){ cat(class(t(x)%*%(x%*%beta - y)), '\n'); as.vector(t(x)%*%(x%*%beta - y)) }
     # g <- function(beta) { lambda*norm(as.matrix(beta),'1') }
     # proxg <- function(beta, tau) {  sign(beta)*(sapply(c(abs(beta) - drop(tau)*lambda), FUN=function(x) {max(x,0)}))}
     # tau1 <- 10
     # t<-system.time(fit  <- fasta(f,gradf,g,proxg,bprev,tau1) )
     # out$beta[,i] <- sol$x[1:d,1]
     # out$b0[i] <- sol$x[d+1,1]
     # bprev <- sol$x
     # cat("==================")
     ct = ct+t[1]
     tt = c(tt,ct)
  }
  out$t = tt
  return(out)
}

test_lognet <- function(data, nlambda = 100, ratio=0.01, fista_it = 20, trialN = 10, skip=c(),prec=2.0*1e-6, penalty = "lasso"){
  cat("ASP-Newton timing:\n")
  picasso.rtime <- rep(0, trialN) 
  picasso.KKTerr <- rep(0, trialN)
  for (i in 1:trialN){
    method = penalty
    if(penalty=="lasso")
      method = "l1"
    t <- system.time(fitp<-picasso(data$X, data$Y,family="binomial", lambda.min.ratio=ratio, method = method, #alg = "proximal",
                                   standardize=FALSE, verbose=FALSE, prec=prec, nlambda=nlambda))
    picasso.rtime[i] <- t[1]
    err <- rep(0, nlambda)
    for (j in 1:nlambda){
      err[j] <- lognet_KKT(data, fitp$beta[,j], fitp$intercept[j], fitp$lambda[j])
    }
    picasso.KKTerr[i] <- mean(err)
  }
  cat("mean running time: \n")
  print(mean(picasso.rtime))
  cat("standard deviation of running time: \n")
  print(sqrt(var(picasso.rtime)))
  cat("mean KKT error: \n")
  print(mean(picasso.KKTerr))
  cat("last KKT error: \n")
  print(err[nlambda])
  if(!is.null(data$true_beta))
  {  
    cat("estimation error: \n")
    print(norm(as.matrix(fitp$beta[,nlambda] - data$true_beta)))
  }
  cat("objective values: \n")
  print(lognet_loss(data, fitp$beta[,nlambda], fitp$intercept[nlambda], fitp$lambda[nlambda]))
  
  if (!("glmnet" %in% skip)){
    cat("glmnet timing:\n")
    rtime <- rep(0, trialN)
    KKTerr <- rep(0, trialN)

    for (i in 1:trialN){
      t <- system.time(fit<-glmnet(data$X, data$Y, family="binomial",
                                   lambda = fitp$lambda,
                                   standardize=FALSE, thresh=prec))
      rtime[i] <- t[1]
      err <- rep(0, nlambda)
      for (j in 1:nlambda){
        err[j] <- lognet_KKT(data, fit$beta[,j], fit$a0[j], fit$lambda[j])
      }
      KKTerr[i] <- mean(err)
    }
    cat("mean running time: \n")
    print(mean(rtime))
    cat("standard deviation of running time: \n")
    print(sqrt(var(rtime)))
    cat("mean KKT error: \n")
    print(mean(KKTerr))
    cat("last KKT error: \n")
    print(err[nlambda])
    if(!is.null(data$true_beta))
    {  
      cat("estimation error: \n")
      print(norm(as.matrix(fit$beta[,nlambda] - data$true_beta)))
    }
    cat("objective values: \n")
    print(lognet_loss(data, fit$beta[,nlambda], fit$a0[nlambda], fit$lambda[nlambda]))
    
  }
  if (!("ncvreg" %in% skip)){
    cat("ncvreg timing:\n")
    rtime <- rep(0, trialN)
    KKTerr <- rep(0, trialN)
    
    method = penalty
    if(penalty=="mcp")
      method = "MCP"
    if(penalty=="scad")
      method = "SCAD"
    for (i in 1:trialN){
      t <- system.time(fit<-ncvreg(data$X, data$Y, family="binomial", lambda.min = ratio*10,penalty=method,max.iter=100000,
                                   eps=prec,nlambda = nlambda))
      rtime[i] <- t[1]
    }
    err <- rep(0, nlambda)
    for (j in 1:nlambda){
      err[j] <- lognet_KKT(data, fit$beta[-1,j], fit$beta[1,j], fit$lambda[j])
    }
    cat("mean running time: \n")
    print(mean(rtime))
    cat("standard deviation of running time: \n")
    print(sqrt(var(rtime)))
    cat("mean KKT error: \n")
    print(mean(err))
    cat("last KKT error: \n")
    print(err[nlambda])
    if(!is.null(data$true_beta))
    {  
      cat("estimation error: \n")
      print(norm(as.matrix(fit$beta[-1,nlambda] - data$true_beta)))
    }
    cat("objective values: \n")
    print(lognet_loss(data, fit$beta[-1,nlambda], fit$beta[1,nlambda], fit$lambda[nlambda]))
  }
  # 
  # if (!("gcdnet" %in% skip)){
  #   cat("gcdnet timing:\n")
  #   rtime <- rep(0, trialN) 
  #   KKTerr <- rep(0, trialN)
  #   for (i in 1:trialN){
  #     t <- system.time(fit<-gcdnet(data$X, data$Y, method="logit", 
  #                                  lambda = fitp$lambda,
  #                                  standardize=FALSE, eps=2*1e-6))
  #     rtime[i] <- t[1]
  #     err <- rep(0, nlambda)
  #     for (j in 1:nlambda){
  #       err[j] <- lognet_KKT(data, fit$beta[,j], fit$b0[j], fit$lambda[j])
  #     }
  #     KKTerr[i] <- mean(err)
  #   }
  #   cat("mean running time: \n")
  #   print(mean(rtime))
  #   cat("standard deviation of running time: \n")
  #   print(sqrt(var(rtime)))
  #   cat("mean KKT error: \n")
  #   print(mean(KKTerr))
  # }
  # 
  # if (!("fista" %in% skip)){
  #   cat("fista timing:\n")
  #   rtime <- rep(0, trialN) 
  #   KKTerr <- rep(0, trialN)
  #   lambdas <- fitp$lambda
  #   for (i in 1:trialN){
  #     t <- system.time(fit<-pathfista(data, lambdas, max_it=fista_it))
  #     rtime[i] <- t[1]
  #     err <- rep(0, nlambda)
  #     for (j in 1:length(lambdas)){
  #       err[j] <- lognet_KKT(data, fit$beta[,j], fit$b0[j], lambdas[j])
  #     }
  #     KKTerr[i] <- mean(err)
  #   }
  #   cat("mean running time: \n")
  #   print(mean(rtime))
  #   cat("standard deviation of running time: \n")
  #   print(sqrt(var(rtime)))
  #   cat("mean KKT error: \n")
  #   print(mean(KKTerr))
  # }
}


test_gausnet <- function(data, nlambda = 100, ratio=0.01, fista_it = 20, trialN = 10, algo="",prec=2.0*1e-6, max.iter=1000, penalty = "lasso"){
  cat("=============")
  if (algo %in% c("asp-newton","greedy","cyclic") ){
    cat(algo, "timing:\n")
    picasso.rtime <- rep(0, trialN) 
    picasso.KKTerr <- rep(0, trialN)
    for (i in 1:trialN){
      method = penalty
      if(penalty=="lasso")
        method = "l1"
      if(algo=="asp-newton")
        t <- system.time(fitp<-picasso(data$X, data$Y,family="gaussian", lambda.min.ratio=ratio, method=method, #alg = "proximal",
                                       standardize=FALSE, verbose=FALSE, prec=prec, nlambda=nlambda, max.ite=max.iter))
      else if(algo=="cyclic")
        t <- system.time(fitp<-oldpicasso(data$X, data$Y,family="gaussian", lambda.min.ratio=ratio, alg = "cyclic",
                                       standardize=FALSE, verbose=FALSE, prec=prec, nlambda=nlambda, max.ite=max.iter))
      else if(algo=="greedy")
        t <- system.time(fitp<-oldpicasso(data$X, data$Y,family="gaussian", lambda.min.ratio=ratio, alg = "greedy",
                                       standardize=FALSE, verbose=FALSE, prec=prec, nlambda=nlambda, max.ite=max.iter))
      picasso.rtime[i] <- t[1]
    }
    
    err <- rep(0, nlambda)
    for (j in 1:nlambda){
      system.time(err[j] <- gausnet_KKT(data, fitp$beta[,j], fitp$intercept[j], fitp$lambda[j]))
      # cat("sparsity: ")
      # print( sum(fitp$beta[,j]==0)/nrow(fitp$beta) )
    }
    
    cat("mean running time: \n")
    print(mean(picasso.rtime))
    cat("standard deviation of running time: \n")
    print(sqrt(var(picasso.rtime)))
    cat("mean KKT error: \n")
    print(mean(err))
    cat("last KKT error: \n")
    print(err[nlambda])
    if(!is.null(data$true_beta))
    {  
      cat("estimation error: \n")
      print(norm(as.matrix(fitp$beta[,nlambda] - data$true_beta)))
    }
    cat("objective values: \n")
    print(gausnet_loss(data, fitp$beta[,nlambda], fitp$intercept[nlambda], fitp$lambda[nlambda]))
  }
  
  
  if (algo == "glmnet"){
    cat("glmnet timing:\n")
    rtime <- rep(0, trialN)
    KKTerr <- rep(0, trialN)
    fitp<-picasso(data$X, data$Y,family="gaussian", lambda.min.ratio=ratio, #alg = "proximal",
                  standardize=FALSE, verbose=FALSE, prec=prec, nlambda=nlambda, max.ite=1)
    
    for (i in 1:trialN){
      t <- system.time(fit<-glmnet(data$X, data$Y, family="gaussian", lambda=fitp$lambdas,
                                   standardize=FALSE, thresh=prec))
      rtime[i] <- t[1]
    }
    err <- rep(0, nlambda)
    for (j in 1:nlambda){
      err[j] <- gausnet_KKT(data, fit$beta[,j], fit$a0[j], fit$lambda[j])
    }
    cat("mean running time: \n")
    print(mean(rtime))
    cat("standard deviation of running time: \n")
    print(sqrt(var(rtime)))
    cat("mean KKT error: \n")
    print(mean(err))
    cat("last KKT error: \n")
    print(err[nlambda])
    if(!is.null(data$true_beta))
    {  
      cat("estimation error: \n")
      print(norm(as.matrix(fit$beta[,nlambda] - data$true_beta)))
    }
    cat("objective values: \n")
    print(gausnet_loss(data, fit$beta[,nlambda], fit$a0[nlambda], fit$lambda[nlambda]))
  }
  
  if (algo == "ncvreg"){
    cat("ncvreg timing:\n")
    rtime <- rep(0, trialN)
    KKTerr <- rep(0, trialN)
    
    method = penalty
    if(penalty=="mcp")
      method = "MCP"
    if(penalty=="scad")
      method = "SCAD"
    for (i in 1:trialN){
      t <- system.time(fit<-ncvreg(data$X, data$Y, family="gaussian", lambda.min = ratio,penalty=method,max.iter=100000,
                                   eps=prec,nlambda = nlambda,returnX=TRUE))
      rtime[i] <- t[1]
    }
    err <- rep(0, nlambda)
    for (j in 1:nlambda){
      err[j] <- gausnet_KKT(data, fit$beta[-1,j], fit$beta[1,j], fit$lambda[j])
    }
    cat("mean running time: \n")
    print(mean(rtime))
    cat("standard deviation of running time: \n")
    print(sqrt(var(rtime)))
    cat("mean KKT error: \n")
    print(mean(err))
    cat("last KKT error: \n")
    print(err[nlambda])
    cat("objective values: \n")
    print(gausnet_loss(data, fit$beta[-1,nlambda], fit$beta[1,nlambda], fit$lambda[nlambda]))
  }

  if (algo == "gcdnet"){
    cat("gcdnet timing:\n")
    rtime <- rep(0, trialN)
    KKTerr <- rep(0, trialN)
    for (i in 1:trialN){
      t <- system.time(fit<-gcdnet(data$X, data$Y, method="ls",lambda.factor=ratio,
                                   ,nlambda = nlambda,
                                   standardize=FALSE, eps=prec))
      rtime[i] <- t[1]
    }
    err <- rep(0, nlambda)
    for (j in 1:nlambda){
      err[j] <- gausnet_KKT(data, fit$beta[,j], fit$b0[j], fit$lambda[j])
    }
    cat("mean running time: \n")
    print(mean(rtime))
    cat("standard deviation of running time: \n")
    print(sqrt(var(rtime)))
    cat("mean KKT error: \n")
    print(mean(err))
    cat("last KKT error: \n")
    print(err[nlambda])
  }

  if (algo == "fista"){
    fitp<-picasso(data$X, data$Y,family="gaussian", lambda.min.ratio=ratio, #alg = "proximal",
                                         standardize=FALSE, verbose=FALSE, prec=prec, nlambda=nlambda, max.ite=20)
    cat("fista timing:\n")
    rtime <- rep(0, trialN)
    KKTerr <- rep(0, trialN)
    lambdas <- fitp$lambda
    for (i in 1:trialN){
      t <- system.time(fit<-pathfista(data, fitp$lambda, max_it=fista_it, lostfamily='square',tol=prec))
      rtime[i] <- t[1]
    }
    err <- rep(0, nlambda)
    for (j in 1:length(lambdas)){
      err[j] <- gausnet_KKT(data, fit$beta[,j], fit$b0[j], lambdas[j])
    }
    KKTerr[i] <- mean(err)
    cat("mean running time: \n")
    print(mean(rtime))
    cat("standard deviation of running time: \n")
    print(sqrt(var(rtime)))
    cat("mean KKT error: \n")
    print(mean(KKTerr))
    cat("last KKT error: \n")
    print(err[nlambda])
    if(!is.null(data$true_beta))
    {  
      cat("estimation error: \n")
      print(norm(as.matrix(fitp$beta[,nlambda] - data$true_beta)))
    }
    cat("objective values: \n")
    print(gausnet_loss(data, fitp$beta[,nlambda], fitp$intercept[nlambda], fitp$lambda[nlambda]))
  }
}


test_poi <- function(data, nlambda = 100, ratio=0.01, fista_it = 20, trialN = 10, skip=c(), prec=2.0*1e-6){
  cat("ASP-Newton timing:\n")
  picasso.rtime <- rep(0, trialN) 
  picasso.KKTerr <- rep(0, trialN)
  for (i in 1:trialN){
    t <- system.time(fitp<-picasso(data$X, data$Y,family="poisson", lambda.min.ratio=ratio, #alg = "proximal",
                                   standardize=FALSE, verbose=FALSE, prec=prec, nlambda=nlambda))
    picasso.rtime[i] <- t[1]
    err <- rep(0, nlambda)
    for (j in 1:nlambda){
      err[j] <- poi_KKT(data, fitp$beta[,j], fitp$intercept[j], fitp$lambda[j])
    }
    picasso.KKTerr[i] <- mean(err)
  }
  cat("mean running time: \n")
  print(mean(picasso.rtime))
  cat("standard deviation of running time: \n")
  print(sqrt(var(picasso.rtime)))
  cat("mean KKT error: \n")
  print(mean(picasso.KKTerr))
  cat("last KKT error: \n")
  print(err[nlambda])
  if(!is.null(data$true_beta))
  {  
    cat("estimation error: \n")
    print(norm(as.matrix(fitp$beta[,nlambda] - data$true_beta)))
  }
  
  
  if (!("glmnet" %in% skip)){
    cat("glmnet timing:\n")
    rtime <- rep(0, trialN)
    KKTerr <- rep(0, trialN)
    
    for (i in 1:trialN){
      t <- system.time(fit<-glmnet(data$X, data$Y, family="poisson",
                                   lambda = fitp$lambda,
                                   standardize=FALSE, thresh=prec))
      rtime[i] <- t[1]
      err <- rep(0, nlambda)
      for (j in 1:nlambda){
        err[j] <- poi_KKT(data, fit$beta[,j], fit$a0[j], fit$lambda[j])
      }
      KKTerr[i] <- mean(err)
    }
    cat("mean running time: \n")
    print(mean(rtime))
    cat("standard deviation of running time: \n")
    print(sqrt(var(rtime)))
    cat("mean KKT error: \n")
    print(mean(KKTerr))
    cat("last KKT error: \n")
    print(err[nlambda])
    if(!is.null(data$true_beta))
    {  
      cat("estimation error: \n")
      print(norm(as.matrix(fit$beta[,nlambda] - data$true_beta)))
    }
    
  }
}

test_sqrtlasso <- function(data, nlambda = 100, ratio=0.01, fista_it = 20, trialN = 10, skip=c(),prec=2.0*1e-6, max.iter=1000){
  
  cat("ASP-Newton timing:\n")
  picasso.rtime <- rep(0, trialN) 
  picasso.KKTerr <- rep(0, trialN)
  for (i in 1:trialN){
    t <- system.time(fitp<-picasso(data$X, data$Y,family="sqrtlasso", lambda.min.ratio=ratio, 
                                   standardize=FALSE, verbose=FALSE, prec=prec, nlambda=nlambda, max.ite=max.iter))
    picasso.rtime[i] <- t[1]
  }
  
  cat("mean running time: \n")
  print(mean(picasso.rtime))
  cat("standard deviation of running time: \n")
  print(sqrt(var(picasso.rtime)))
  cat("objective values: \n")
  print(sqrtlasso_loss(data, fitp$beta[,nlambda], fitp$intercept[nlambda], fitp$lambda[nlambda]))  
  
  
  cat("scalreg timing:\n")
  rtime <- rep(0, trialN) 
  for (i in 1:trialN){
    t <- system.time(fit<-scalreg(data$X, data$Y,lam0=fitp$lambda[nlambda]))
    rtime[i] <- t
  }
  
  cat("mean running time: \n")
  print(mean(rtime))
  cat("standard deviation of running time: \n")
  print(sqrt(var(rtime)))
  cat("objective values: \n")
  print(mean(fit$residuals^2)^0.5 + fitp$lambda[nlambda]*sum(abs(fit$coefficients)))
  
  
  cat("flare timing:\n")
  rtime <- rep(0, trialN) 
  for (i in 1:trialN){
    t <- system.time(fit<-slim(data$X, data$Y,nlambda = nlambda, lambda.min.ratio=ratio))
    rtime[i] <- t
  }
  
  cat("mean running time: \n")
  print(mean(rtime))
  cat("standard deviation of running time: \n")
  print(sqrt(var(rtime)))
  cat("objective values: \n")
  print(sqrtlasso_loss(data, fit$beta[,nlambda], fit$intercept[nlambda], fit$lambda[nlambda]))
}