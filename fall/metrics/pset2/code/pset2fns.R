#####################################################
#               Metrics Problem Set 2
#                     Tanya Rajan
# Description: This code defines functions for both
# of the computational questions on this problem set.
# The main functions of interest are the OLS function
# reg and the IV function iv. 
#####################################################

########### Helper functions ##############

# evaluating arguments in text
callit <- function(...){
  arg <- parse(text=paste0(..., collapse=""))
  eval(parse(text=arg))
}

# pulling apart matrices
pull <- function(m, dim){
  if (dim == 1){
    return(lapply(1:nrow(m), function(i) m[i,])) # extract rows
  }
  if (dim == 2){
    return(lapply(1:ncol(m), function(j) m[,j])) # isolates cols
  }
}

########### Regression Functions ##############

### leverage calculation for jackknife ###
levs <- function(h, pi, z, x){
  h<-as.numeric(h)
  return((z%*%pi - h*x)/(1-h))
}

### OLS regression function ###

reg <- function(Y, X, se=""){
  # dimensions
  n <- nrow(as.matrix(X))
  k <- ncol(as.matrix(X))
  
  # solving for beta
  Xols <- cbind(1,X)
  XXinv <-  solve(t(Xols)%*%Xols, tol=1e-40)
  beta <- XXinv%*%(t(Xols)%*%Y)
  
  # standard errors
  u2 <- ((Y- Xols%*%beta)^2)
  covmat <- mean(u2)*XXinv
  stderr <- sqrt(diag(covmat))
  covered<- (beta[2] + 1.96*stderr[2] >= 1 & beta[2] - 1.96*stderr[2] <= 1)
  
  # heteroskedastic se
  if (se == "HC1"){
    covmat<-(n/(n-k))*XXinv%*%t(Xols)%*%diag(as.vector(u2))%*%Xols%*%XXinv
    stderr <- sqrt(diag(covmat))
  }
  
  # returns
  return(list(b=beta, se=stderr, cov=covered, covmat=covmat))
}



### 2SLS Function ###
iv <- function(Y, X, Z, est="", se=""){
  # dimensions
  n <- nrow(as.matrix(X))
  k <- ncol(as.matrix(X))
  
  # first stage
  Zols <- cbind(1, Z)
  Xols <- cbind(1, X)
  ZZinv <- solve(t(Zols)%*%Zols, tol=1e-40)
  Pi <- ZZinv%*%(t(Zols)%*%X)
  Xhat <- cbind(1, Zols%*%Pi )
  
  # jackknife first stage
  if (est=="jack"){
    zlist <- pull(Zols, 1)
    xlist <- pull(as.matrix(X, ncol=1), 1)
    hlist <- lapply(zlist, function(x) t(x)%*%ZZinv%*%x)
    Xhat <- t(rbind(1, mapply(levs, h=hlist, z=zlist, x=xlist, MoreArgs=list(Pi))))
  }
  
  # second stage
  XXinv2 <- solve(t(Xhat)%*%Xols, tol=1e-20)
  beta <- solve(t(Xhat)%*%Xols, tol=1e-20)%*%(t(Xhat)%*%Y)
  
  # standard errors
  XXhatinv <- solve(t(Xhat)%*%Xhat, tol=1e-20)
  XXinv <- solve(t(Xols)%*%Xols, tol=1e-20)
  u2 <- (Y - Xols %*% beta)^2
  res2 <-sum(u2)/(n-k)
  covmat <- res2*XXhatinv
  stderr <- sqrt(diag(covmat))
  covered<- (beta[2] + 1.96*stderr[2] >= 1 & beta[2] - 1.96*stderr[2] <= 1)
  
  # heteroskedastic se
  if (se == "HC1"){
    covmat<-XXhatinv%*%t(Xhat)%*%diag(as.vector(u2))%*%Xhat%*%XXhatinv
    stderr <- sqrt(diag(covmat))
  }
  
  # returns
  return(list(b=beta, se=stderr, cov=covered))
}


#### Logit Likelihood Function ###
fnL <- function(theta, yvec,xvec){
  z <- xvec%*%theta
  v1 <- log(exp(z)/(1 + exp(z)))
  v0 <- log(1/(1 + exp(z)))
  inner <- yvec*v1 + (1-yvec)*v0
  return(sum(inner))
}


### nonparametric sieve regression ###
sieve <- function(K,y,x, extra){
  sieve_dta <- as.data.table(cbind(y,1))
  
  # set up sieve polynomial and regress
  for (i in 1:K){
    sieve_dta[,paste0("x",i):=x^i]
  }
  Y<- as.matrix(y,ncol=1)
  Xols<- as.matrix(cbind(sieve_dta[,y:=NULL], extra))
  ols<- solve(t(Xols)%*%Xols, tol = 1e-100)%*%(t(Xols)%*%Y)
  out <- Xols%*%ols
  return(out)
}



### LARF ###
larf <- function(y, d, z, x, xpred, K, fs=""){
    # setup 
    xlog <- cbind(1, x)
    
    # first step: estimating tau0
    if (fs == "parametric"){
      thetat <- rep(.1, ncol(xlog))
      logit <-maxLik(fnL, start=thetat, yvec=z, xvec=xlog)
      tau0 <- exp(xlog%*%logit$estimate)/(1 + exp(xlog%*%logit$estimate))
    }
    if (fs == "nonparametric"){
      sy<- xpred[,1]
      sx<- xpred[,-1]
      tau0 <- sieve(K,z,sy,sx)
    }
    
    # kappa 
    kappa <- 1 - ((d * (1-z))/(1 - tau0)) - (((1-d)*z)/tau0)
    
    # estimator
    xfull <- cbind(1,d,x)
    kxxinv <- solve(t(xfull)%*%diag(as.vector(kappa))%*%xfull, tol=1e-40)
    kxy <- t(xfull)%*%diag(as.vector(kappa))%*%y
    theta <- kxxinv%*%kxy
  
  # returns
  return(list(b=theta, se=stderr))
}


### Bootstrapping ###
# Bootstrapping Standard errors
bs <- function(data, reps){
  set.seed(1234)
  
  # initializing output objects
  namer <- names(data)
  outlet <- array(dim=c(nrow(data), ncol(data), reps))
  dtab <- as.data.table(data)
  dtab[,id:=seq(1:nrow(dtab))]
  
  # drawing row numbers and merging for each rep
  bsints <- matrix(sample(1:nrow(data), nrow(data)*reps, replace=TRUE), ncol=reps)
  for (r in 1:reps){
    bstest <-as.data.table(bsints[,r])
    names(bstest) <- "id"
    rdat <- merge(bstest, dtab, by = "id")
    rmat <- as.matrix(rdat[,id:=NULL])
    outlet[,,r]<- rmat
    dimnames(outlet)[[2]] <- namer
  }
  return(outlet)
}


### Hypothesis Testing ###

# Anderson-Rubin test statistic
artest <- function(betas, bnew, bidx, y,x, z){
  # 2SLS estimates, replacing beta with new beta
  betas[bidx]<-bnew
  
  # regression of uhat on Z to get g(b), Omega
  upred <- (y- cbind(1,x)%*%betas)
  ar.reg0 <- reg(upred, z,se="HC1")
  ar.covmat0 <- ar.reg0$covmat
  ar.gb0 <- ar.reg0$b
  
  # calculating test statistic
  ar.stat0 <- t(ar.gb0)%*%solve(ar.covmat0)%*%ar.gb0
  return(ar.stat0 < qchisq(.95, df=ncol(z)))
}




