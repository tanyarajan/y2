#####################################################
#               Metrics Problem Set 3
#                     Tanya Rajan
# Description: This code defines functions for both
# of the computational questions on this problem set.
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

# Probit (but this is wrong, I think)
fnP <- function(theta, yvec, xvec){
  z <- xvec%*%theta
  v1 <- log(pnorm(z))
  v0 <- log(1 - pnorm(z))
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

