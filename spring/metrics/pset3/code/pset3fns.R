#####################################################
#               Metrics Problem Set 1
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

# Euclidean norm
enorm <- function(x){
  sqrt(sum((x)^2))
}

########### Regression Functions ##############

### leverage calculation for jackknife ###
levs <- function(h, pi, z, x){
  h<-as.numeric(h)
  return((z%*%pi - h*x)/(1-h))
}

### OLS regression function ###

reg <- function(Y, X, se="", c=NULL, const="yes"){
  # dimensions
  n <- nrow(as.matrix(X))
  k <- ncol(as.matrix(X))
  
  # solving for beta
  Xols <- cbind(1,X)
  if (const=="none"){Xols <- X}
  XXinv <-  solve(t(Xols)%*%Xols, tol=1e-100)
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
  
  # cluster-robust
  if (se == "cluster"){
    clist <- as.numeric(unique(c))
    omega <- matrix(0, dim(Xols)[2], dim(Xols)[2])
    XX.c <- matrix(0, dim(Xols)[2], dim(Xols)[2])
    for (cl in clist){
      xg <- Xols[c==cl,]
      u <- (Y- Xols%*%beta)[c==cl,]
      omega <- omega + t(xg)%*%u%*%t(u)%*%xg
      XX.c <- XX.c + t(xg)%*%xg
    }
    cN <- length(clist)
    div <- 1/cN
    adjust <- (cN/(cN - 1))*((n-1)/(n - k))
    covmat <- adjust*solve(XX.c, tol=1e-20)%*%(div*omega)%*%solve(div*XX.c, tol=1e-20)
    stderr <- sqrt(diag(covmat))
  }
  
  # returns
  return(list(b=beta, se=stderr, cov=covered, covmat=covmat))
}

