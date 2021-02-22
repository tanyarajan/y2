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

reg <- function(Y, X, se="", c=NULL){
  # dimensions
  n <- nrow(as.matrix(X))
  k <- ncol(as.matrix(X))
  
  # solving for beta
  Xols <- cbind(1,X)
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


### wild bootstrap ###
wildboot <- function(gendata,reps, true){
  # Running constrained regression on all except rel4
  selx <- c("E3", "E4", "E5", "T2", "T3", "T4", "T5", 
            "rel1", "rel2", "rel3", "rel4", "rel5", "rel6")
  yreg <- as.matrix(gendata$Y, ncol=1)
  xreg <- as.matrix(subset(gendata, select=selx))
  xconst <- xreg[,which(colnames(xreg)!="rel4")]
  betaout <- reg(yreg,xconst)$b
  beta0 <- c(betaout[1:11], true, betaout[12:13])
  upred <- yreg - cbind(1,xreg)%*%beta0
  cluster <- as.matrix(gendata$person, ncol=1)
  
  
  # Function to run regression on adjusted
  adj.reg <- function(x, b0, uhat, cluster){
    rad.y <- cbind(1, x)%*%b0 + uhat
    wild.reg <- reg(rad.y, x, se="cluster", c=cluster)
    wild.coeff <- rev(rev(wild.reg$b)[1:6])
    wild.se <- rev(rev(wild.reg$se)[1:6])
    wild.wald <- (wild.coeff[4] - true)/(wild.se[4])
    return(wild.wald)
  }
  
  # Rademacher adjustment applied over all bootstrap samples
  set.seed(1234)
  rad.draw <- matrix(sample(c(-1,1), size=length(upred)*reps, prob = rep(1/2, 2), replace=T), ncol=reps)
  wild.out <- apply(rad.draw, 2, function(i){adj.reg(xreg, b0=beta0, upred*i, cluster)})
  
  # Normal clustered regression
  reg.out <- reg(yreg, xreg, se="cluster", c=cluster)
  reg.coeff <- rev(rev(reg.out$b)[1:6])
  reg.se <- rev(rev(reg.out$se)[1:6])
  wald <- (reg.coeff[4] - true)/(reg.se[4])
  
  return(abs(wald) > quantile(abs(wild.out), .975))
}

