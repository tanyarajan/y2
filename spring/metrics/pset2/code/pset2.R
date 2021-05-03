#####################################################
#               Metrics Problem Set 2
#                     Tanya Rajan
# Description: This code runs lasso and ridge 
# regressions over 10,000 Monte Carlo draws 
# generated from the DGP specified in pset 1.
#####################################################

# setup
rm(list=ls())
library(data.table)
library(ggplot2)
library(knitr)
library(kableExtra)
library(MASS)
library(abind)
library(doParallel)
library(gsubfn)
library(pracma)
library(mvnfast)
library(gridExtra)

# set working directory (change to run)
filepath <- "/Users/tanyarajan/Documents/git/y2/spring/metrics/pset2" 
setwd(filepath) 

# set seed
seeder = 123
set.seed(seeder, kind = "L'Ecuyer-CMRG")

# read in code file defining functions
source("code/pset2fns.R")

# setting cores for parallelization
cor<-floor(detectCores(all.tests=FALSE)*.4)
if (is.na(cor) | cor==0){cor <- 1}
registerDoParallel(cores=cor)

#####################################################
#                   Data Generation                 #
#####################################################

drawdata <- function(N, P, M, vcv="ident", rho=NULL){
  # mean 
  mu = rep(0, P+1)
  
  # different variance matrix (either identity or w/ rho covariance)
  if (vcv == "ident"){sigma=eye(P+1)}
  else if (vcv== "auto"){
    sigma1 = eye(P)*(rho - 1)
    sigma2 = matrix(rho, nrow=P, ncol=P)
    sigma = sigma2 - sigma1
    sigma = cbind(sigma, rep(0, P))
    sigma = rbind(sigma, rep(0, P+1))
    sigma[P+1, P+1] = 1
  }
  
  # drawing x, epsilon, and y
  set.seed(seeder)
  out <- lapply(1:M, function(m){rmvn(N, mu, sigma)})
  xep = abind(out, along=3)
  y = xep[,1,] - xep[,2,] + xep[,P+1,]
  x = xep[,-(P+1),]
  return(list(x=x, y=y, ep=xep[,(P+1),]))
}

#####################################################
#           LASSO Objective and Updating            #
#####################################################
# parsing data
parse.data <- function(data){
  y <- data[,1]
  x <- data[,-1]
  return(list(y=y, x=x))
}


# lasso objective
lasso.objective <- function(b, data, lambda){
  list[y,x] <- parse.data(data)
  return(sum((y - x%*%b)^2) + lambda*sum(abs(b)))
}


# lasso update
lasso.update <- function(b, data, lambda){
  # setup
  list[y,x] <- parse.data(data)
  p <- dim(x)[2]
  N <- dim(x)[1]
  
  # calculating bj
  bvec <- b
  for (k in 1:p){
    xnot = x[,-k]
    bnot = bvec[-k]
    betahat = (1/N)*t(x[,k])%*%(y - xnot%*%bnot)
    signhat = (abs(betahat) == betahat)*2 - 1
    meat <- max(abs(betahat) - .5*lambda, 0)
    bvec[k] <- signhat * meat
  }
  return(unlist(bvec))
}

# standardization function
standard <- function(vec){
  return((vec - mean(vec))/std(vec))
}

# lasso 
lasso <- function(b_initial=rep(0,setP), data, lambda, eps=1e-06, max=1000){
  # setup
  list[y,x] <- parse.data(data)
  
  # standardize every regressor in data
  stdX <- apply(x, 2, standard)
  stddata <- cbind(y, stdX)
  
  # loop using stopping criteria
  iter = 1
  diff = 100
  objective = list()
  b_new = b_initial
  while (iter < max & diff > eps){
    b_old = b_new
    b_new = lasso.update(b_old, stddata, lambda)
    objective[iter] <- lasso.objective(b_new, stddata, lambda)
    iter = iter + 1
    diff = max(abs(b_new - b_old))
  }
  if (iter >= max){status<-"Max Iterations Reached"}
  if (diff <= eps){status<-"Tolerance Satisfied"}
  return(list(estimate=b_new, objective=unlist(objective), status=status))
}


#####################################################
#                  Ridge Regression                 #
#####################################################
ridge <- function(data, lambda){
  # setup
  list[y,x] <- parse.data(data)
  
  # ridge estimator
  bhat <- solve(t(x)%*%x + lambda*eye(dim(x)[2]), tol=1e-40)%*%t(x)%*%y
  return(bhat)
  
}


#####################################################
#                 Question 1 Answers                #
#####################################################

# set parameters for Monte Carlo draws
setM <- 10000
setP <- 90
setN <- 100
setL <- 20 

# drawing the sample
sample1 <- drawdata(setN, setP, setM)
data1 <- abind(sample1$y, sample1$x, along=2)

# running lasso across all samples
lasso.loop <- foreach (j=1:setM) %dopar%{
  lasso(data=data1[,,j], lambda = setL)$estimate[1:2] 
}
lasso_results <- matrix(unlist(lasso.loop), ncol=2, byrow=T)
sum(lasso_results[,1])
sum(lasso_results[,2])


# running ols across all samples
ols.loop <- foreach (j=1:setM) %dopar%{
  list[y,x] <- parse.data(data1[,,j])
  reg(y,x, const="none")$b[1:2] 
}
ols_results <- matrix(unlist(ols.loop), ncol=2, byrow=T)
colMeans(ols_results)


# running ridge across all samples
ridge.loop <- foreach (j=1:setM) %dopar%{
  ridge(data=data1[,,j], lambda = setL)[1:2] 
}
ridge_results <- matrix(unlist(ridge.loop), ncol=2, byrow=T)
colMeans(ridge_results)


# finding minimum lambda that returns a vector of zero estimates
bcheck <- 1
lam.loop <- .01
while (bcheck > 0){
  bcheck <- abs(sum(lasso(data=data1[,,100], lambda=lam.loop)$estimate))
  lam.loop <- lam.loop + .01
}
lam.max <- lam.loop - .01

# calculating lasso for all the lambda values between 0.01lambda and max lambda
grid <- (seq(1:1000)/1000)*lam.max
bres <- matrix(0, nrow=1, ncol=5)
for (i in 1:length(grid)){
  bres <- rbind(bres, lasso(data=data1[,,100], lambda=grid[i])$estimate[1:5])
}
bres <- bres[-1,]

# plotting beta_1 through beta_5 as a function of lambda
ldta <- as.data.table(cbind(grid, bres))
names(ldta) <- c("lambda", "b_1", "b_2", "b_3", "b_4", "b_5")
ldta.long <- melt(ldta, "lambda")
ggplot(data=ldta.long, aes(x=lambda, y=value)) +
  geom_line(aes(color=variable)) + theme_bw() + xlab(expression(lambda)) +
  ylab("Estimate")
ggsave("figures/pset2_beta_by_lambda.png")


#####################################################
#             Data Generation for Q2                #
#####################################################
dynamic.dgp <- function(beta, rho, sigma, t){
  x <- matrix(0, nrow = t+1, ncol=50)
  ep <- mvrnorm(t, rep(0, 50), sigma)
  x[1,] <- mvrnorm(1, rep(0, 50), (1/(1-rho^2))*sigma)
  tx <- 2
  while (tx - 1  < t){
    x[tx, ] <- rho*x[tx - 1,] + ep[tx]
    tx<- tx + 1
  }
  y <- x%*%beta + rnorm(t+1, 0, 1)
  return(list(y = y, x = x))
}


#####################################################
#                 Cross Validation                  #
#####################################################
cross.validation <- function(data, lambda_seq, tau, h){
  # setup
  list[t1, t2] <- tau

  # sample splitting
  maxt <- length(data$y)
  cnt <- 1
  yend <- 1
  out <- matrix(0, ncol=length(lambda_seq))
  while (yend < maxt){
    # testing and training data
    train.dta <- cbind(data$y[(cnt + h) : (cnt + t1 - 1 + h)], data$x[cnt:(cnt + t1 - 1),])
    ymax <- cnt + t1 + t2 -1 + h
    if (ymax > maxt){ymax <- maxt}
    test.x <- data$x[(cnt + t1) : (ymax-h),]
    test.y <- data$y[(cnt + t1 + h) : ymax]
    
    # looping through each lambda
    lres <- foreach(l = 1:length(lambda_seq)) %dopar% {
      bhat<-lasso(b_initial=rep(0,50), data = train.dta, lambda = lambda_seq[l])$estimate
      pred.y <- test.x%*%bhat
      return((test.y-pred.y)^2)
    }
    lres <- matrix(unlist(lres), ncol=length(lambda_seq), byrow=F)
    
    # return MSFE
    out <- rbind(out, lres)
    cnt <- cnt + 10
    
    yend <- ymax
  }
  
  #return(colMeans(out[-1,]))
  return(colMeans(out))
}

hi<-cross.validation(data.2b, c(.2,.3), setTau, 1)
print(dim(hi))


#####################################################
#                 Question 2 Answers                #
#####################################################
# setup
setR <- .9
setT <- 200
setTau <- c(100, 10)

# generate data for part b
sig.2b <- (1-setR^2)*eye(50)
b.2b <- rep(0, 50)
b.2b[1] <- 5
data.2b <- dynamic.dgp(b.2b, setR, sig.2b, setT)

# generate data for part c
b.2c <- rnorm(50, 0, 1)
data.2c <- dynamic.dgp(b.2c, setR, sig.2b, setT)

# generate data for part d
sig.tilde <- (matrix(.8, ncol=10, nrow=10) + diag(rep(.2, 10)))*(1-setR^2)
sig.2d <- kronecker(eye(5), sig.tilde)
data.2d <- dynamic.dgp(b.2b, setR, sig.2d, setT)
  
# cross validation
lam20 <- (seq(0:20))/4
res.2b <- cross.validation(data.2b, lam20, setTau, 1)
res.2c <- cross.validation(data.2c, lam20, setTau, 1)
res.2d <- cross.validation(data.2d, lam20, setTau, 1)

# plotting
cvdta <- as.data.table(cbind(lam20, res.2b, res.2c, res.2d))
names(cvdta) <- c("lambda", "part2b", "part2c", "part2d")
p1 <- ggplot(data=cvdta, aes(x=lambda)) + geom_point(aes(y=part2b)) + 
  theme_bw() + xlab(expression(lambda)) + ylab("MSFE") + ggtitle("Part 2b") +
  geom_vline(xintercept = lam20[which.min(res.2b)], color="red")
p2 <- ggplot(data=cvdta, aes(x=lambda)) + geom_point(aes(y=part2c)) + 
  theme_bw() + xlab(expression(lambda)) + ylab("MSFE") + ggtitle("Part 2c") +
  geom_vline(xintercept = lam20[which.min(res.2c)], color="red")
p3 <- ggplot(data=cvdta, aes(x=lambda)) + geom_point(aes(y=part2d)) + 
  theme_bw() + xlab(expression(lambda)) + ylab("MSFE") + ggtitle("Part 2d") + 
  geom_vline(xintercept = lam20[which.min(res.2d)], color="red")

grid.arrange(p1, p2, p3, ncol=2)
ggsave("figures/pset2_msfe.png")


# results at the minimizing MSFE in full sample
lasso(b_initial=rep(0,50), data=data.2b, lambda=lam20[which.min(res.2b)])
lasso(b_initial=rep(0,50), data=data.2c, lambda=lam20[which.min(res.2c)])
lasso(b_initial=rep(0,50), data=data.2d, lambda=lam20[which.min(res.2d)])


