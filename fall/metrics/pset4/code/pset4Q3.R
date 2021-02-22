#####################################################
#               Metrics Problem Set 4
#                      Question 3
#                     Tanya Rajan
# Description: This file reproduces figures 10, 11 of
# Kellogg, Mogstad, Pouliot, and Torgovitsky (2020)
# using data from Abadie and Gardeazabal (2003).
#####################################################

# setup
rm(list=ls())
require(data.table)
require(ggplot2)
require(pracma)
require(nloptr)
require(optimx)
require(gurobi)

# set working directory (change to run)
filepath <- "/Users/tanyarajan/Documents/git/psets/metrics_fall2020/pset4" 
setwd(filepath) 

# read in code file defining functions
source("code/pset4fns.R")

#####################################################
#                Setting up Data                    #
#####################################################
# read in data
datain <- as.data.table(read.csv("data/basque.csv"))
datain <- datain[regionno!=1,]

# gdp data only
gdponly <- subset(datain,select=c("gdpcap","year", "regionno"))
Xsc <- gdponly[year <= 1969 & regionno != 17,]
X1 <- gdponly[year <= 1969 & regionno == 17,]

# matching and outcome variables
Xcontrol <- as.matrix(dcast(Xsc, regionno~year, value.var="gdpcap")[,regionno:=NULL])
Xtreat <- t(as.matrix(dcast(X1, regionno~year, value.var="gdpcap")[,regionno:=NULL]))
Ycontrol <- as.matrix(dcast(gdponly, regionno~year, value.var="gdpcap")[regionno!=17,][,regionno:=NULL])
Ytreat <- as.matrix(dcast(gdponly, regionno~year, value.var="gdpcap")[regionno==17,][,regionno:=NULL])

# number of control regions
Nsc <- nrow(Xcontrol)

#####################################################
#               Weighting Functions                 #
#####################################################

# Synthetic control weights function
sc.weighter <- function(x, x1){
  # objective function
  minimizer <- function(w,x,x1){
    xtreat <- matrix(t(x1), ncol=1)
    xcontr <- matrix(t(x), nrow = ncol(x))
    return((t(xtreat - xcontr%*%w)%*%(xtreat - xcontr%*%w)))
  }

 # optimization with constraints
  wguess <- rep(1/Nsc, Nsc)
  tol        <- .1 # higher tolerance helps here
  const      <- rbind(rep(1, length(wguess)), rep(-1, length(wguess)), diag(Nsc))
  equal      <- c(1-tol, -1-tol, rep(0, Nsc)) #equatlity and bound constraints
  weights <- constrOptim(wguess, minimizer, x=x, x1=x1, 
                            ui=const, ci=equal, gr=NULL, 
                            control=list(reltol=1e-100, abstol=1e-100))$par
  return(weights)
}


# Penalized SC
scpen.weighter <- function(x, x1, pi){
  # objective function
  minimizer.pen <- function(w, pi, x,x1){
    xtreat1 <- matrix(t(x1), ncol=1)
    xcontr1 <- matrix(t(x), nrow = ncol(x))
    V <- diag(ncol(x))
    term1 <- (t(xtreat1- xcontr1%*%w)%*%V%*%(xtreat1 - xcontr1%*%w))
    
    x1 <- matrix(x1, ncol=1)
    xtreat2 <- repmat(x1, 1, nrow(x))
    xcontr2 <- matrix(t(x), nrow = ncol(x))
    term2 <- w%*%colSums((xtreat2 - xcontr2)^2)
    return((1-pi)*term1 + pi*sum(unlist(term2)))
  }
    
  # optimization with constraints
  wguess <- rep(1/Nsc, Nsc)
  tol        <- .06 # higher tolerance helps here!
  const      <- rbind(rep(1, length(wguess)), rep(-1, length(wguess)), diag(Nsc))
  equal      <- c(1-tol, -1-tol, rep(0, Nsc))
  weights <- constrOptim(wguess, minimizer.pen, x=x, x1=x1, pi=pi, 
                         ui=const, ci=equal, gr=NULL, 
                         control=list(reltol=1e-100, abstol=1e-100))$par
  return(weights)
}



# Matching weights function
match.weighter <- function(x, x1, k){
  x1 <- matrix(x1, ncol=1)
  xtreat <- repmat(x1, 1, nrow(x))
  xcontr <- matrix(t(x), nrow = ncol(x))
  dist <- colSums((xtreat - xcontr)^2)
  weights <- (dist <= sort(dist)[k])/k
  return(weights)
}


#####################################################
#         Functions for Estimators                  #
#####################################################

# synthetic control (penalized and non-penalized)
sc <- function(y,x,x1,pen=NULL, pi=NULL){
  if (is.null(pen)){w<-sc.weighter(x, x1)}
  if (!is.null(pen)){w<-scpen.weighter(x,x1,pi)}
  return(as.vector(w%*%y))
}

# matching
nnmatch <- function(y,x,x1,k){
  w<-match.weighter(x, x1, k)
  return(as.vector(w%*%y))
}

# MASC
masc <- function(y,x, x1, tuner){
  k <- tuner[1] 
  phi <- tuner[2]
  ma <- nnmatch(y, x, x1, k)
  sc <- sc(y, x, x1)
  return(phi*ma + (1-phi)*sc)
}

#####################################################
#          Cross-Validation Functions               #
#####################################################

# Cross validation
xv.masc <- function(x, x1, k, folds=8:14){
  numerator <- list()
  denom <- list()
  error <- list()
  ests <- list()
  
  # Get Y values for the treatment
  idx <- 1 
  for (t in folds){
    # treatment year one period in the future
    y0 <- x[,(t+1)]
    y1 <- x1[(t+1),]
    
    # computing phi
    gamma.ma <- nnmatch(y0, x[,1:t], x1[1:t,], k)
    gamma.sc <- sc(y0, x[,1:t], x1[1:t,])
    numerator[[idx]]<- (gamma.ma - gamma.sc)*(y1 - gamma.sc)
    denom[[idx]]<- (gamma.ma - gamma.sc)^2
    
    # saving estimators
    ests[[idx]] <- c(y1, gamma.ma, gamma.sc)
    idx <- idx + 1
  }
  
  # calculate phi and bound it btwn 0 and 1
  phi <- sum(unlist(numerator))/sum(unlist(denom))
  if (phi > 1){phi <- 1}
  if (phi < 0){phi <- 0}
  
  outerr <- lapply(1:(idx-1), function(i){ests[[i]][1] - 
      (phi*ests[[i]][2] + (1-phi)*ests[[i]][3])})
  error<- sum(unlist(outerr)^2)
  
  # get error
  return(list(phi = phi, error = sum(unlist(error))))
}

# cross validate matching
xv.nn <- function(x, x1, k, folds=8:14){
  error <- list()
  ests <- list()
  
  # Get Y values for the treatment
  idx <- 1 
  for (t in folds){
    # treatment year one period in the future
    y0 <- x[,(t+1)]
    y1 <- x1[(t+1),]
    gamma.ma <- nnmatch(y0, x[,1:t], x1[1:t,], k)

    # saving estimators
    ests[[idx]] <- y1-gamma.ma
    idx <- idx + 1
  }
  
  # calculate errors
  outerr <- lapply(1:(idx-1), function(i){ests[[i]]})
  error<- sum(unlist(outerr)^2)
  
  # get error
  return(error)
}

# cross validate scpen
xv.scpen <- function(x, x1, pi, folds=8:14){
  error <- list()
  ests <- list()
  
  # Get Y values for the treatment
  idx <- 1 
  for (t in folds){
    # treatment year one period in the future
    y0 <- x[,(t+1)]
    y1 <- x1[(t+1),]
    gamma.scpen <- sc(y0, x[,1:t], x1[1:t,], pen="T", pi=pi)
    
    # saving estimators
    ests[[idx]] <- y1-gamma.scpen
    idx <- idx + 1
  }
  
  # calculate errors
  outerr <- lapply(1:(idx-1), function(i){ests[[i]]})
  error<- sum(unlist(outerr)^2)
  
  # get error
  return(error)
}


#####################################################
#                  Run Everything                   #
#####################################################

#### Cross Validate to Get Optimal Tuners ####
# x-val scpenn
grid.pi <- 1:100/100
err.cvscpen <- lapply(grid.pi, function(g){xv.scpen(x=Xcontrol, x1=Xtreat, g)})
opt.pi <- grid.pi[which(err.cvscpen==min(unlist(err.cvscpen)))]

# x-val masc
grid.m <- 1:10
phi.cvmasc <- lapply(grid.m, function(g){xv.masc(Xcontrol, Xtreat, g)$phi})
err.cvmasc <- lapply(grid.m, function(g){xv.masc(Xcontrol, Xtreat, g)$error})
opt.k.masc <- which(err.cvmasc==min(unlist(err.cvmasc)))
opt.phi <- phi.cvmasc[[opt.k.masc]]

# x-val nn
err.cvnn <- lapply(grid.m, function(g){xv.nn(Xcontrol, Xtreat, g)})
opt.k.nn <- which(err.cvnn==min(unlist(err.cvnn)))[[1]]



#### Calculate Weights and Estimates ####
trtyr<-which(colnames(Ycontrol)=="1969")

out.masc <- masc(Ycontrol, Xcontrol[,1:trtyr], Xtreat[1:trtyr], c(opt.k.masc, opt.phi))
out.sc <- as.vector(sc(Ycontrol, Xcontrol[,1:trtyr], Xtreat[1:trtyr]))
out.scpen <- as.vector(sc(Ycontrol, Xcontrol[,1:trtyr], Xtreat[1:trtyr], pen="T", pi=opt.pi))
out.nn <- as.vector(nnmatch(Ycontrol, Xcontrol[,1:trtyr], Xtreat[1:trtyr], opt.k.nn))
  
#### Plot outcomes ####
plot.align <- function(out){
  data <- as.data.table(cbind(out, t(Ytreat), as.numeric(colnames(Ycontrol))))
  names(data) <- c("SC", "Basque", "Year")
  dataL <- melt(data, id.vars="Year")
  ggplot(data=dataL, aes(x=Year, y=value, color=variable, group=variable)) + 
    geom_point() + geom_line() + theme_bw()
}
plot.align(out.sc)
plot.align(out.nn)
plot.align(out.masc)
plot.align(out.scpen)

#### Plot differences ####
diff.masc <- as.numeric(Ytreat) - out.masc
plotdiff <- as.data.table(cbind(diff.masc*1000, as.numeric(colnames(Ycontrol))))
names(plotdiff) <- c("Difference", "Year")
meaneff <- plotdiff[Year>=1969, mean(Difference)]
plotdiffL <- melt(plotdiff, id.vars="Year")
difgraph<-ggplot(data=plotdiffL, aes(x=Year, y=value)) + 
  geom_line() + theme_bw() + geom_hline(yintercept=0, linetype="dotted") +
  geom_hline(yintercept=meaneff, color="coral") + 
  geom_vline(xintercept=1969.5, linetype="dotted", color="cornflowerblue") +
  ylab("GDP per Capita") + annotate("text", x=c(1960,1969), 
                                    y=c(meaneff+50, -1550), 
                                    label=c("Mean Effect", "Treatment"))
ggsave(difgraph, width = 7, height = 5, file = "figures/tr_masc.png")


### Plot difference from MASC ####
#out.scpen <- as.vector(sc(Ycontrol, Xcontrol[,1:trtyr], Xtreat[1:trtyr], pen="T", pi=.99))
start<-which(colnames(Ycontrol)=="1969")
end <- ncol(Ycontrol)
diff.ma <- out.masc[start:end]-out.nn[start:end]
diff.sc <- out.masc[start:end]-out.sc[start:end]
diff.scpen <- out.masc[start:end]-out.scpen[start:end]
years <- as.numeric(colnames(Ycontrol))[start:end]
plots <- as.data.table(cbind(diff.ma*1000, diff.sc*1000, diff.scpen*1000, years))
names(plots) <- c("MA", "SC", "Penalized", "Year")
performance<-ggplot(data=plots, aes(x=Year)) + 
  geom_line(aes(y=MA), color="darkgoldenrod") + 
  geom_line(aes(y=SC), color="cornflowerblue") +
  geom_line(aes(y=Penalized), color="coral") + theme_bw() +
  ylab("Difference in GDP (Estimator - MASC)")
ggsave(performance, width = 7, height = 5, file = "figures/tr_perf.png")


