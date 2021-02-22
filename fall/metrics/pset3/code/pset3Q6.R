#####################################################
#               Metrics Problem Set 3
#                      Question 6
#                     Tanya Rajan
# Description: This code simulates data using the 
# Mogstad and Torgovitsky (2018) Section 2.4 setup.
# It then sets up and solves the optimization problem
# to find upper and lower bounds for the ATT target
# parameter.
#####################################################

# setup
rm(list=ls())
require(data.table)
require(ggplot2)
require(slam)
require(gurobi)

# set working directory (change to run)
filepath <- "/Users/tanyarajan/Documents/git/psets/metrics_fall2020/pset3" 
setwd(filepath) 

# read in code file defining functions
source("code/pset3fns.R")

# setting cores for parallelization
cor<-floor(detectCores(all.tests=FALSE)*.7)
if (is.na(cor) | cor==0){cor <- 1}
registerDoParallel(cores=cor)


#####################################################
#           Generating Simulated Data               #
#####################################################
set.seed(1234)
N  <- 50000
ud <- seq(1:N)/N

# MTRs from Section 2.4
m0  <- 0.9 - 1.1*ud + 0.3*(ud^2)
m1  <- 0.35 - 0.3*ud - 0.05*(ud^2)
mte <- -0.55 + 0.8*ud - 0.35*(ud^2)

# Generating Zs
Z <- sample(1:4, N, replace=T)
data <- as.data.table(cbind(ud,Z))

# Generating Ds
pzlist <- c(.12, .29, .48, .78)
for (z in 1:4){
  prob = pzlist[z]
  nz = nrow(data[Z==z])
  data[Z==z, D:=(ud < prob)]
  data[Z==z, pz:= prob]
}

# Generating Y
data[,Y :=m1*D + m0*(1-D)]

#####################################################
#       Calculating Optimization Constraints        #
#####################################################

#### Values for beta_s ####
# IV s-function 
data[,s_iv := (Z - mean(Z))/cov(D,Z)]

# IV beta_s
bs.iv <- (mean(data$Y*data$Z) - mean(data$Y)*mean(data$Z))/
  (mean(data$D*data$Z) - mean(data$D)*mean(data$Z))

# 2SLS s-function
data[,z1 := (Z==1)][,z2 := (Z==2)][,z3 := (Z==3)][,z4 := (Z==4)]
xtilde  <- t(cbind(1, data$D))
ztilde  <- t(cbind(data$z1, data$z2, data$z3, data$z4))
pi      <- xtilde%*%t(ztilde)%*%solve(ztilde%*%t(ztilde))
tsls_c1 <- solve(pi%*%((1/N)*ztilde%*%t(xtilde)), tol=1e-20)%*%pi%*%ztilde
data[,s_tsls := tsls_c1[2,]]

# 2SLS beta_s
tsls_c2 <- solve(pi%*%((1/N)*ztilde%*%t(xtilde)), tol=1e-20)%*%(pi%*%((1/N)*ztilde%*%data$Y))
bs.tsls <- tsls_c2[2,]


#### Parametric and non-parametric functions to compute gamma_s ####
# Bernstein Polynomial
bernstein <- function(u,k,K,pz){
  out<-sapply(k:K, function(i){(-1)^(i-k)*choose(K, i)*choose(i,k)*u^i})
  if(length(u) == 1){return(out)}
  else{return(rowSums(out))}
}

# Constant spline (cutpoints are the propensity scores)
spline <- function(u, k, K, pz){
  qs <- c(0, pz, 1)
  idx <- k+1 # to account for 0 index not existing
  spind <- (u >= qs[idx] & u < qs[idx+1] )
  return(as.numeric(spind))
}


#### Computing gamma_s ####
integrator <- function(sfn, k, K, pz, func){
  
  # setting bounds of integration
  low1 = rep(0,4); up1 = pz
  low0 = pz; up0 = rep(1,4)
  
  # integrating the bases
  fn <- callit(func)
  bk1<-sapply(1:4, function(i){integrate(fn, lower = low1[i], upper = up1[i], 
                                         k=k, K=K, pz=pz)$value})
  bk0<-sapply(1:4, function(i){integrate(fn, lower = low0[i], upper = up0[i], 
                                         k=k, K=K, pz=pz)$value})
  
  # calculating final gamma
  s <- callit(sfn)
  gs0 = mean(s*bk0); gs1 = mean(s*bk1);
  if (sfn == "s_att"){gs0 <- -gs1}
  if (sfn == "s_atu"){gs1 <- -gs0}
  return(list(gs0 = gs0, gs1 = gs1))
}


# Function to calculate gamma_s terms
calc.gs <- function(K){
  out <- list()
  out$piv    <- sapply(0:K, function(i){integrator("s_iv", i, K, pzlist, "bernstein")})
  out$ptsls  <- sapply(0:K, function(i){integrator("s_tsls", i, K, pzlist, "bernstein")})
  out$patt   <- sapply(0:K, function(i){integrator("s_att", i, K, pzlist, "bernstein")})
  out$npiv   <- sapply(0:5, function(i){integrator("s_iv", i, K, pzlist, "spline")})
  out$nptsls <- sapply(0:5, function(i){integrator("s_tsls", i, K, pzlist, "spline")})
  out$npatt  <- sapply(0:5, function(i){integrator("s_att", i, K, pzlist, "spline")})
  return(out)
}

#####################################################
#           Solving the Linear Program              #
#####################################################
# Subsetting to unique Z values and weights
wdata  <- data
wdata  <- unique(wdata, by=c("Z"))
wdata  <- wdata[order(Z)]

# Calculating weights for all methods
s_iv   <- wdata$s_iv
s_tsls <- wdata$s_tsls
ED     <- nrow(data[D==1])/nrow(data)
s_att  <- 1/ED
s_atu  <- 1/(1-ED)

# Calculating all gammas
gammas <- lapply(0:20, calc.gs)

# Function to solve the optimization problem
optimizer <- function(K, gammas, sense, type, addl="none"){
  idx = K+1 # accounting for gamma_0 index = 1
  
  # parametric constraints and objective
  if (type == "par"){
    ivconst    <- c(unlist(gammas[[idx]]$piv[1,]), unlist(gammas[[idx]]$piv[2,]))
    tslsconst  <- c(unlist(gammas[[idx]]$ptsls[1,]), unlist(gammas[[idx]]$ptsls[2,]))
    objective  <- c(unlist(gammas[[idx]]$patt[1,]), unlist(gammas[[idx]]$patt[2,]))
  }

  # nonparametric constraints and objective
  if (type == "nonpar"){
    K= 5
    idx = K+1
    ivconst    <- c(unlist(gammas[[idx]]$npiv[1,]), unlist(gammas[[idx]]$npiv[2,]))
    tslsconst  <- c(unlist(gammas[[idx]]$nptsls[1,]), unlist(gammas[[idx]]$nptsls[2,]))
    objective  <- c(unlist(gammas[[idx]]$npatt[1,]), unlist(gammas[[idx]]$npatt[2,]))
  }
  
  # setting constraints and operators
  constraints  <- rbind(ivconst, tslsconst)
  rhsvals      <- c(bs.iv, bs.tsls)
  operators    <- c('=', '=')
  
  # adding shape constraints if specified
  if (addl == "decreasing"){
    add0 <- matrix(0, nrow=K, ncol=K+1)
    addmat <- add0
    for (i in 1:K){
      addmat[i,i] <- 1
      addmat[i,(i)+1] <- -1
    }
    shape_const0 <- cbind(addmat, add0)
    shape_const1 <- cbind(add0, addmat)
    shaper <- rbind(shape_const0, shape_const1)
    
    # updating matrices
    constraints <- rbind(constraints, shaper)
    rhsvals     <- c(rhsvals, rep(0, (2*K)))
    operators   <- c(operators, rep(">", (2*K)))
  }
  
  # building the model
  model <- list()
  model$A          <- constraints
  model$obj        <- objective
  model$modelsense <- sense
  model$rhs        <- rhsvals
  model$sense      <- operators
  model$lb         <- 0
  model$ub         <- 1
  return(gurobi(model)$objval)
}

# Running through all specifications
p_max     <- sapply(1:19, optimizer, gammas=gammas, sense='max', type='par')
p_min     <- sapply(1:19, optimizer, gammas=gammas, sense='min', type='par')
pdec_max  <- sapply(1:19, optimizer, gammas=gammas, sense='max', type='par', addl='decreasing')
pdec_min  <- sapply(1:19, optimizer, gammas=gammas, sense='min', type='par', addl='decreasing')
np_max    <- sapply(1:19, optimizer, gammas=gammas, sense='max', type='nonpar')
np_min    <- sapply(1:19, optimizer, gammas=gammas, sense='min', type='nonpar')
npdec_max <- sapply(1:19, optimizer, gammas=gammas, sense='max', type='nonpar', addl='decreasing')
npdec_min <- sapply(1:19, optimizer, gammas=gammas, sense='min', type='nonpar', addl='decreasing')

# Calculating ATT
att <- mean(m1[data$D==1] - m0[data$D==1])
attvec <- rep(att, length(p_max))

# Creating data table to graph
g  <-as.data.table(cbind(seq(1:19), attvec, p_max, p_min, pdec_max, pdec_min, 
                         np_max, np_min, npdec_max, npdec_min))
gL <- melt(g, id.vars="V1")
gL[,shape:=as.numeric(grepl("dec",variable))][variable=="attvec",shape:=2]
gL[,shape:=as.character(shape)]
gL[,nonpar:=grepl("np",variable)]

#### Graphing ####
# Color and legend options
colors    <- c( "black",rep(c("dodgerblue", "dodgerblue", "coral", "coral"),2))
linetypes <- c("dotted", "solid", "solid","solid", "solid", 
               "dashed", "dashed", "dashed", "dashed")
ll        <- list()
ll$name   <- ""
ll$breaks <- c("attvec", "p_max","pdec_max", "np_max", "npdec_max")
ll$labels <- c("ATT", "P","P, decr.", "NP", "NP, decr.")

# Graph
graph <- ggplot(data = gL, aes(x=V1, y=value, color=variable, linetype=variable)) + 
  scale_color_manual(values=colors, name=ll$name, breaks=ll$breaks, labels=ll$labels) + 
  scale_linetype_manual(values=linetypes, name=ll$name, breaks=ll$breaks, labels=ll$labels) +
  geom_point() + geom_line() + theme_bw() + ylab("Upper and Lower Bounds") + 
  xlab("Polynomial Degree") + theme(legend.position="bottom", 
                                    legend.key.size=unit(2, "lines"), 
                                    text=element_text(size=16),
                                    legend.text=element_text(size=13)) 
ggsave(graph, width = 7, height = 5, file = "figures/tr_bounds.png")


