#####################################################
#         Labor Pset 6: Hryshko Replication
#####################################################

# setup
rm(list=ls())
require(data.table)
require(ggplot2)
require(gsubfn)
require(knitr)
require(kableExtra)

# set working directory (change to run)
filepath <- "/Users/tanyarajan/Documents/git/psets_heckmanlabor/pset6" 
setwd(filepath) 

# load data
df <- as.data.table(read.csv(file = "data/nlsy_clean.csv"))

# Function to extract data
extract_data <- function(datain){
  data <- copy(datain)
  y <- data$y
  y_lag1 = data$ly1
  y_lag2 = data$ly2
  y_lead1 = data$fy1
  y_lead2 = data$fy2
  yforward = as.matrix(data[, c("caseid", "T", "y", "ly1", "ly2", "fy1", "fy2"):= NULL])

  # need to figure out yforward for first moment
  return(list(y=y, y_lag1=y_lag1, y_lag2=y_lag2, 
              y_lead1=y_lead1, y_lead2=y_lead2, yforward=yforward))
}

df2 <- copy(df)
test = as.matrix(df2[, c("caseid", "T", "y", "ly1", "ly2", "fy1", "fy2"):= NULL])

list[a,b,c,d,e,f] <- extract_data(df)

#####################################################
#                 Moment Conditions                 #
#####################################################
# split vector omega into individual parameters
parse_omega <- function(omega){
  s_b   = omega[1] 
  s_xi  = omega[2]
  theta = omega[3]
  s_ep  = omega[4]
  s_u   = omega[5]
  return(list(s_b=s_b, s_xi=s_xi, theta=theta, s_ep=s_ep, s_u=s_u))
}

# Moment 1: from covariance gamma_k
moment1 <- function(omega, y, yforward, Tmax){
  list[s_b, s_xi, theta, s_ep, s_u] <- parse_omega(omega)
  
  gmat <- matrix(NA, ncol=(Tmax-3), nrow=dim(yforward)[1])
  for (i in 3:Tmax){
    idx = i - 3
    gmat[,idx] <- y*yforward[,idx]
  }

  # currently just using gamma3, but technically shoudl be for gamma_k, k>=3
  return(gmat - s_b)
}

# Moment 2: from summation moment (below eqn 9 in the paper)
moment2 <- function(omega, y, y_lag1, y_lag2, y_lead1, y_lead2){
  list[s_b, s_xi, theta, s_ep, s_u] <- parse_omega(omega)
  # need to fill in
  gsum <- y_lag2 + y_lag1 + y + y_lead1 + y_lead2
  gmom <-y*gsum
  return(gmom - 5*s_b - s_xi)
}

# Moment 3: from covariance gamma_2
moment3 <- function(omega, y, y_lag2){
  list[s_b, s_xi, theta, s_ep, s_u] <- parse_omega(omega)
  g2 <- y*y_lag2
  return(g2 - s_b + theta*s_ep)
}

# Moment 4: from covariance gamma_1
moment4 <- function(omega, y, y_lag1){
  list[s_b, s_xi, theta, s_ep, s_u] <- parse_omega(omega)
  g1 <- y*y_lag1
  return(g1 - s_b + (theta - 1)^2 * s_ep + s_u)
}

# Moment 5: from covariance gamma_0
moment5 <- function(omega, y){
  list[s_b, s_xi, theta, s_ep, s_u] <- parse_omega(omega)
  g0 <- y*y
  return(g0 - s_xi - s_b - (1 + (1-theta)^2 + theta^2)*s_ep - 2*s_u)
}


#####################################################
#                 Gradient Functions                #
#####################################################
# Gradient 1
grad1 <- function(omega, Tmax){
  m1num <- Tmax - 3
  list[s_b, s_xi, theta, s_ep, s_u] <- parse_omega(omega)
  return(t(matrix(c(1, 0, 0, 0, 0), nrow=5, ncol=m1num))) 
}

# Gradient 2
grad2 <- function(omega){
  list[s_b, s_xi, theta, s_ep, s_u] <- parse_omega(omega)
  return(c(-5, -1, 0, 0, 0))
}

# Gradient 3
grad3 <- function(omega){
  list[s_b, s_xi, theta, s_ep, s_u] <- parse_omega(omega)
  return(c(-1, 0, s_ep, 0, 0))
}

# Gradient 4
grad4 <- function(omega){
  list[s_b, s_xi, theta, s_ep, s_u] <- parse_omega(omega)
  return(c(-1, 0, 2*(theta-1)*s_ep, (theta-1)^2, 1))
}

# Gradient 5
grad5 <- function(omega){
  list[s_b, s_xi, theta, s_ep, s_u] <- parse_omega(omega)
  return(c(-1, -1, 2*(1-theta)*s_ep + 2*theta, 1 + (1-theta)^2 + theta^2, -2))
}


#####################################################
#              GMM Estimation (EWMD)                #
#####################################################
# stacking moments
gmm_stack <- function(omega, data, Tmax){
  # extract data
  list[y, y_lag1, y_lag2, y_lead1, y_lead2, y_lead3] <- extract_data(data)

  # form moments
  m1 <- moment1(omega, y, y_lead3, Tmax)
  m2 <- moment2(omega, y, y_lag1, y_lag2, y_lead1, y_lead2)
  m3 <- moment3(omega, y, y_lag2)
  m4 <- moment4(omega, y, y_lag1)
  m5 <- moment5(omega, y)
  mstack = cbind(m1, m2, m3, m4, m5)
  mstack[is.na(mstack)] <- 0
  return(mstack)
}

gmm_stack(c(1,2,3,4,5), df, 15)



#-----------------------------------
#      Parameter Estimates 
#-----------------------------------
# gmm objective
gmm_objective <- function(omega, data, model, Tmax){
  
  # setup
  if (model == 1){omega[2] <- 0}
  moments <- colMeans(gmm_stack(omega, data, Tmax), na.rm=T)
  
  # GMM estimation
  return(t(moments)%*%moments)
}

# optimization
omega_guess = c(1, 0, 1, 1, 1) # obviously change to make this better
results1 <-optim(omega_guess, gmm_objective, data=df, Tmax=15, model=1, method="L-BFGS-B", lower = c(0,0,-1,0,0), hessian = TRUE)
results2 <-optim(omega_guess, gmm_objective, data=df, Tmax=15, model=2, method="L-BFGS-B", lower = c(0,0,-1,0,0), hessian = TRUE)


# displaying coefficients
coef1 <- results1$par
coef2 <- results2$par

#-----------------------------------
#        Standard Errors
#-----------------------------------

# Sigma function at center of asymptotic variance formula
compute_meat <- function(omega, data, model, Tmax){
  # W matrix
  mstack <- gmm_stack(omega, data, Tmax)
  sqsum <- matrix(0, dim(mstack)[2], dim(mstack)[2])
  for (i in 1:dim(data)[1]){
    sqsum <- sqsum + mstack[i,]%*%t(mstack[i,])
  }
  meat <- (1/dim(data)[1])*sqsum
  return(meat)
}

# computing standard errors
gmm_se <- function(omega, data, model, Tmax){
  # meat
  meat<-compute_meat(omega, data, model, Tmax)
  
  # gradient
  G <- rbind(grad1(omega, Tmax), grad2(omega), grad3(omega), grad4(omega), grad5(omega))
  if (model == 1){
    G <- rbind(grad1(omega, Tmax), grad2(omega), grad3(omega), grad4(omega), grad5(omega))
    G <- G[,-2]
    }

  # variance
  bread <- G%*%solve(t(G)%*%G, tol=1e-20)
  variance <- t(bread)%*%meat%*%bread
  return(sqrt(diag(variance)/dim(data)[1]))

}

se1 <- gmm_se(results1$par, df, 1, 15)
se2 <- gmm_se(results2$par, df, 2, 15)

             
fmtse1 <- c(0,0,0,0)
for (j in 1:length(se1)){
  fmtse1[j] <- as.character(paste0("(",round(se1[j], 4),")"))
}
fmtse1 <- append(fmtse1, "-", after=1)

fmtse2 <- c(0,0,0,0,0)
for (j in 1:length(se2)){
  fmtse2[j] <- as.character(paste0("(",round(se2[j], 4),")"))
}

outtable <- rbind(round(coef1, 4), fmtse1, round(coef2,4), fmtse2)
rownames(outtable) <- c("HIP Model", "  ", "HIP + RW Model", " ")
#outtable <- cbind(names, outtable)
colnames(outtable)<-c( "$\\sigma_\\beta$","$\\sigma_\\xi$","$\\theta$","$\\sigma_\\varepsilon$", "$\\sigma_u$")

kable(outtable, format="latex", booktabs=TRUE, escape=FALSE) %>%
  write(file="tables/tr_table2.tex")



