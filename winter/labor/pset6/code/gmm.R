#####################################################
#         Labor Pset 6: Hryshko Replication
#####################################################

# setup
rm(list=ls())
require(data.table)
require(ggplot2)
require(gsubfn)

# set working directory (change to run)
filepath <- "/Users/tanyarajan/Documents/git/psets_heckmanlabor/pset6" 
setwd(filepath) 

# load data
df <- as.data.table(read.csv(file = "data/nlsy_clean.csv"))

# Function to extract data
extract_data <- function(data){
  y <- data$y
  y_lag1 = data$ly1
  y_lag2 = data$ly2
  y_lead1 = data$fy1
  y_lead2 = data$fy2
  # need to figure out yforward for first moment
  yforward = 0
  return(list(y=y, y_lag1=y_lag1, y_lag2=y_lag2, 
              y_lead1=y_lead1, y_lead2=y_lead2, yforward=yforward))
}


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
moment1 <- function(omega, y, yforward){
  list[s_b, s_xi, theta, s_ep, s_u] <- parse_omega(omega)
  # need to fill in
}

# Moment 2: from summation moment (below eqn 9 in the paper)
moment2 <- function(omega, y, y_lag, y_lag2, y_lead1, y_lead2){
  list[s_b, s_xi, theta, s_ep, s_u] <- parse_omega(omega)
  # need to fill in
  gsum <- y_lag2 + y_lag1 + y + y_lead1 + y_lead2
  gmom <- mean(y*gsum, na.rm=TRUE)
  return(gmom - 5*s_b - s_xi)
}

# Moment 3: from covariance gamma_2
moment3 <- function(omega, y, y_lag2){
  list[s_b, s_xi, theta, s_ep, s_u] <- parse_omega(omega)
  g2 <- mean(y*y_lag2, na.rm=TRUE)
  return(g2 - s_b + theta*s_ep)
}

# Moment 4: from covariance gamma_1
moment4 <- function(omega, y, y_lag1){
  list[s_b, s_xi, theta, s_ep, s_u] <- parse_omega(omega)
  g1 <- mean(y*y_lag1, na.rm=TRUE)
  return(g1 - s_b + (theta - 1)^2 * s_ep + s_u)
}


# Moment 5: from covariance gamma_0
moment5 <- function(omega, y){
  list[s_b, s_xi, theta, s_ep, s_u] <- parse_omega(omega)
  g0 <- mean(y*y, na.rm=TRUE)
  return(g0 - s_xi - s_b - (1 + (1-theta)^2 + theta^2)*s_ep - 2*s_u)
}


#####################################################
#                 GMM Objective Fn                  #
#####################################################

gmm_objective <- function(omega, data){
  # extract data
  list[y, y_lag1, y_lag2, y_lead1, y_lead2, yforward] <- extract_data(data)
  
  # form moments
  m1 <- moment1(omega, y, yforward)
  m2 <- moment2(omega, y, y_lag1, y_lag2, y_lead1, y_lead2)
  m3 <- moment3(omega, y, y_lag2)
  m4 <- moment4(omega, y, y_lag1)
  m5 <- moment5(omega, y)
  mstack = rbind(m1, m2, m3, m4, m5)
  
  # objective (weighting matrix = identity)
  return(dot(mstack, mstack))
}


# optimization
omega_guess = c(1,2,3,4,5) # obviously change to make this better
optim(omega_guess, gmm_objective, method="BFGS", lower = 0, hessian = TRUE)


# Optim function arguments (for reference):
#optim(par, fn, gr = NULL, …,  method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
#                 "Brent"), lower = 0, upper = Inf, control = list(), hessian = FALSE)
