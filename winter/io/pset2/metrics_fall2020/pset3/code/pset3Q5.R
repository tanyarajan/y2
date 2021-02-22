#####################################################
#               Metrics Problem Set 3
#                      Question 5
#                     Tanya Rajan
# Description: This file uses data from Angrist and
# Evans (1998) to calculate MTRs from 5 different 
# specifications. It then uses the MTRs to compute
# and plot MTEs over U (with X fixed at its mean).
# Finally, it uses MTRs to calculate target parameters
# of interest: ATE, ATT, ATUT, and LATE
#####################################################

# setup
rm(list=ls())
require(data.table)
require(ggplot2)
require(pracma)
require(doParallel)
require(kableExtra)
require(maxLik)
require(ggpubr)
require(future.apply)

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
#                Setting up Data                    #
#####################################################
# read in data
datain <- as.data.table(read.csv("data/angrist_evans_clean.csv"))

# variables of interest
Y <- as.matrix(datain$worked, ncol=1)
D <- as.matrix(datain$more2kids, ncol=1)
X <- as.matrix(datain[,list(age, ageat1st, agekid1, agekid2, boy1st, 
                            boy2nd, black, hispanic, otherrace)])
Z <- list()
Z[[1]] <- datain$samesex
Z[[2]] <- datain$twins
Z[[3]] <- cbind(datain$samesex, datain$twins)



D00 <- datain[samesex==0 & twins==0, more2kids]
D01 <- datain[samesex==0 & twins==1, more2kids]
D10 <- datain[samesex==1 & twins==0, more2kids]
D11 <- datain[samesex==1 & twins==1, more2kids]


#####################################################
#           Defining MTR Specifications             #
#####################################################
ED <- length(D[which(D==1)])/length(D)
spec1 <- function(p1, p0, u, x, xu){
  # m1
  a1 = p1[1]
  b1 = 2*p1[2]
  g1 = p1[3:(2+ncol(x))]
  m1 <- a1 + b1*u + x%*%g1
  
  # m0
  a0 = p0[1] - p0[2]
  b0 = 2*p0[2]
  g0 = p0[3:(2+ncol(x))]
  m0 <- a0 + b0*u + x%*%g0
  
  # return
  return(list(mtr1 = m1, mtr0 = m0))
}


spec2 <- function(p1, p0, u, x, xu){
  
  # m1 
  a1 = p1[2] + p1[1]
  b1 = 2*(p1[3] + p1[4])
  g = p1[5:(4+ncol(x))] 
  m1 <- a1 + b1*u + x%*%g
  
  # m2
  a0 = p1[1]
  b0 = 2*p1[3]
  m0 <- a0 + b0*u + x%*%g
  
  # return
  return(list(mtr1 = m1, mtr0 = m0))
}

spec3 <- function(p1, p0, u, x, xu){
  # m1
  a1 = p1[1]
  b1 = 2*p1[2]
  g1 = p1[3:(2+ncol(x))]
  d1 = 2*p1[(3+ncol(x)):(2+2*ncol(x))]
  m1 <- a1 + b1*u + x%*%g1 + xu%*%d1
  
  # m0
  a0 = p0[1] - p0[2]
  b0 = 2*p0[2]
  g0 = p0[3:(2+ncol(x))] - p0[(3+ncol(x)):(2+2*ncol(x))]
  d0 = 2*p0[(3+ncol(x)):(2+2*ncol(x))]
  m0 <- a0 + b0*u + x%*%g0 + xu%*%d0
  
  # return
  return(list(mtr1 = m1, mtr0 = m0))
}

spec4 <- function(p1, p0, u, x, xu){
  u2 = u^2
  
  # m1
  a1 = p1[1]
  b11 = 2*p1[2]
  b12 = 3*p1[3]
  g1 = p1[4:(3+ncol(x))]
  m1 <- a1 + b11*u + b12*u2 + x%*%g1
  
  # m0
  a0 = p0[1] - p0[2]
  b01 = 2*(p0[2] - p0[3])
  b02 = 3*p0[3]
  g0 = p0[4:(3+ncol(x))]
  m0 <- a0 + b01*u + b02*u2 + x%*%g0
  
  #return
  return(list(mtr1 = m1, mtr0 = m0))
}

spec5 <- function(p1, p0, u, x, xu){
  u2 = u^2
  u3 = u^3
  
  # m1
  a1 = p1[1]
  b11 = 2*p1[2]
  b12 = 3*p1[3]
  b13 = 4*p1[4]
  g1 = p1[5:(4+ncol(x))]
  m1 <- a1 + b11*u + b12*u2 + b13*u3 + x%*%g1
  
  # m0
  a0 = p0[1] - p0[2]
  b01 = 2*(p0[2] - p0[3])
  b02 = 3*(p0[3] - p0[4])
  b03 = 4*p0[4]
  g0 = p0[5:(4+ncol(x))]
  m0 <- a0 + b01*u + b02*u2 + b03*u3 + x%*%g0
  
  #return
  return(list(mtr1 = m1, mtr0 = m0))
}



#####################################################
#                Estimating MTRs                    #
#####################################################
#### Finding MTR coefficients ####
mtr_coeffs <- function(y, d, x, z){
  # calculate propensity score (using own logit function)
  #Zlog <- cbind(1, z, x)
  #theta.g <- rep(.01, ncol(Zlog))
  #logit <- maxLik(fnL, start = theta.g, yvec = d, xvec=cbind(Zlog))
  #p <- exp(Zlog%*%logit$estimate)/(1 + exp(Zlog%*%logit$estimate))
  
  # calculating using R's function to speed it up for now
  logit <- glm(d ~ z + x, family = "binomial")
  p <- exp(cbind(1, z, x)%*%logit$coefficients)/(1 + exp(cbind(1, z,x)%*%logit$coefficients))
  
  # preparing variables
  xp <- x*repmat(p, 1, ncol(x))
  p2 <- p^2
  p3 <- p^3
  
  # specifications
  reggr <- list()
  reggr[[1]] <- cbind(p, x)
  reggr[[2]] <- cbind(p, x) # need to fix
  reggr[[3]] <- cbind(p, x, xp)
  reggr[[4]] <- cbind(p, p2, x)
  reggr[[5]] <- cbind(p, p2, p3, x)
  
  # running the regressions
  pars1 <- list()
  pars0 <- list()
  u1 <- seq(1:length(y[which(d==1)]))
  u0 <- seq(1:length(y[which(d==0)]))
  for (i in 1:5){
    dep1<- y[which(d==1)]
    indep1 <- reggr[[i]][which(d==1),]
    dep0 <- y[which(d==0)]
    indep0 <- reggr[[i]][which(d==0),]
    pars1[[i]] <- reg(dep1, indep1)$b
    pars0[[i]] <- reg(dep0, indep0)$b
    if (i==2){ # run saturated specification
      pd <- p*d
      pars1[[i]] <- reg(y, cbind(d, p, pd, x))$b
      pars0[[i]] <- pars1[[i]]
    }
  }
  
  # return
  return(list(p1 = pars1, p0 = pars0, pz=p, pzcoeff=logit$coefficients))
  
}

#### Calculating MTRs over UD ####
mtr <- function(u, x, pars1, pars0){
  # feeding uniform u and evaluating at mean of x
  xval <- repmat(x, nrow(u), 1)
  xu <-  xval*repmat(u, 1, length(x))
  
  # looping through specifications and saving output
  mtr = array(dim=c(nrow(u), 2, 5))
  mte = array(dim=c(nrow(u), 1, 5))
  for (i in 1:5){
    out <- callit("spec",i)(pars1[[i]], pars0[[i]], u, xval, xu)
    mtr[,1,i] <- out$mtr0
    mtr[,2,i] <- out$mtr1
    mte[,1,i] <- out$mtr1 - out$mtr0
  }
  return(list(mtr = mtr, mte=mte))
}


#### Graphing MTE results ####
graph.mte <- function(i, mtrs, u, min, max){
  # creating dataset with all mtrs
  plotdta <- as.data.table(cbind(mtrs, u))
  names(plotdta) <- c("MTR0", "MTR1", "u")
  
  # calculating MTEs
  plotdta[,`:=`(MTE = MTR1-MTR0)]
  
  # plotting
  plotdta <- plotdta[order(u)]
  pL <- melt(plotdta, id.vars="u")
  lt <- c("dashed", "dotted", "solid")
  out <- ggplot(aes(y=value, x= u, linetype=variable), data=pL) + 
    geom_line() + geom_line() + theme_bw() + ylim(min , max) + 
    ggtitle(paste0("Specification ", i)) + ylab("MTE") + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_linetype_manual(values=lt)
  
  # return
  return(out)
}



#### Computing Target parameters ####
te_mtr <- function(mtrs, pzu, pziv, spec){
  # creating dataset with all mtrs
  pdata <- as.data.table(cbind(mtrs[,,1], mtrs[,,2], mtrs[,,3], mtrs[,,4], mtrs[,,5]))
  names(pdata) <- c("mtr10", "mtr11", "mtr20", "mtr21", "mtr30", "mtr31", 
                    "mtr40", "mtr41", "mtr50", "mtr51")
  
  # calculating MTEs
  pdata[,`:=`(mte1=mtr11 - mtr10, mte2=mtr21 - mtr20, mte3 = mtr31 - mtr30, 
              mte4=mtr41-mtr40, mte5 = mtr51-mtr50)]
  
  # calculating proportions for ATU/ATT weights
  pd1 = length(D[which(D == 1)])/length(D)
  pd0 = length(D[which(D == 0)])/length(D)
  
  # ATE, ATU, ATT weights
  pdata[,atew:=1]
  pdata[,plessu:=pzu]
  pdata[,pgreatu:=1-plessu]
  pdata[,atuw:=plessu/pd0]
  pdata[,attw:=pgreatu/pd1]
  if (spec != 3){pdata[,latew:=pziv]}
  if (spec == 3){pdata[,latew1:=pziv[[1]]][,latew2:=pziv[[2]]]} # 2 instruments
  
  # calculating parameters per person
  pdata[,`:=`(ate1=mte1*atew, ate2=mte2*atew, ate3=mte3*atew, ate4=mte4*atew, ate5=mte5*atew)]
  pdata[,`:=`(atu1=mte1*atuw, atu2=mte2*atuw, atu3=mte3*atuw, atu4=mte4*atuw, atu5=mte5*atuw)]
  pdata[,`:=`(att1=mte1*attw, att2=mte2*attw, att3=mte3*attw, att4=mte4*attw, att5=mte5*attw)]
  if (spec != 3){
    pdata[,`:=`(late1=mte1*latew, late2=mte2*latew, late3=mte3*latew, late4=mte4*latew, late5=mte5*latew)]
  }
  if (spec == 3){
    pdata[,`:=`(late11=mte1*latew1, late21=mte2*latew1, late31=mte3*latew1, 
                late41=mte4*latew1, late51=mte5*latew1)]
    pdata[,`:=`(late12=mte1*latew2, late22=mte2*latew2, late32=mte3*latew2, 
                late42=mte4*latew2, late52=mte5*latew2)]
    w.inst1<-reg(D, Z[[3]])$b[2,]
    w.inst2<-reg(D, Z[[3]])$b[3,]
    psi <- w.inst1/(w.inst1 + w.inst2)
    pdata[,`:=`(late1=late11*pi + late12*(1-pi), late2=late21*pi + late22*(1-pi), 
                late3=late31*pi + late32*(1-pi), late4=late41*pi + late42*(1-pi),
                late5=late51*pi + late52*(1-pi))]
  }
  
  # taking means
  ATE  <- colMeans(cbind(pdata$ate1, pdata$ate2, pdata$ate3, pdata$ate4, pdata$ate5))
  ATU  <- colMeans(cbind(pdata$atu1, pdata$atu2, pdata$atu3, pdata$atu4, pdata$atu5))
  ATT  <- colMeans(cbind(pdata$att1, pdata$att2, pdata$att3, pdata$att4, pdata$att5))
  LATE <- colMeans(cbind(pdata$late1, pdata$late2, pdata$late3, pdata$late4, pdata$late5))
  
  #return(pdata)
  output <- rbind(ATE, ATU, ATT, LATE)
  colnames(output) <-  c("Spec. 1","Spec. 2","Spec. 3","Spec. 4","Spec. 5")
  return(output)
}


#####################################################
#                 Running all Code                  #
#####################################################
# Setting up ud
Nu <- 100
ud <- as.matrix(seq(1:Nu)/Nu, ncol=1)

run_mtr <- function(ud, y,d,x,z,spec){
  # getting parameters
  mtrpar <- mtr_coeffs(y,d,x,z)
  
  # calculating mtr over ud with x at its mean
  xbar <- colMeans(x)
  mtrout <- mtr(ud, xbar, mtrpar$p1, mtrpar$p0)
  mte <- mtrout$mte
  
  # calculating pz for x at its mean
  xrep <- repmat(xbar, length(y),1)
  pzcoeff <- mtrpar$pzcoeff
  pz <- exp(cbind(1, z, xrep)%*%pzcoeff)/(1 + exp(cbind(1, z,xrep)%*%pzcoeff))
  
  # weights for estimators
  plan(multicore)
  pzu.comp<-future_apply(ud,1,function(u){length(which(pz<u))/length(pz)})
  
  # iv weights
  pzbounds <- c(NA, NA)
  if (is.null(ncol(z))){
    pzbounds[1] <- mean(pz[z==0])
    pzbounds[2] <- mean(pz[z==1])
    pziv.comp <- (ud < pzbounds[2] & ud > pzbounds[1])/(pzbounds[2] - pzbounds[1])
  }
  if (!is.null(ncol(z))){ # 2sls case
    pzbounds[1] <- mean(pz[z[,1]==0])
    pzbounds[2] <- mean(pz[z[,1]==1])
    pzbounds[3] <- mean(pz[z[,2]==0])
    pzbounds[4] <- mean(pz[z[,2]==1])
    pziv1 <- (ud < pzbounds[2] & ud > pzbounds[1])/(pzbounds[2] - pzbounds[1])
    pziv2 <- (ud < pzbounds[4] & ud > pzbounds[3])/(pzbounds[4] - pzbounds[3])
    pziv.comp <- c(pziv1, pziv2)
  }
  
  # treatment effects
  finalout <-te_mtr(mtrout$mtr, pzu.comp, pziv.comp, spec)
  
  # graphing
  mtes <-cbind(mtrout$mte[,,1],mtrout$mte[,,2],mtrout$mte[,,3],mtrout$mte[,,4],
               mtrout$mte[,,5],mtrout$mtr[,,1],mtrout$mtr[,,2],mtrout$mtr[,,3],
               mtrout$mtr[,,4],mtrout$mtr[,,5])
  ymin <- min(mtes)
  ymax <- max(mtes)
  g    <- lapply(1:5, function(i){graph.mte(i, mtrout$mtr[,,i], ud, ymin, ymax)})
  mte.inst <- ggarrange(g[[1]],g[[2]],g[[3]],g[[4]],g[[5]], nrow=3, ncol=2)
  ggsave(mte.inst, width = 9, height = 11, file = paste0("figures/tr_mte_inst",spec,".png"))
  
  return(finalout)
}

# write to tex
tes1 <- run_mtr(ud, Y, D, X, Z[[1]], 1)
kable(tes1, digits=4, format="latex", booktabs=TRUE) %>%
  write(file="tables/tr_tes1.tex")

tes2 <- run_mtr(ud, Y, D, X, Z[[2]], 2)
kable(tes2, digits=4, format="latex", booktabs=TRUE) %>%
  write(file="tables/tr_tes2.tex")

tes3 <- run_mtr(ud, Y, D, X, Z[[3]], 3)
kable(tes2, digits=4, format="latex", booktabs=TRUE) %>%
  write(file="tables/tr_tes3.tex")


# 2SLS estimates
iv(Y, cbind(D,X), cbind(Z[[1]], X))$b[2,]
iv(Y, cbind(D,X), cbind(Z[[1]], X))$se[2]
iv(Y, cbind(D,X), cbind(Z[[2]], X))$b[2,]
iv(Y, cbind(D,X), cbind(Z[[2]], X))$se[2]
iv(Y, cbind(D,X), cbind(Z[[3]], X))$b[2,]
iv(Y, cbind(D,X), cbind(Z[[3]], X))$se[2]

