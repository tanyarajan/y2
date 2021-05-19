#####################################################
#               Metrics Problem Set 3
#                     Tanya Rajan
# Description: This code runs tree code
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
filepath <- "/Users/tanyarajan/Documents/git/y2/spring/metrics/pset3" 
setwd(filepath) 

# set seed
seeder = 123

# read in code file defining functions
#source("code/pset3fns.R")

# setting cores for parallelization
cor<-floor(detectCores(all.tests=FALSE)*.4)
if (is.na(cor) | cor==0){cor <- 1}
registerDoParallel(cores=cor)

#####################################################
#                  Splitting Leaves                 #
#####################################################

# example dataset
y.c <- c(0, 1, 1, 0, 1, 0, 1, 1, 0, 0)
x.c <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
x2.c <- c(2.1, 2.3, 2.5, 2.9, 3.5, 4, 6, 10, 15, 22)
data.d <- cbind(y.c, x.c, x2.c)
data.c <- cbind(y.c, x.c)

# sse function
sse <- function(data, c, k){
  # setup
  p <- dim(data)[2] - 1
  n <- dim(data)[1]
  y <- data[,1]
  x <- matrix(data[,-1], nrow=n)
  
  # sse
  ymin <- sum(y[x[,k] <= c])/length(x[,k][x[,k] <= c])
  yplus <- sum(y[x[,k] > c])/length(x[,k][x[,k] > c])
  if (is.na(yplus)){yplus=0}
  return(sum((y-mean(y))^2) - sum((y - ymin*(x[,k]<=c) - yplus*(x[,k]>c))^2))
}

sse(df, 2, 1)


# growing the tree
tree.grow <- function(data){
  # setup
  p <- dim(data)[2] - 1
  n <- dim(data)[1]
  y <- data[,1]
  x <- matrix(data[,-1], nrow=n)
  
  # running thru partitions
  outmat <- matrix(NA, (n-1), p)
  for (ci in 1:(n-1)){
    for (k in 1:p){
      outmat[ci, k] <- sse(data, x[ci,k], k)
    }
  }
  #print(outmat)
  opt.ids <- which(outmat==max(outmat), arr.ind=T)
  if (dim(opt.ids)[1] > 1){opt.idx <- opt.ids[1,]}
  else{opt.idx <- opt.ids}
  kstar <- opt.idx[2]
  opt.part <- c(kstar, x[opt.idx[1],kstar])
  ypart <- (x[,kstar] <= x[opt.idx[1],kstar])
  
  return(list(split=opt.part, sse=outmat[opt.idx[1]], part=ypart))
}

list[s1, s2, s3] <- tree.grow(df)


# split matrix
#split.c <- as.data.table(matrix(c(1, 0, 1, 8), ncol=4))
#names(split.c) <- c("id", "depth", "kstar", "c")

# assign matrix
#n.c <- dim(data.c)[1]
#assign.c <- as.data.table(cbind(c(2,2,2,2,2,2,2,2,3,3), rep(1,n.c)))
#names(assign.c) <- c("leaf", "depth")


#####################################################
#                 Updating the Tree                 #
#####################################################

# updating the tree
tree.update <- function(data, split, assign, min.size.for.split, max.depth){
  # setup
  p <- dim(data)[2] - 1
  n <- dim(data)[1]
  y <- data[,1]
  x <- matrix(data[,-1], nrow=n)

  # assign matrix
  leaf.sizes <- copy(assign[, .N, by=list(leaf,depth)])

  # end if all leaves are small enough
  #if (sum(leaf.sizes[,3] > min.size.for.split) == 0){ return(list(a=555,b=555))}
  
  # splitting remaining leaves and choosing split with max sse
  out.sse <- matrix(0, ncol=5)
  outvecs <- list()
  cnt <- 1
  for (i in 1:dim(leaf.sizes)[1]){
    if (leaf.sizes[i,3] >= min.size.for.split && leaf.sizes[i,2] < max.depth){
      list[sp, sse, part] <- tree.grow(data[assign[,1]==as.numeric(leaf.sizes[i,1]),])
      out.sse <- rbind(out.sse,c(leaf.sizes$leaf[i], leaf.sizes$depth[i], sp[1], sp[2], sse))
      outvecs[[cnt]] <- part
      cnt = cnt+1
    }
  }
  if (dim(out.sse)[1] == 1){return(list(a=555, b=555))}
  maxsse <- which(out.sse==max(out.sse), arr.ind=T)
  maxrow <- maxsse[1]

  # updating split matrix
  rower <- unlist(out.sse[maxrow, -5])
  split_new <- as.data.table(rbind(as.matrix(split), t(rower)))

  # updating assign matrix
  bool <- outvecs[[(maxrow-1)]]
  assign_new <- copy(assign)
  dp <- rower[2]+1
  assign_new <- assign_new[leaf==rower[1], depth:=dp]
  leaves <- assign_new[leaf==rower[1],]$leaf
  leaves[bool==T] <- rower[1]*2
  leaves[bool==F] <- rower[1]*2 + 1
  assign_new <- assign_new[leaf==rower[1], leaf:=leaves]
  
  # return
  return(list(split=split_new, assign=assign_new))
  
}

#hi<-tree.update(data.c, split.c, assign.c, 3, 3)

#b<-tree.update(data.c, hi$split, hi$assign, 3, 3)

#test<-tree.update(data.c, b$split, b$assign, 3, 3)

#####################################################
#              Growing the Whole Tree               #
#####################################################

tree <- function(data, min.size.for.split=10, max.depth=10){
  # setup
  p <- dim(data)[2] - 1
  n <- dim(data)[1]
  y <- data[,1]
  x <- matrix(data[,-1], nrow=n)
  
  # run tree grow once and generate initial split and assign matrices
  list[sp1, sse1, part1]<-tree.grow(data)
  split <- as.data.table(matrix(c(1, 0, sp1[1], sp1[2]), ncol=4))
  names(split) <- c("leaf", "depth", "kstar", "c")
  init_leaves <- ((part1+1)%%2)+2
  assign <- as.data.table(cbind(init_leaves, rep(1,n)))
  names(assign) <- c("leaf", "depth")
  
  # loop through until we hit end condition
  status <- 1
  assign_new <- assign
  split_new <- split
  counter = 1
  while (split_new != 555){
    assign_old <- assign_new
    split_old <- split_new
    list[split_new, assign_new] <- tree.update(data, split_old, assign_old, min.size.for.split, max.depth)
    #print(split_new)
    #if (counter %% 10 == 0){print(counter)}
    counter = counter+1
  }
  
  # leaves matrix
  leaves <- copy(assign_old[, .N, by=list(leaf,depth)])
  #print(leaves)
  means <- c(rep(0, dim(leaves)[1]))
  for (i in 1:dim(leaves)[1]){
    means[i] <- mean(y[assign_old$leaf == as.numeric(leaves[i,1])], na.rm=T) 
  }
  leaves <- cbind(leaves, means)
  
  
  return(list(assign=assign_old, split=split_old, leaves=leaves))
}

hi<-tree(data.c, min.size.for.split=3, max.depth=3)


#####################################################
#                 Prediction Step                   #
#####################################################

tree.predict <- function(data, tree){
  list[ass, sp, le] <- tree
  dta <- as.data.table(copy(data))
  n <- dim(data)[1]
  x <- as.matrix(dta[,-1], nrow=n)
  dta[,id:=seq(1:n)]
  
  # initial leaves
  xinit <- unlist(sp[[1,3]])
  bool <- ((unlist(x[,xinit]) <= sp[[1,4]]))
  leaf <- ((bool+1)%%2)+2
  init_depth <- as.matrix(rep(1,n), ncol=1)
  data.assign <- as.data.table(cbind(seq(1:n),leaf))
  names(data.assign) <- c("id", "leaf")

  for (i in 2:dim(sp)[1]){
    xidx <- sp[[i,3]]
    thresh <- sp[[i,4]]
    old.leaf <- sp[[i,1]]
    data.assign$leaf[(x[,xidx] <= thresh & data.assign$leaf==old.leaf)] <- old.leaf*2 
    data.assign$leaf[(x[,xidx] > thresh & data.assign$leaf==old.leaf)] <- old.leaf*2 + 1
  }

  out<-merge(data.assign, le, by=c("leaf"), all.x=T)
  return(out[order(id)]$means)
  
}

tree.predict(data.c, hi)


#####################################################
#                 DGP Generation                    #
#####################################################

# data generation
setN <- 1000
data.gen<-mvrnorm(setN, rep(0,3), eye(3))

y <- apply(data.gen, 1, function(x){3*min(x[1], x[2]) + x[3]})
y2 <- apply(data.gen, 1, function(x){3*x[1] - 3*x[2] + x[3]})
x <- data.gen[,-3]
df <- cbind(y,x)
df2 <- cbind(y2, x)


# ols and tree on full data
mse.ols<-mean((y - predict.lm(lm(y ~ x - 1)))^2)
options(warn = -1)
mse.tree<-mean((y-tree.predict(df, tree(df)))^2)



#####################################################
#                   10-Fold CV                      #
#####################################################
# split into 10 parts
dt<-as.data.table(df)
names(dt) <- c("y", "x1", "x2")
dt[,group:=sample.int(10, setN, replace=T)]

dt2<-as.data.table(df2)
names(dt2) <- c("y", "x1", "x2")
dt2[,group:=sample.int(10, setN, replace=T)]

# testing and training function
test.train <- function(data, g, min.size.for.split=10, max.depth=10){
  train <- as.matrix(copy(dt[group != g,][,group:=NULL]))
  #print(train)
  test <- as.matrix(copy(dt[group == g][,group:=NULL]))
  train.tree <- tree(train, min.size.for.split, max.depth)
  train.ols <- lm(train[,1] ~ train[,2] + train[,3]  - 1)
  test.pred.tree <- tree.predict(test, train.tree)
  test.pred.ols <- predict.lm(train.ols)
  mse.tree <- mean((y-test.pred.tree)^2)
  mse.ols <- mean((y-test.pred.ols)^2)
  print(paste0("done with ", g))
  return(list(mse.tree=mse.tree, mse.ols=mse.ols))
}

# answers to part e
xval.e<-lapply(seq(1:10), function(x){test.train(dt,x)})
xval.e.out <- unlist(xval.e)
xval.e.ols <- mean(xval.e.out[names(xval.e.out) == "mse.ols"])
xval.e.tree <- mean(xval.e.out[names(xval.e.out) == "mse.tree"])


# answers to part f
xval.f<-lapply(seq(1:10), function(x){test.train(dt2,x)})
xval.f.out <- unlist(xval.f)
xval.f.ols <- mean(xval.f.out[names(xval.f.out) == "mse.ols"])
xval.f.tree <- mean(xval.f.out[names(xval.f.out) == "mse.tree"])

