library(MASS)
library(rpart)
library(glmnet)
library(rpart.plot)
library(partykit)
library(intervals)
library(ranger)

fname <- file('stdin')
open(fname)
n = as.numeric(readLines(fname, n=1)) # {200,400,800,1600}
transf = as.numeric(readLines(fname, n=1)) #{1,2}
etasetting = as.numeric(readLines(fname, n=1)) # {0,1,2}
sim.round = as.numeric(readLines(fname, n=1)) # {1,2,3,4,5}
close(fname)

Xsetting = 1
if(Xsetting == 1){
  testdata <- readRDS('src/rf1128-testdata-unif.rds')
}
if(Xsetting == 2){
  testdata <- readRDS('src/rf1128-testdata-gaus.rds')
}
n0 = testdata$n0
X0 = testdata$X0

p = 6

filename = paste0('Treatment1204-eXEst-MethodRepro-Rpart-n=',n,'-(A)IPW=',transf,'-etasetting=',etasetting,'-simround=',sim.round,'.RDS')
print(filename)

kappa <- function(X){
  # X is the n by p matrix
  te.effect <- function(x){
    # x is the vector here
    if(x[3]<=0){
      if((x[1]<=0) & (x[2]<=0)) return(-2.4)
      if((x[1]<=0) & (x[2]>0)) return(-3.8)
      if((x[1]>0) & (x[2]<=0)) return(-1.2)
      if((x[1]>0) & (x[2]>0)) return(-0)
    }
    else if(x[3]>0){
      if((x[4]<=0) & (x[5]<=0)) return(2.5)
      if((x[4]<=0) & (x[5]>0)) return(3.5)
      if((x[4]>0) & (x[5]<=0)) return(4.6)
      if((x[4]>0) & (x[5]>0)) return(1.4)
    }
  }
  return(apply(X, MARGIN=1, FUN=te.effect))
}
eta <- function(X){
  base.effect <- function(x){
    if(etasetting==0){
      return(0)
    }
    else if(etasetting==1){
      return(0.5*x[1]+x[2])
    }
    else if(etasetting==2){
      return(0.5*x[1]^2 + x[2])
    }
  }
  return(apply(X, MARGIN = 1, FUN=base.effect))
}

b <- c(c(-2,-1,0,1,2)/5, 0)
e <- function(X, b){
  val = exp(X%*%b)/(1+ exp(X%*%b))
  return(as.vector(val))
}

nsim = 20
save.info = list(NA)
for(i.sim in 1:nsim){
  print(paste0('i.sim=',i.sim,'/',nsim))
  set.seed((sim.round-1)*nsim+i.sim)
  if(Xsetting == 1){
    X = matrix(runif(n*p, min=-1,max=1), nrow=n, ncol=p)
  }else if(Xsetting == 2){
    X = mvrnorm(n, mu=rep(0, 6), diag(p))
  }
  colnames(X) = paste0('X',1:p)
  eX = e(X,b)
  W = rbinom(n, size=1, eX)
  Y = eta(X) + (1/2)*(2*W-1)*kappa(X) + rnorm(n, sd=0.1)
  Y.pure = eta(X) + (1/2)*(2*W-1)*kappa(X)
  # eX estimation
  fit = glm(W~X, family = 'binomial')
  eXhat = predict(fit, type='response')
  
  # mu_0 and mu_1 estimate using cross-fitting
  copyid = sample(1:2, size=n, replace=T)
  mu0 = mu1 = rep(NA, n)
  for(copy in copyid){
    map0 = (W==0)&(copyid==copy)
    data0 = data.frame(cbind(X[map0,], Y[map0])); colnames(data0) = c(colnames(X),'Y')
    fmla0 = as.formula(paste('Y ~ ',paste(colnames(X), collapse = '+')))
    rf0 = ranger(fmla0, data=data0)
    map0.pred = (copyid==(3-copy)) # 1-> 2; 2->1
    mu0[map0.pred] = predict(rf0, data=X[map0.pred,])$predictions
    
    map1 = (W==1)&(copyid==copy)
    data1 = data.frame(cbind(X[map1,], Y[map1])); colnames(data1) = c(colnames(X),'Y')
    fmla1 = as.formula(paste('Y ~ ',paste(colnames(X), collapse = '+')))
    rf1 = ranger(fmla1, data=data1)
    map1.pred = (copyid==(3-copy)) # 1-> 2; 2->1
    mu1[map1.pred] = predict(rf1, data=X[map1.pred,])$predictions
  }
  etahat = (mu0 + mu1)/2
  
  if(transf==1){
    # propensity score transformed
    Y_prox = W*Y/eXhat - (1-W)*Y/(1-eXhat)
  }else if(transf==2){
    # AIPW
    Y_prox = W/eXhat*(Y-mu1) - (1-W)/(1-eXhat)*(Y-mu0) + (mu1 - mu0)
  }
  
  
  # repro -------------------------------------------------------------------
  M = 500
  repro.list = repro.list.small = list(NA)
  for(m in 1:M){
    if(m%% 100==1) print(paste0('repro-',m,'/',M))
    eps.repro = rnorm(n, sd=0.1)
    Y.repro = Y - eps.repro
    if(transf==1){
      Y_prox.repro = W*Y.repro/eXhat - (1-W)*Y.repro/(1-eXhat)
    }else if(transf==2){
      Y_prox.repro = W/eXhat*(Y.repro-mu1) - (1-W)/(1-eXhat)*(Y.repro-mu0) + (mu1 - mu0)
    }
    
    data.repro = as.data.frame(cbind(X, Y_prox.repro))
    colnames(data.repro) = c(colnames(X),'Y')
    fmla = as.formula(paste('Y~', paste(colnames(X), collapse = '+')))
    tree.repro = rpart(fmla, data=data.repro, control=rpart.control(cp=0.001))
    opcp = tree.repro$cptable[,1][which.min(tree.repro$cptable[,4])]
    tree.repro = prune(tree.repro, cp=opcp)
    
    term.nodes = unique(tree.repro$where)
    mean.vec = sd.vec = mean.rob.vec = lower.rob.vec = upper.rob.vec = rep(NA, length(term.nodes))
    names(mean.vec) = names(sd.vec) = names(mean.rob.vec) = names(lower.rob.vec) = names(upper.rob.vec) = term.nodes
    for(node in term.nodes){
      map = (tree.repro$where == node)
      y = Y_prox.repro[map]
      y.mean = mean(y)
      y.sd = sqrt(var(y)/(length(y)))
      mean.vec[as.character(node)] = y.mean
      sd.vec[as.character(node)] = y.sd
      
      # robust estimation
      map0 = (tree.repro$where == node)&(W==0)
      temp0 = Y.repro[map0] - etahat[map0]
      # temp0 = Y[map0] - eta(X[map0,])
      map1 = (tree.repro$where == node)&(W==1)
      temp1 = Y.repro[map1] - etahat[map1]
      # temp1 = Y[map1] - eta(X[map1,])
      thres = qbinom(c(0.025, 0.975), size=length(y), prob=1/2)
      temp = c(temp0*(-2), temp1*2)
      y.mean.rob = median(temp)
      y.ci.rob = quantile(temp, probs = thres/length(y))
      
      # # robust estimation - 2
      # map0 = (tree.repro$where == node)&(W==0)
      # temp0 = Y.repro[map0] - etahat[map0]
      # thres0 = qbinom(c(0.025, 0.975), size=sum(map0), prob=1/2)
      # y.ci.rob0 = quantile(temp0*(-2), probs = thres0/sum(map0))
      # 
      # map1 = (tree.repro$where == node)&(W==1)
      # temp1 = Y.repro[map1] - etahat[map1]
      # thres1 = qbinom(c(0.025, 0.975), size=sum(map1), prob=1/2)
      # y.ci.rob1 = quantile(temp1*2, probs = thres1/sum(map1))
      # 
      # temp = c(temp0*(-2), temp1*2)
      # y.mean.rob = median(temp)
      # y.ci.rob = c(min(y.ci.rob0[1],y.ci.rob1[1]),
      #              max(y.ci.rob0[2],y.ci.rob1[2]))
      
      mean.rob.vec[as.character(node)] = y.mean.rob
      lower.rob.vec[as.character(node)] = y.ci.rob[1]
      upper.rob.vec[as.character(node)] = y.ci.rob[2]
    }
    repro.list[[m]] = list(tree.repro = tree.repro,
                           term.nodes = term.nodes,
                           mean.vec = mean.vec,
                           sd.vec = sd.vec,
                           mean.rob.vec = mean.rob.vec,
                           lower.rob.vec = lower.rob.vec,
                           upper.rob.vec = upper.rob.vec)
    repro.list.small[[m]] = list(term.nodes = term.nodes,
                                 mean.vec = mean.vec,
                                 sd.vec = sd.vec,
                                 mean.rob.vec = mean.rob.vec,
                                 lower.rob.vec = lower.rob.vec,
                                 upper.rob.vec = upper.rob.vec)
  }
  
  
  # evaluation - asymp & rob------------------------------------------------------
  ## evaluation on the test data observations ##
  Y0_prox.pred.mat = Y0_prox.lower.mat = Y0_prox.upper.mat = matrix(NA, nrow=n0, ncol=M)
  Y0_prox.rob.pred.mat = Y0_prox.rob.lower.mat = Y0_prox.rob.upper.mat = matrix(NA, nrow=n0, ncol=M)
  for(m in 1:M){
    tree.repro = repro.list[[m]]$tree.repro
    mean.vec = repro.list[[m]]$mean.vec
    sd.vec = repro.list[[m]]$sd.vec
    mean.rob.vec = repro.list[[m]]$mean.rob.vec
    lower.rob.vec = repro.list[[m]]$lower.rob.vec
    upper.rob.vec = repro.list[[m]]$upper.rob.vec
    
    y_prox.pred.node = predict(as.party(tree.repro), newdata = data.frame(X0), type='node')
    # standard asymp
    y_prox.pred.val = (sapply(y_prox.pred.node, FUN=function(x) mean.vec[as.character(x)]))
    y_prox.pred.sd = (sapply(y_prox.pred.node, FUN=function(x) sd.vec[as.character(x)]))
    Y0_prox.pred.mat[,m] = y_prox.pred.val
    Y0_prox.lower.mat[,m] = y_prox.pred.val - qnorm(0.975)*y_prox.pred.sd
    Y0_prox.upper.mat[,m] = y_prox.pred.val + qnorm(0.975)*y_prox.pred.sd
    
    # robust
    Y0_prox.rob.pred.mat[,m] = (sapply(y_prox.pred.node, FUN=function(x) mean.rob.vec[as.character(x)]))
    Y0_prox.rob.lower.mat[,m] = (sapply(y_prox.pred.node, FUN=function(x) lower.rob.vec[as.character(x)]))
    Y0_prox.rob.upper.mat[,m] = (sapply(y_prox.pred.node, FUN=function(x) upper.rob.vec[as.character(x)]))
  }
  # summary - asymp
  summ.Y0_prox = matrix(NA, nrow=n0, ncol=3); colnames(summ.Y0_prox) = c('val','cov','len')
  Y0_prox.single.CI = matrix(NA, nrow=M, ncol=2)
  for(i0 in 1:n0){
    if(i0%% 200==1) print(paste0('interval union -',i0,'/',n0))
    summ.Y0_prox[i0,1] = mean(Y0_prox.pred.mat[i0,])
    for(m in 1:M){
      Y0_prox.single.CI[m,] = c(Y0_prox.lower.mat[i0,m], Y0_prox.upper.mat[i0,m])
    }
    Y0_prox.single.union = as.matrix(interval_union(Intervals(Y0_prox.single.CI)))
    summ.Y0_prox[i0,2] = sum((Y0_prox.single.union[,2] >= kappa(X0)[i0])*
                               (Y0_prox.single.union[,1] <= kappa(X0)[i0]))
    summ.Y0_prox[i0,3] = sum(Y0_prox.single.union[,2] - Y0_prox.single.union[,1])
  }
  # summary - robust
  summ.Y0_prox.rob = matrix(NA, nrow=n0, ncol=3); colnames(summ.Y0_prox.rob) = c('val','cov','len')
  Y0_prox.rob.single.CI = matrix(NA, nrow=M, ncol=2)
  for(i0 in 1:n0){
    if(i0%% 200==1) print(paste0('robust interval union -',i0,'/',n0))
    summ.Y0_prox.rob[i0,1] = mean(Y0_prox.rob.pred.mat[i0,])
    for(m in 1:M){
      Y0_prox.rob.single.CI[m,] = c(Y0_prox.rob.lower.mat[i0,m], Y0_prox.rob.upper.mat[i0,m])
    }
    Y0_prox.rob.single.union = as.matrix(interval_union(Intervals(Y0_prox.rob.single.CI)))
    summ.Y0_prox.rob[i0,2] = sum((Y0_prox.rob.single.union[,2] >= kappa(X0)[i0])*
                                   (Y0_prox.rob.single.union[,1] <= kappa(X0)[i0]))
    summ.Y0_prox.rob[i0,3] = sum(Y0_prox.rob.single.union[,2] - Y0_prox.rob.single.union[,1])
  }
  
  save.info[[i.sim]] <- list(repro.list=repro.list.small,
                             summ.Y0_prox=summ.Y0_prox,
                             summ.Y0_prox.rob=summ.Y0_prox.rob)
}
saveRDS(save.info, filename)
