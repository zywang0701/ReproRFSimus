# Partitioned Linear Models -----------------------------------------------
# Y_k = X_k^T beta_k + sigma_k N(0,1) 
# Compare the following methods:
# 1. truth: the test's data group labels are known 
# 2. post: MOB function $Y~X$ without noise inside
# 3. oracle: MOB function $Y~X+eps0$ with truth noise
# 4. repro: MOB function $Y~X+eps$ with repro noise
# The original code is in simu0119
# we also consider the filtering 5\% in simu0205
library(rpart)
library(MASS)
library(intervals)
library(party)
source('helper.R')


## Partitioning --------------------------------------------------------------
if(group.setting == 1){
  get.group = function(x){
    if((x[1]<=0) & (x[2]<=0)) return(1)
    if((x[1]<=0) & (x[2]>0)) return(2)
    if((x[1]>0) & (x[2]<=0)) return(3)
    if((x[1]>0) & (x[2]>0)) return(4)
  }
}else if(group.setting == 2){
  get.group = function(x){
    if(x[1]<=-0.2){
      if(x[2]<=-0.1) return(1)
      if(x[2]>-0.1) return(2)
    }
    if(x[1]>-0.2){
      if(x[2]<=0.1) return(3)
      if(x[2]>0.1) return(4)
    }
  }
}else if(group.setting == 3){
  get.group = function(x){
    if(x[1]<=-0.5){
      return(1)
    }
    if(x[1]>-0.5){
      if(x[2]<=0) return(2)
      if(x[2]>0){
        if(x[3]<=0) return(3)
        if(x[3]>0) return(4)
      }
    }
  }
}


## Coef in each partition --------------------------------------------------
if(coef.setting==1){
  coef1 = rep(1,6)
  coef2 = coef1 + rep(0.1, 6)
  coef3 = coef1 + rep(-0.5,6)
  coef4 = coef1 + rep(0.15,6)
}
if(coef.setting==2){
  coef1 = seq(1,6)/5
  coef2 = coef1 + rep(0.1, 6)
  coef3 = coef1 - seq(0,5)/10
  coef4 = coef1 + c(0.1, 0.2, 0.3, -0.1, -0.2, -0.3)
}
if(coef.setting==3){
  coef1 = rep(1,6)
  coef2 = seq(0,2.5,by=0.5)
  coef3 = -rep(0.5,6)
  coef4 = c(0, 0.5, 1, -0.5, -1, 1)
}
B = rbind(coef1,coef2,coef3,coef4)

## Noise Level in each partition -------------------------------------------
if(noi.str.setting == 1){
  noi.str = c(2,1.5,0.5,0.8)
}else if(noi.str.setting == 2){
  noi.str = c(1.5, 1.2, 0.9, 0.5)
}else if(noi.str.setting == 3){
  noi.str = c(1.2, 1.1, 1, 0.8)
}


## Setup ------------------------------------------------------
L = 4
n = 1000
n.test = 2000
p = 5
M = 500 # repro times
group.setting = 1
coef.setting = 1
noi.str.setting = 1
# covariates X covariance
Sigma = diag(p)

## Main Code ---------------------------------------------------------------
nsim = 100
# info.pred.list: results for predictive inference
# info.cond.list: results for conditional inference
# uni_nodes.list: results for partitions
info.pred.list = info.cond.list = uni_nodes.list = rep(list(NA), nsim)
for(i.sim in 1:nsim){
  print(i.sim)
  ## training 
  X = mvrnorm(n, mu=rep(0,p), Sigma)
  colnames(X) = paste0('X',1:p,seq="")
  group = apply(X, 1, get.group) # the group label for each observation Xi
  Xint = cbind(1, X)
  Y_star = rowSums(Xint*B[group,])
  noise = rep(NA, n)
  for(i in 1:n){
    noise[i] = rnorm(1, sd=noi.str[group[i]])
  }
  Y = Y_star + noise
  ## test
  X.test=mvrnorm(n.test, mu=rep(0,p), Sigma)
  colnames(X.test)=paste0('X',1:p,seq="")
  group.test = apply(X.test, 1, get.group)
  Xint.test = cbind(1, X.test)
  Y_star.test = rowSums(Xint.test * B[group.test,])
  noise.test = rep(NA, n.test)
  for(i in 1:n.test){
    noise.test[i] = rnorm(1, sd=noi.str[group.test[i]])
  }
  Y.test = Y_star.test + noise.test
  
  ## truth
  out.truth = helper(group, group.test, X, Y, X.test)
  ######## post ##########
  data1 = as.data.frame(cbind(X, Y)); colnames(data1) = c(colnames(X),"Y")
  tree.post <- mob(Y ~ X1+X2+X3+X4+X5 | X1+X2+X3+X4+X5,
                   data = data1, model  = linearModel)
  train.nodes.post = predict(tree.post, as.data.frame(X), type='node')
  test.nodes.post = predict(tree.post, as.data.frame(X.test), type='node')
  out.post = helper(train.nodes.post, test.nodes.post, X, Y, X.test)
  ######## oracle ########
  data0 = as.data.frame(cbind(X, Y, noise)); colnames(data0) = c(colnames(X), "Y", "eps")
  tree0 = mob(Y ~ X1+X2+X3+X4+X5+eps | X1+X2+X3+X4+X5,
              data = data0, model  = linearModel)
  plot(tree0)
  train.nodes0 = predict(tree0, as.data.frame(cbind(X, eps=noise)), type='node')
  test.nodes0 = predict(tree0, as.data.frame(cbind(X.test, eps=noise.test)), type='node')
  out0 = helper(train.nodes0, test.nodes0, X, Y, X.test)
  ######## repro #######
  ci.gen.low = ci.gen.up = ci.cond.gen.low = ci.cond.gen.up = matrix(NA, nrow=n.test, ncol=M)
  uni_nodes.vec = rep(NA, M)
  for(m in 1:M){
    print(paste0('repro-',m,'/M'))
    noise.gen = rnorm(n)
    data.gen = as.data.frame(cbind(X, Y, noise.gen))
    colnames(data.gen) = c(colnames(X),"Y","eps")
    tree.gen <- mob(Y ~ X1+X2+X3+X4+X5+eps | X1+X2+X3+X4+X5,
                    data = data.gen, model  = linearModel)
    train.nodes = predict(tree.gen, as.data.frame(cbind(X,eps=noise.gen)), type='node')
    uni_nodes.vec[m] = length(unique(train.nodes))
    noise.test.gen = rnorm(n.test)
    test.nodes = predict(tree.gen, as.data.frame(cbind(X.test, eps=noise.test.gen)), type='node')
    out = helper(train.nodes, test.nodes, X, Y, X.test)
    ci.gen.low[,m] = out$ci.test[,1]
    ci.gen.up[,m] = out$ci.test[,2]
    ci.cond.gen.low[,m] = out$ci.cond.test[,1]
    ci.cond.gen.up[,m] = out$ci.cond.test[,2]
  }
  
  summ.mat = matrix(NA, nrow=n.test, ncol=2); colnames(summ.mat) = c('cov','len')
  summ.cond.mat = matrix(NA, nrow=n.test, ncol=2); colnames(summ.cond.mat) = c('cov','len')
  ci.i.mat = ci.i.cond.mat = matrix(NA, nrow=M, ncol=2)
  for(i in 1:n.test){
    for(m in 1:M){
      ci.i.mat[m,] = c(ci.gen.low[i,m], ci.gen.up[i,m])
      ci.i.cond.mat[m,] = c(ci.cond.gen.low[i,m], ci.cond.gen.up[i,m])
    }
    ci.i.summ = as.matrix(interval_union(Intervals(ci.i.mat)))
    ci.i.cond.summ = as.matrix(interval_union(Intervals(ci.i.cond.mat)))
    summ.mat[i,1] = sum((ci.i.summ[,2] > Y.test[i])*(ci.i.summ[,1] < Y.test[i]))
    summ.cond.mat[i,1] = sum((ci.i.cond.summ[,2] > Y_star.test[i])*(ci.i.cond.summ[,1] < Y_star.test[i]))
    summ.mat[i,2] = sum(ci.i.summ[,2] - ci.i.summ[,1])
    summ.cond.mat[i,2] = sum(ci.i.cond.summ[,2] - ci.i.cond.summ[,1])
  }
  cov.repro = len.repro = cov.cond.repro = len.cond.repro = rep(NA, L)
  for(l in 1:L){
    map.test.l = (group.test == l)
    cov.repro[l] = mean(summ.mat[map.test.l,1])
    len.repro[l] = mean(summ.mat[map.test.l,2])
    cov.cond.repro[l] = mean(summ.cond.mat[map.test.l,1])
    len.cond.repro[l] = mean(summ.cond.mat[map.test.l,2])
  }
  
  ####### Summary #######
  report.summary.pred = matrix(NA, nrow=L, ncol=2*4)
  colnames(report.summary.pred) = c('truth-C','truth-L','post-C','post-L','orac-C','orac-L','prop-C','prop-L')
  report.summary.cond = matrix(NA, nrow=L, ncol=2*4)
  colnames(report.summary.cond) = c('truth-C','truth-L','post-C','post-L','orac-C','orac-L','prop-C','prop-L')
  report.summary.pred[,1] = check(group.test, Y.test, out.truth$ci.test)$cov
  report.summary.pred[,2] = check(group.test, Y.test, out.truth$ci.test)$len
  report.summary.pred[,3] =  check(group.test, Y.test, out.post$ci.test)$cov
  report.summary.pred[,4] =  check(group.test, Y.test, out.post$ci.test)$len
  report.summary.pred[,5] =  check(group.test, Y.test, out0$ci.test)$cov
  report.summary.pred[,6] =  check(group.test, Y.test, out0$ci.test)$len
  report.summary.pred[,7] = cov.repro
  report.summary.pred[,8] = len.repro
  
  report.summary.cond[,1] = check.cond(group.test, Y_star.test, out.truth$ci.cond.test)$cov
  report.summary.cond[,2] = check.cond(group.test, Y_star.test, out.truth$ci.cond.test)$len
  report.summary.cond[,3] =  check(group.test, Y_star.test, out.post$ci.cond.test)$cov
  report.summary.cond[,4] =  check(group.test, Y_star.test, out.post$ci.cond.test)$len
  report.summary.cond[,5] =  check(group.test, Y_star.test, out0$ci.cond.test)$cov
  report.summary.cond[,6] =  check(group.test, Y_star.test, out0$ci.cond.test)$len
  report.summary.cond[,7] = cov.cond.repro
  report.summary.cond[,8] = len.cond.repro
  
  info.pred.list[[i.sim]] = report.summary.pred
  info.cond.list[[i.sim]] = report.summary.cond
  uni_nodes.list[[i.sim]] = uni_nodes.vec
}
filename = paste0('simu0119','-group',group.setting,'-coef',coef.setting,
                  '-noi',noi.str.setting,'-Sigma',Sigma.setting,'-M',M,
                  '-round',sim.round,'.RData')
save(info.pred.list, info.cond.list, uni_nodes.list, file=filename)