# Partitioned Mean Models -----------------------------------------------
# Y_k = mu_k + sigma_k N(0,1) 
# Compare the following methods:
# 1. repro100: all K
# 2. repro95: remove 5\% outlier of K
# the original code is in simu0206 and simu0206-2
library(rpart)
library(MASS)
library(intervals)
library(randomForest)
library(party)
library(intervals)
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
}else if(group.setting == 4){
  get.group = function(x){
    if(x[1]<=-0.5){return(1)}
    if(x[1]>-0.5){
      if((x[2]<=0)&(x[3]<=0.5)){
        return(2)
      }else if((x[4]<=0)&(x[5]<=0.5)){
        return(3)
      }else{
        return(4)
      }
    }
  }
}

## Noise Level in each partition -------------------------------------------
noi.str = c(1.5, 1.2, 1, 0.8)

## Group Mean in each partition -------------------------------------------
group_mean = c(1, 1.2, 0.8, 1.5)

## Setup
L = 4
n = 1000
p = 5
if(Sigma.setting == 1){
  Sigma = diag(p)
}else if(Sigma.setting == 2){
  Sigma = toeplitz(0.8^{seq(0,4,length.out=p)})
}
M = ntree = 250


## Main Code ---------------------------------------------------------------
MSE.info = matrix(NA, nrow=nsim, ncol=3)
colnames(MSE.info) = c('rf','repro100','repro95')
for(i.sim in 1:nsim){
  print(paste0('=====>i.sim=',i.sim,'/',nsim))
  X = mvrnorm(n, mu=rep(0,p), Sigma)
  colnames(X) = paste0('X',1:p,seq="")
  group = apply(X, 1, get.group)
  Xint = cbind(1, X)
  Y_star = group_mean[group]
  noise = rep(NA, n)
  for(i in 1:n){
    noise[i] = rnorm(1, sd=noi.str[group[i]])
  }
  Y = Y_star + noise
  
  set.seed(2021)
  n.test=1000
  X.test=mvrnorm(n.test, mu=rep(0,p), Sigma)
  colnames(X.test)=paste0('X',1:p,seq="")
  group.test = apply(X.test,1 ,get.group)
  Y_star.test = group_mean[group.test]
  
  rf = randomForest(X, Y, ntree=ntree, keep.forest = TRUE)
  out.rf = predict(rf, newdata=X.test, predict.all = TRUE)
  
  # M = 500
  pred = matrix(NA, nrow=n.test, ncol=M)
  uni_nodes.vec = rep(NA, M)
  for(m in 1:M){
    print(paste0('repro-',m,'/',M))
    noise.gen = rnorm(n)
    data.gen = as.data.frame(cbind(X, Y, noise.gen))
    colnames(data.gen) = c(colnames(X),"Y","eps")
    tree.gen <- mob(Y~eps | X1+X2+X3+X4+X5,
                    data = data.gen, model  = linearModel)
    train.nodes = predict(tree.gen, as.data.frame(cbind(X,eps=noise.gen)), type='node')
    uni_nodes.vec[m] = length(unique(train.nodes))
    noise.test.gen = rnorm(n.test)
    test.nodes = predict(tree.gen, as.data.frame(cbind(X.test, eps=noise.test.gen)), type='node')
    out = helper(train.nodes, test.nodes, X, Y, X.test)
    pred[,m] = out$pred.test
  }
  prob.nodes = sort(prop.table(table(uni_nodes.vec)), decreasing=T)
  alpha = 0.05
  nodes.filter = as.numeric(names(prob.nodes[1:(which(cumsum(prob.nodes) >= 0.999-alpha)[1])]))
  M.filter = (1:M)[uni_nodes.vec %in% nodes.filter]
  pred.filter = pred[,M.filter]
  
  # MSE for every data point
  error.mat = matrix(NA, nrow=n.test, ncol=3); colnames(error.mat) = c('rf','repro100','repro95')
  for(i.test in 1:n.test){
    error.mat[i.test,1] = mean(out.rf$individual[i.test,]) - Y_star.test[i.test]
    error.mat[i.test,2] = mean(pred[i.test,]) - Y_star.test[i.test]
    error.mat[i.test,3] = mean(pred.filter[i.test,]) - Y_star.test[i.test]
  }
  MSE.info[i.sim,] = apply(error.mat, MARGIN=2, FUN=function(x) mean(x^2))
}
filename = paste0('simu0206','-group',group.setting,
                  '-noi',noi.str.setting,'-Sigma',Sigma.setting,
                  '-round',sim.round,'.RData')
save(MSE.info, file=filename)

