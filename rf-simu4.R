# Generic Function Fast/Slow Comparison ------------------------------------------
# This code is to compare random forests and reproRF-fast, reproRF-slow
# reproRF-fast: use rpart package to fit each tree
# reproRF-slow: use mob to fit each tree
# the original code is located in "simu0223"
library(MASS)
library(rpart)
library(randomForest)
library(party)

sim.round = 1 #as.numeric(readLines(f, n=1)) # {1,2,3,4,5}
setting = 1 #as.numeric(readLines(f, n=1)) # {1,2,3}
p = 5 #as.numeric(readLines(f, n=1)) # {5,10,20}
n = 200 #as.numeric(readLines(f, n=1)) # {100,200,500}

## specify functions f(x) --------------------------------------------------
if(setting == 1){
  f <- function(X){
    return(10*sin(10*pi*X[,1]))
  }
}
if(setting == 2){
  f <- function(X){
    return(10*sin(pi*X[,1]*X[,2])+20*(X[,3] - 0.05)^2+10*X[,4]+5*X[,5])
  }
}
if(setting == 3){
  f.x <- function(x){
    if(x[4]<0.383){
      if(x[2]<0.2342){
        val = 8.177
      }else{
        if(x[1]<0.2463){
          val = 8.837
        }else{
          val = 13.15
        }
      }
    }else{
      if(x[1]<0.47){
        if(x[5]<0.2452){
          val = 10.99
        }else{
          if(x[3]>=0.2234){
            val = 13.87
          }else{
            val = 18.03
          }
        }
      }else{
        if(x[2]<0.2701){
          val = 15.02
        }else{
          if(x[5]<0.5985){
            val = 18.61
          }else{
            val = 21.74
          }
        }
      }
    }
    return(val)
  }
  f <- function(X){
    apply(X, MARGIN=1, FUN=f.x)
  }
}


nsim = 50
MSE.mat = matrix(NA, nrow=nsim, ncol=3); colnames(MSE.mat) = c('rf','fast','slow')
time.cost.mat = matrix(NA, nrow=nsim, ncol=3); colnames(time.cost.mat) = c('rf','fast','slow')
for(i.sim in 1:nsim){
  print(paste0('isim-',i.sim,'/',nsim))
  X = matrix(runif(n*p), nrow=n, ncol=p)#mvrnorm(n, rep(0,p), diag(p))
  Y = f(X)+rnorm(n)
  
  n0 = 50000
  X0 = matrix(runif(n0*p), nrow=n0, ncol=p)#mvrnorm(n0,rep(0,p), diag(p))
  Y0.star = f(X0)
  Y0 = Y0.star+rnorm(n0)
  
  colnames(X) = paste0('X',1:p)
  colnames(X0) = colnames(X)
  
  ## random forests ##
  start.time = Sys.time()
  ntree = M = 200
  rf = randomForest(X, Y, ntree=ntree, keep.forest = TRUE)
  pred.rf = predict(rf, newdata=X0, predict.all = TRUE)
  end.time = Sys.time()
  cost.rf = difftime(end.time, start.time, units='secs')
  
  ## repro ##
  pred.fast = pred.slow = matrix(NA, nrow=n0, ncol=M)
  time.cost = c(0,0)
  for(m in 1:M){
    print(paste0('=======>repro-',m,'/',M))
    noise.gen = rnorm(n)
    
    ## repro fast ##
    start.time = Sys.time()
    Y.gen = Y - noise.gen
    data.gen = data.frame(cbind(X, Y.gen))
    colnames(data.gen) = c(colnames(X), 'Y')
    fmla = as.formula(paste('Y ~ ',paste(colnames(X), collapse = '+')))
    fit.fast.m = rpart(fmla, data=data.gen)
    out.fast.m = predict(fit.fast.m, data.frame(X0))
    pred.fast[,m] = out.fast.m
    end.time = Sys.time()
    time.cost[1] = time.cost[1] + difftime(end.time, start.time, units='secs')
    
    ## repro slow ##
    start.time = Sys.time()
    data.gen = as.data.frame(cbind(X, Y, noise.gen))
    colnames(data.gen) = c(colnames(X),"Y","eps")
    fmla = as.formula(paste0('Y~', paste(colnames(X), collapse = '+'), '+eps', '|', paste(colnames(X), collapse = '+')))
    tree.gen <- mob(fmla,data = data.gen, model = linearModel)
    noise0.gen = rnorm(n0)
    fit.slow.m = predict(tree.gen, newdata=data.frame(cbind(X0, eps=noise0.gen)))
    pred.slow[,m] = fit.slow.m
    end.time = Sys.time()
    time.cost[2] = time.cost[2] + difftime(end.time, start.time, units='secs')
  }
  
  error.mat = matrix(NA, nrow=n0, ncol=3); colnames(error.mat) = c('rf','fast','slow')
  for(i0 in 1:n0){
    error.mat[i0,1] = mean(pred.rf$individual[i0,]) - Y0.star[i0]
    error.mat[i0,2] = mean(pred.fast[i0,]) - Y0.star[i0]
    error.mat[i0,3] = mean(pred.slow[i0,]) - Y0.star[i0]
  }
  MSE.mat[i.sim,] = apply(error.mat, MARGIN=2, FUN=function(x) mean(x^2))
  time.cost.mat[i.sim,] = c(cost.rf, time.cost)
}
filename = paste0('simu0223','-setting',setting,'-p',p,'-n',n,
                  '-round',sim.round,'.RData')
save(MSE.mat, time.cost.mat, file=filename)

