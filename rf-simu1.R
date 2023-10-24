# Generic Function Investigation ------------------------------------------
# This code is to compare random forests and repro RF for generic functions
# the original code with RData file is named "simu_0227"
library(rpart)
library(randomForest)
library(MASS)

setting = 4 # {1,2,3,4,5,6}
p = 5
n = 1000 # sample size for source data
n0 = 5000 # sample size for target data

## specify functions f(x) --------------------------------------------------
if(setting == 1){
  f <- function(X){
    return(1*sin(10*pi*X[,1]))
  }
}
if(setting == 2){
  f <- function(X){
    return(1*sin(pi*X[,1]*X[,2])+2*(X[,3] - 0.05)^2+1*X[,4]+0.5*X[,5])
  }
}
if(setting == 3){
  f <- function(X){
    return(sin(pi*X[,1])+sin(pi*X[,2])+(X[,3]-0.5)^2+X[,4]+0.5*X[,5])
  }
}
if(setting == 4){
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
if(setting == 5){
  f <- function(X){
    return(10*sin(10*pi*X[,1]))
  }
}
if(setting == 6){
  f <- function(X){
    return(10*sin(pi*X[,1]*X[,2])+20*(X[,3] - 0.05)^2+10*X[,4]+5*X[,5])
  }
}


# Main Code ---------------------------------------------------------------
nsim = 5
MSE.mat = matrix(NA, nrow=nsim, ncol=3)
for(i.sim in 1:nsim){
  if(i.sim %% 10 == 1){
    print(paste0('isim-',i.sim,'/',nsim))
  }
  X = mvrnorm(n, rep(0,p), diag(p))
  X0 = mvrnorm(n0,rep(0,p), diag(p))
  Y = f(X)+rnorm(n)
  Y0.star = f(X0) # without noise target label
  Y0 = Y0.star+rnorm(n0) 
  
  colnames(X) = paste0('X',1:p)
  colnames(X0) = colnames(X)
  
  ## random forests ##
  start.time = Sys.time()
  ntree = M = 500
  rf = randomForest(X, Y, ntree=ntree, keep.forest = TRUE)
  pred.rf = predict(rf, newdata=X0, predict.all = TRUE)
  end.time = Sys.time()
  cost.rf = difftime(end.time, start.time, units='secs')
  
  ## repro ##
  sigma.hat = sqrt(mean((Y-pred.rf$aggregate)^2))
  pred.oracle = pred.hat = matrix(NA, nrow=n0, ncol=M)
  for(m in 1:M){
    if(m %% 20 == 1){
      print(paste0('--------->repro-',m,'/',M))
    }
    noise0.gen = rnorm(n)
    noise1.gen = rnorm(n, sd=sigma.hat)
    ## repro fast oracle ##
    Y.gen = Y - noise0.gen
    data.gen = data.frame(cbind(X, Y.gen))
    colnames(data.gen) = c(colnames(X), 'Y')
    fmla = as.formula(paste('Y ~ ',paste(colnames(X), collapse = '+')))
    fit.m = rpart(fmla, data=data.gen)
    out.m = predict(fit.m, data.frame(X0))
    pred.oracle[,m] = out.m
    ## repro fast hat ##
    Y.gen = Y - noise1.gen
    data.gen = data.frame(cbind(X, Y.gen))
    colnames(data.gen) = c(colnames(X), 'Y')
    fmla = as.formula(paste('Y ~ ',paste(colnames(X), collapse = '+')))
    fit.m = rpart(fmla, data=data.gen)
    out.m = predict(fit.m, data.frame(X0))
    pred.hat[,m] = out.m
  }
  
  error.mat = matrix(NA, nrow=n0, ncol=3)
  for(i0 in 1:n0){
    error.mat[i0,1] = mean(pred.rf$individual[i0,]) - Y0.star[i0]
    error.mat[i0,2] = mean(pred.oracle[i0,]) - Y0.star[i0]
    error.mat[i0,3] = mean(pred.hat[i0,]) - Y0.star[i0]
  }
  MSE.mat[i.sim,] = colMeans((error.mat)^2)
}


## Output ------------------------------------------------------------------
# The output is MSE measured compared with Y0.star
# 1st column: random forest
# 2nd column: repro with known sigma^2
# 3rd column: repro with estimated sigma^2
MSE.mat
