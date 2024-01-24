set.seed(0)
n0 = 1000
p = 6


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

b <- c(c(-2,-1,0,1,2)/5, 0)
e <- function(X, b){
  val = exp(X%*%b)/(1+ exp(X%*%b))
  return(as.vector(val))
}

if(Xsetting==1){
  X0 = matrix(runif(n0,min=-1,max=1), nrow=n0, ncol=p)
  filename = 'rf1128-testdata-unif.rds'
}
if(Xsetting==2){
  X0 = MASS::mvrnorm(n0, mu=rep(0,p), Sigma = diag(p))
  filename = 'rf1128-testdata-gaus.rds'
}
colnames(X0) = paste0('X',1:p)
eX0 = e(X0,b)
W0 = rbinom(n0, size=1, eX0)
Y0 = (1/2)*(2*W0-1)*kappa(X0) + rnorm(n0, sd=0.1)
Y0.pure = (1/2)*(2*W0-1)*kappa(X0)
Y0_prox = W0*Y0/eX0 - (1-W0)*Y0/(1-eX0)
Y0_prox.pure = W0*Y0.pure/eX0 - (1-W0)*Y0.pure/(1-eX0)


saveRDS(list(n0 = n0,
             X0 = X0,
             W0 = W0,
             Y0 = Y0,
             Y0.pure = Y0.pure,
             Y0_prox = Y0_prox,
             Y0_prox.pure = Y0_prox.pure), filename)