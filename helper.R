helper <- function(train.nodes, test.nodes, x.train, y.train, x.test){
  uni_nodes = sort(unique(train.nodes))
  n_nodes = length(uni_nodes)
  beta.list = C.list = rep(list(NA), n_nodes)
  pred.test = se.test = se.cond.test = rep(NA, length(test.nodes))
  sigmahat.vec = rep(NA,n_nodes); names(sigmahat.vec) = uni_nodes
  for(i.node in 1:n_nodes){
    node = uni_nodes[i.node]
    map = (train.nodes == node)
    sublm = lm(y.train[map]~x.train[map,,drop=F])
    beta = coef(sublm); names(beta)=NULL
    sigmasq.hat = sum(resid(sublm)^2)/(sum(map)-1-ncol(x.train))
    C = vcov(sublm); names(C)=NULL
    beta.list[[i.node]] = beta; C.list[[i.node]] = C; sigmahat.vec[i.node] = sigmasq.hat
    
    map.test = (test.nodes == node)
    if(sum(map.test)==0) next
    if(sum(map.test)==1){
      pred.test[map.test] = c(1, x.test[map.test,])%*%beta
      se.test[map.test] = sqrt(c(1,x.test[map.test,])%*%C%*%c(1,x.test[map.test,]) + sigmasq.hat)
      se.cond.test[map.test] = sqrt(c(1,x.test[map.test,])%*%C%*%c(1,x.test[map.test,]))
    }else{
      pred.test[map.test] = cbind(1, x.test[map.test,])%*%beta
      se.test[map.test] = sqrt(diag(cbind(1, x.test[map.test,])%*%C%*%t(cbind(1, x.test[map.test,])))
                               +sigmasq.hat)
      se.cond.test[map.test] = sqrt(diag(cbind(1, x.test[map.test,])%*%C%*%t(cbind(1, x.test[map.test,]))))
    }
  }
  ci.test = cbind(pred.test - qnorm(0.975)*se.test,
                  pred.test + qnorm(0.975)*se.test)
  ci.cond.test = cbind(pred.test - qnorm(0.975)*se.cond.test,
                       pred.test + qnorm(0.975)*se.cond.test)
  return(list(pred.test = pred.test,
              se.test = se.test,
              ci.test = ci.test,
              ci.cond.test = ci.cond.test,
              beta.list = beta.list,
              C.list = C.list,
              sigmas = sqrt(sigmahat.vec)))
}

check <- function(group.test, Y.test, ci.test){
  L = length(unique(group.test))
  cov = len = rep(NA, L)
  for(l in 1:L){
    map.test.l = (group.test == l)
    cov.l = mean((Y.test[map.test.l] < ci.test[map.test.l,2])*
                   (Y.test[map.test.l] > ci.test[map.test.l,1]))
    cov[l] = cov.l
    len[l] = mean(ci.test[map.test.l,2] - ci.test[map.test.l,1])
  }
  return(list(cov=cov,
              len=len))
}

check.cond <- function(group.test, cond.test, ci.cond.test){
  L = length(unique(group.test))
  cov = len = rep(NA, L)
  for(l in 1:L){
    map.test.l = (group.test == l)
    cov.l = mean((cond.test[map.test.l] < ci.cond.test[map.test.l,2])*
                   (cond.test[map.test.l] > ci.cond.test[map.test.l,1]))
    cov[l] = cov.l
    len[l] = mean(ci.cond.test[map.test.l,2] - ci.cond.test[map.test.l,1])
  }
  return(list(cov=cov,
              len=len))
}