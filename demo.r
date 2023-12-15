
## A toy example to show how to use 'RGCR' package.

rm(list=ls());gc()
library(RGCR)
seed.rand = 314
set.seed(seed.rand)
W = getAdjMat('demoDATA//SW32')
n = ncol(W)
W.deg = W %*% rep(1,n)
W.squared = getAdjMat.d(W,2)
eig.value =  eigen(W.squared)$vectors[,1]
delta = seq(0,1,0.05)
len.delta = length(delta)
p = 0.5
nruns = 200

y0.ugan = getY0.ugander('demoDATA//SW32')
ugan.beta = 0.5
ugan.gamma = 2.5

ep1h.3net.com.uniform = read.csv('demoEPest//ep1h-3net-com-uniform.csv')
ep0h.3net.com.uniform = read.csv('demoEPest//ep0h-3net-com-uniform.csv')
ep1h.1hop.com.uniform = read.csv('demoEPest//ep1h-1hop-com-uniform.csv')
ep0h.1hop.com.uniform = read.csv('demoEPest//ep0h-1hop-com-uniform.csv')

data.names = c('Ugander\'s (misspecified)')
design.names = c('3net-com-uniform','1hop-com-uniform')
estimator.names = c('ht')

res.3net = data.frame(matrix(0,nruns,length(delta)))
res.1hop = data.frame(matrix(0,nruns,length(delta)))

for(nrun in 1:nruns) {
  cat('nrun=',nrun,'...\n')
  set.seed((nrun-1)*len.delta)
  
  ## 3net-com-uniform-ht
  z = RGCR(W, 
           cluster.type = '3net', 
           rand.type='com', 
           weight.type='uniform',
           p=p,
           prior=NA,
           W.squared=W.squared,
           eig.value=eig.value)$treat
  valid = getValidNode.delta(W,z,delta)
  valid.1 = valid$valid.1
  valid.0 = valid$valid.0
  ratio1 = as.matrix(valid.1/ep1h.3net.com.uniform)
  ratio0 = as.matrix(valid.0/ep0h.3net.com.uniform)
  
  y = getY.Ugander(W,z,y0.ugan,beta=ugan.beta,gamma=ugan.gamma)
  res.3net[nrun,] = (t(y) %*% ratio1)/n - rev((t(y) %*% ratio0)/n)
  
  ## 1hop-com-uniform-ht
  z = RGCR(W, 
           cluster.type = '1hop', 
           rand.type='com', 
           weight.type='uniform',
           p=p,
           prior=NA,
           W.squared=W.squared,
           eig.value=eig.value)$treat
  valid = getValidNode.delta(W,z,delta)
  valid.1 = valid$valid.1
  valid.0 = valid$valid.0
  ratio1 = as.matrix(valid.1/ep1h.1hop.com.uniform)
  ratio0 = as.matrix(valid.0/ep0h.1hop.com.uniform)
  
  y = getY.Ugander(W,z,y0.ugan,beta=ugan.beta,gamma=ugan.gamma)
  res.1hop[nrun,] = (t(y) %*% ratio1)/n - rev((t(y) %*% ratio0)/n)
  
}

tau.ugan = mean(getY.Ugander(W,z=rep(1,n),y0.ugan,beta=ugan.beta,gamma=ugan.gamma)-y0.ugan)
mse.ugan = c()
mse.ugan = rbind(mse.ugan,
                 apply((res.3net-tau.ugan)^2,2,mean),
                 apply((res.1hop-tau.ugan)^2,2,mean))
var.ugan = c()
var.ugan = rbind(var.ugan,
                 apply((res.3net-tau.ugan),2,var),
                 apply((res.1hop-tau.ugan),2,var))
bias2.ugan = mse.ugan - var.ugan

par(mfrow=c(1,3))
plot(delta,mse.ugan[1,],ylim = c(0,12),xlab='Delta for working models', ylab='mse',main='MSE: 3-net v.s. 1-hop',col='red',pch=15,lty=1)
lines(delta,mse.ugan[1,],col='red',type='l',pch=15,lty=1)
lines(delta,mse.ugan[2,],col='blue',type='p',pch=16,lty=2)
lines(delta,mse.ugan[2,],col='blue',type='l',pch=16,lty=2)
legend("topright", inset=.05, c('3net','1hop'),pch=c(15,16),lty=c(1,2), col=c("red", "blue"))
plot(delta,var.ugan[1,],ylim = c(0,12),xlab='Delta for working models', ylab='var',main='Variance: 3-net v.s. 1-hop',col='red',pch=15,lty=1)
lines(delta,var.ugan[1,],col='red',type='l',pch=15,lty=1)
lines(delta,var.ugan[2,],col='blue',type='p',pch=16,lty=2)
lines(delta,var.ugan[2,],col='blue',type='l',pch=16,lty=2)
legend("topright", inset=.05, c('3net','1hop'),pch=c(15,16),lty=c(1,2), col=c("red", "blue"))
plot(delta,bias2.ugan[1,],ylim = c(0,12),xlab='Delta for working models', ylab='bias2',main='Bias2: 3-net v.s. 1-hop',col='red',pch=15,lty=1)
lines(delta,bias2.ugan[1,],col='red',type='l',pch=15,lty=1)
lines(delta,bias2.ugan[2,],col='blue',type='p',pch=16,lty=2)
lines(delta,bias2.ugan[2,],col='blue',type='l',pch=16,lty=2)
legend("topright", inset=.05, c('3net','1hop'),pch=c(15,16),lty=c(1,2), col=c("red", "blue"))





