

## This file only needs to be run once. This may take some time (10min or so). 
## The results are stored in the 'demoEPest' folder. 
## I have run it, so you can just go to demo.r.

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
p = 0.5
K = 2
L = 10

# 3net-com-uniform
cluster.type = '3net'
rand.type = 'com'
weight.type = 'uniform'
simu.type = paste0(cluster.type,'-',rand.type,'-',weight.type)
set.seed(seed.rand)
EPh = estExposure.MC(W,
                     cluster.type=cluster.type,
                     rand.type=rand.type,
                     weight.type=weight.type,
                     p=p, K=K, L=L,
                     delta = delta,
                     W.squared=W.squared,
                     eig.value=eig.value)
write.csv(EPh$ep1h,file=paste0('demoEPest//ep1h-',simu.type,'.csv'),row.names = F)
write.csv(EPh$ep0h,file=paste0('demoEPest//ep0h-',simu.type,'.csv'),row.names = F)

# 3net-ind-uniform
cluster.type = '3net'
rand.type = 'ind'
weight.type = 'uniform'
simu.type = paste0(cluster.type,'-',rand.type,'-',weight.type)
set.seed(seed.rand)
EPh = estExposure.MC(W,
                     cluster.type=cluster.type,
                     rand.type=rand.type,
                     weight.type=weight.type,
                     p=p, K=K, L=L,
                     delta = delta,
                     W.squared=W.squared,
                     eig.value=eig.value)
write.csv(EPh$ep1h,file=paste0('demoEPest//ep1h-',simu.type,'.csv'),row.names = F)
write.csv(EPh$ep0h,file=paste0('demoEPest//ep0h-',simu.type,'.csv'),row.names = F)

# 1hop-com-uniform
cluster.type = '1hop'
rand.type = 'com'
weight.type = 'uniform'
simu.type = paste0(cluster.type,'-',rand.type,'-',weight.type)
set.seed(seed.rand)
EPh = estExposure.MC(W,
                     cluster.type=cluster.type,
                     rand.type=rand.type,
                     weight.type=weight.type,
                     p=p, K=K, L=L,
                     delta = delta,
                     W.squared=W.squared,
                     eig.value=eig.value)
write.csv(EPh$ep1h,file=paste0('demoEPest//ep1h-',simu.type,'.csv'),row.names = F)
write.csv(EPh$ep0h,file=paste0('demoEPest//ep0h-',simu.type,'.csv'),row.names = F)

# 1hop-ind-uniform
cluster.type = '1hop'
rand.type = 'ind'
weight.type = 'uniform'
simu.type = paste0(cluster.type,'-',rand.type,'-',weight.type)
set.seed(seed.rand)
EPh = estExposure.MC(W,
                     cluster.type=cluster.type,
                     rand.type=rand.type,
                     weight.type=weight.type,
                     p=p, K=K, L=L,
                     delta = delta,
                     W.squared=W.squared,
                     eig.value=eig.value)
write.csv(EPh$ep1h,file=paste0('demoEPest//ep1h-',simu.type,'.csv'),row.names = F)
write.csv(EPh$ep0h,file=paste0('demoEPest//ep0h-',simu.type,'.csv'),row.names = F)
