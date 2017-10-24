noneuc <- raster(nrows=3, ncols=2) 
noneuc <- setValues(noneuc,rep(1,6)) # all "noneuc" equal to 1
plot(noneuc)

#create a Transition object from the raster
tr <- transition(r,function(x) exp(mean(x)),8)
tr2 <- transition(r,function(x) exp(1/mean(x)),8)

#create two sets of coordinates
sP1 <- cbind(c(-100,-100,100,100),c(50,-50,50,-50))

# dists
LCd = costDistance(tr,sP1)
LCdmat = as.matrix(LCd);LCdmat
LCd2 = costDistance(tr2,sP1)
LCdmat2 = as.matrix(LCd2);LCdmat2

noneuca <- raster(nrows=3, ncols=2) 
noneuca <- setValues(noneuc,rep(0,6)) # all "noneuc" equal to 0 
plot(noneuca)
#create a Transition object from the raster
tra <- transition(noneuca,function(x) exp(mean(x)),8)
tr2a <- transition(noneuca,function(x) exp(1/mean(x)),8)
# dists
LCda = costDistance(tra,sP1)
LCdmata = as.matrix(LCda);LCdmata
LCd2a = costDistance(tr2a,sP1)
LCdmat2a = as.matrix(LCd2a);LCdmat2a

noneucs = c(unique(noneuc), unique(noneuca))
LCs = c(LCdmat[1,2], LCdmata[1,2])
plot(noneucs,LCs,xlab="noneuc",ylab="LC distance")
# Conclusion: with exp(mean(x)), "noneuc" is proportional to conductivity! (High noneuc => low LC distance)

noneucs = c(unique(noneuc), unique(noneuca))
LCs = c(LCdmat2[1,2], LCdmat2a[1,2])
plot(noneucs,LCs,xlab="noneuc",ylab="LC distance")
# Conclusion: with exp(1/mean(x)), "noneuc" is proportional to resistance (High noneuc => high LC distance)
