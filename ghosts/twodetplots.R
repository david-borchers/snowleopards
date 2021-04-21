library(secr)
library(maptools)
library(sp) # package to manipulate and plot spatial objects

#----------------------- Load the RData objects made by Make_TNN_RData.r --------------------------
load("./Analysis4paper/TNN_boundaries.RData") # Tostboundary,Noyonboundary,Nemegtboundary
load("./Analysis4paper/TNN_masks.RData") # TostMask,NoyonMask,NemegtMask
load("./Analysis4paper/TNN_caphists.RData") # Tost_ch,Noyon_ch,Nemegt_sh,TNN_ch)
#----------------------- ----------------------------------------------- --------------------------
summary(Tostboundary)
summary(TostMask)
summary(Tost_ch)
Tost_traps = traps(Tost_ch)
plot(TostMask)
plot(Tostboundary,add=TRUE)
plot(Tost_traps,add=TRUE)

# Fit a simple model
Tostfit.null <- secr.fit(Tost_ch, detectfn="HHN", mask=TostMask,model=list(D~1, lambda0~1, sigma~1))

detectpars = detectpar(Tostfit.null)
l0 = detectpars$lambda0
sig = detectpars$sigma
distmat =  edist(TostMask, Tost_traps)
coords = as.data.frame(TostMask)

p. = pdot(TostMask, Tost_traps, detectfn="HHN", detectpar=list(lambda0=l0, sigma=sig), noccasions=1) 
tm = usage(Tost_traps)
mypdot = function(d,lambda0,sigma,tm) return(1 - exp(-sum(tm*lambda0*exp(-d^2/sigma^2))))
p.once.HHN = function(d,lambda0,sigma,tm) {
  lambda = sum(tm*lambda0*exp(-d^2/sigma^2))
  return(lambda*exp(-lambda))
}
En.HHN = function(d,lambda0,sigma,tm) return(sum(tm*lambda0*exp(-d^2/sigma^2)))
p1 = p.once.HHN(distmat[1,],l0,sig,tm)
En = En.HHN(distmat[1,],l0,sig,tm)
myp. = apply(distmat,1,mypdot,lambda0=l0,sigma=sig,tm=tm)
p1s = apply(distmat,1,p.once.HHN,lambda0=l0,sigma=sig,tm=tm)
Ens = apply(distmat,1,En.HHN,lambda0=l0,sigma=sig,tm=tm)

p2. = p. - p1s
coverage. = sum(p.)/length(p.)
coverage2. = sum(p2.)/length(p.)
pcmore = 100*(coverage./coverage2.-1)
round(pcmore)

spts = SpatialPoints(coords=coords)
sp. = SpatialPixelsDataFrame(spts, data=data.frame(p.=p.))
sp2. = SpatialPixelsDataFrame(spts, data=data.frame(p2.=p2.))
bias = 100*(1-coverage2./coverage.)
round(bias,1)

par(mfrow=c(2,1))
plot(sp.,main=paste("p(det); coverage=",round(100*coverage.),"%",sep=""),what="image")
plot(Tost_traps,add=TRUE)
plot(sp2.,main=paste("p(det at least twice); coverage=",round(100*coverage2.),"%",sep=""),what="image")
plot(Tost_traps,add=TRUE)

#plot(sp.,main=paste("p(det); coverage=",round(100*coverage.),"%",sep=""),what="scale")


dsp. = SpatialPixelsDataFrame(spts, data=data.frame(dp.=p.-p2.))
plot(dsp.,main="Heatmap of prob(detect>1)-prob(detect)")
sp1. = SpatialPixelsDataFrame(spts, data=data.frame(p1s=p1s))
plot(sp1.,main="Heatmap of prob(detect once)")


