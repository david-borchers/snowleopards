
library(secr)
library(fields)
library(maptools)
library(gdistance)
library(igraph)
library(fields)
library(raster)

# get the fits
load("./Tost/Tost-nonEuc-fits2.RData")
# load plotting functions
source("lcplots.r")

# Plot riggedness, with numbered traps overlaid
# --------------------------
plotcovariate(TostSurface.D.nonU,covariate="stdGC",asp=1,contour=FALSE,col=terrain.colors(40))
text(Tost.cams,labels=as.character(1:40),cex=0.75)

# Predict and plot log(density estiamte)
# --------------------------
TostSurface.nonU<-predictDsurface(Tost.hhn.DHab.nonU, se.D=TRUE, cl.D=TRUE)
covariates(TostSurface.nonU)$logD = log(covariates(TostSurface.nonU)$D.0)
pdat=plotcovariate(TostSurface.nonU,covariate="logD",asp=1,contour=FALSE)


# Plot some least-cost paths
# --------------------------
# Get parameter estimates and make noneuc covariate and put on mask:
lambda0=exp(coef(Tost.hhn.DHab.nonU)["lambda",1]) # on the real scale
sigma=exp(coef(Tost.hhn.DHab.nonU)["sigma",1]) # on the real scale
alpha=coef(Tost.hhn.DHab.nonU)["noneuc.stdGC",1] # on the beta scale
covariates(TostMask1)$noneuc=exp(-alpha*covariates(TostMask1)$stdGC) # add noneuc covariate
# Specify points from and to which to plot paths (the ones here are all traps)
from=matrix(
  c(
    Tost.cams[7,1],Tost.cams[7,2],
    Tost.cams[12,1],Tost.cams[12,2],
    Tost.cams[40,1],Tost.cams[40,2],
    Tost.cams[40,1],Tost.cams[40,2],
    Tost.cams[40,1],Tost.cams[40,2],
    Tost.cams[40,1],Tost.cams[40,2],
    Tost.cams[40,1],Tost.cams[40,2],
    Tost.cams[12,1],Tost.cams[12,2]
    ), byrow=TRUE,ncol=2
  )
to=matrix(
  c(
    Tost.cams[29,1],Tost.cams[29,2],
    Tost.cams[29,1],Tost.cams[29,2],
    Tost.cams[29,1],Tost.cams[29,2],
    Tost.cams[33,1],Tost.cams[33,2],
    Tost.cams[39,1],Tost.cams[39,2],
    Tost.cams[24,1],Tost.cams[24,2],
    Tost.cams[20,1],Tost.cams[20,2],
    Tost.cams[36,1],Tost.cams[36,2]
  ), byrow=TRUE,ncol=2
)
# define cost function
cfun=function(x) mean(x) 
# Plot least-cost paths
plot_lcpath(from,to,TostMask1,costfun=cfun,plotcovariate="stdGC",lwd=2,linecol="white",col=terrain.colors(40))




# Some plots of prob of being detected by given traps:
par(mfrow=c(2,2))
j=7 
predicted.pj=pdot(TostMask1,Tost.cams[j,],noccasions=1,detectfn='HHN',
                  detectpar=list(lambda0=lambda0,sigma=sigma),userdist=userdfn1)
covariates(TostMask1)$predicted.pj=predicted.pj
plotcovariate(TostMask1,covariate="predicted.pj",contour=FALSE,asp=1)
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,col="white")) 
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,cex=0.25,col="black"))
j=32 
predicted.pj=pdot(TostMask1,Tost.cams[j,],noccasions=1,detectfn='HHN',
                  detectpar=list(lambda0=lambda0,sigma=sigma),userdist=userdfn1)
covariates(TostMask1)$predicted.pj=predicted.pj
plotcovariate(TostMask1,covariate="predicted.pj",contour=FALSE,asp=1)
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,col="white")) 
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,cex=0.25,col="black"))
j=36
predicted.pj=pdot(TostMask1,Tost.cams[j,],noccasions=1,detectfn='HHN',
                  detectpar=list(lambda0=lambda0,sigma=sigma),userdist=userdfn1)
covariates(TostMask1)$predicted.pj=predicted.pj
plotcovariate(TostMask1,covariate="predicted.pj",contour=FALSE,asp=1)
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,col="white")) 
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,cex=0.25,col="black"))
j=39 
predicted.pj=pdot(TostMask1,Tost.cams[j,],noccasions=1,detectfn='HHN',
                  detectpar=list(lambda0=lambda0,sigma=sigma),userdist=userdfn1)
covariates(TostMask1)$predicted.pj=predicted.pj
plotcovariate(TostMask1,covariate="predicted.pj",contour=FALSE,asp=1)
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,col="white")) 
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,cex=0.25,col="black"))

