library(secr)
library(fields)
#Running SECR for Tost 2012
all.data.Tost<-read.capthist(captfile = "./Tost/Tost_capthist2012.csv", trapfile = "./Tost/Tost_cams_rugged2012.csv", detector="count", fmt = "trapID", trapcovnames = c("Effort",	"Topo",	"Altidute",	"Rgd"))
boundaryTost=readShapeSpatial("./Tost//Habitat/TostStudy_Area.shp")
plot(boundaryTost)
plot(x=all.data.Tost, add=TRUE)
TostMask=make.mask(traps(all.data.Tost), spacing=500, buffer = 25000, type="trapbuffer", poly=boundaryTost)
SLCost.Tost<-readShapePoly("./Tost//Habitat/Tost_Rgd500m.shp")  #ruggedness pixels averaged over 500m radius

head(SLCost.Tost)
head(TostMask)

TostMask1<-addCovariates(TostMask, SLCost.Tost)

# Standarize Rgd (this makes fits a bit more stable)
# --------------------------------------------------
summary(covariates(traps(all.data.Tost)))
covariates(traps(all.data.Tost))$stdRgd = scale(covariates(traps(all.data.Tost))$Rgd)
summary(covariates(traps(all.data.Tost)))
head(covariates(traps(all.data.Tost)))

# Also standarize stdGRIDCODE for good measure
# -----------------------------------------
summary(covariates(TostMask1))
covariates(TostMask1)$stdGC = scale(covariates(TostMask1)$GRIDCODE)
names(covariates(TostMask1))

head(covariates(TostMask1))
summary(covariates(TostMask1))
summary(covariates(traps(all.data.Tost)))

plot(TostMask1, covariate="stdGC", contour=FALSE, col=terrain.colors(10), legend = FALSE)
head(covariates(traps(all.data.Tost)))
head(covariates(TostMask1))


Tost.cams=traps(all.data.Tost)

#Tost.hhn<-secr.fit(all.data.Tost, model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", mask=TostMask1)
#Tost.hhn.detrug<-secr.fit(all.data.Tost, model=list(D~1, lambda0~Rgd, sigma~Rgd), detectfn="HHN", mask=TostMask1)
#Tost.hhn.Dxy<-secr.fit(all.data.Tost, model=list(D~x+y, lambda0~1, sigma~1), detectfn="HHN", mask=TostM)
Tost.hhn.DHab<-secr.fit(all.data.Tost, model=list(D~stdGC, lambda0~1, sigma~1), detectfn="HHN", mask=TostMask1)
#Tost.hhn.Dx<-secr.fit(all.data.Tost, model=list(D~x, lambda0~1, sigma~1), detectfn="HHN", mask=TostMask1)

#AIC(Tost.hhn.Dx, Tost.hhn.DHab, Tost.hhn.Dxy)
coefficients(Tost.hhn.DHab)
TostSurface<-predictDsurface(Tost.hhn.DHab, se.D=TRUE, cl.D=TRUE)
plot(TostSurface,asp=1,contour=FALSE)
plotcovariate(TostSurface,covariate="stdGC",asp=1,contour=FALSE)
plot(Tost.cams,add=TRUE)

Nhat1<-region.N(Tost.hhn.DHab) #Estimates the population N of the animals within the region defined by mask
Nhat1



# Non-Euclidian fits
# ==================
# This taken straight from secr vignette:
userdfn1 <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc')
  require(gdistance)
  Sraster <- raster(mask, 'noneuc')
  ## conductance is inverse of friction
  trans <- transition(Sraster, transitionFunction = function(x) 1/mean(x),  directions = 16)
  trans <- geoCorrection(trans)
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}
# variant not using inverse of mean
userdfn2 <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc')
  require(gdistance)
  Sraster <- raster(mask, 'noneuc')
  ## conductance is inverse of friction
  trans <- transition(Sraster, transitionFunction = function(x) mean(x),  directions = 16)
  trans <- geoCorrection(trans)
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}
source("scrplotting.r")

Tost.hhn.DHab.nonU<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                             model=list(D~stdGC, lambda0~1, sigma~1, 
                                        noneuc ~ stdGC -1), 
                             details = list(userdist = userdfn1),
                             start = list(noneuc = 1))
coefficients(Tost.hhn.DHab.nonU)
TostSurface.nonU<-predictDsurface(Tost.hhn.DHab.nonU, se.D=TRUE, cl.D=TRUE)
plot(TostSurface.nonU,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Tost.cams,add=TRUE)
plotcovariate(TostSurface.nonU,covariate="stdGC",asp=1,contour=FALSE)
plot(Tost.cams,add=TRUE)

Nhat1.nonU<-region.N(Tost.hhn.DHab.nonU)
Nhat1.nonU

# Try with flat density and non-Euclidian:
Tost.hhn.D.nonU<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                             model=list(D~1, lambda0~1, sigma~1, 
                                        noneuc ~ stdGC -1), 
                             details = list(userdist = userdfn1),
                             start = list(noneuc = 1))

coefficients(Tost.hhn.D.nonU)
TostSurface.D.nonU<-predictDsurface(Tost.hhn.D.nonU, se.D=TRUE, cl.D=TRUE)
plot(TostSurface.D.nonU,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Tost.cams,add=TRUE)
plotcovariate(TostSurface.D.nonU,covariate="stdGC",asp=1,contour=FALSE)
plot(Tost.cams,add=TRUE)

Nhat1D1.nonU<-region.N(Tost.hhn.D.nonU)
Nhat1D1.nonU



# Compare with and without non-Euclidian distance:
# -----------------------------------------------

# First save objects so don't have to refit:
save(Tost.hhn.D.nonU,Tost.hhn.DHab.nonU,Tost.hhn.DHab,file="./Tost/Tost-nonEuc-fits2.RData")
# load fitted objects:
load("./Tost/Tost-nonEuc-fits2.RData")

# Compare AICs:
AIC(Tost.hhn.D.nonU,Tost.hhn.DHab.nonU,Tost.hhn.DHab)

# get density range so plot on same scale
Dlim=range(covariates(TostSurface.nonU)$D.0,covariates(TostSurface)$D.0)
# Do plots next to each other:
windows() #Opens a separate window for plotting maps
par(mfrow=c(2,1))
plot.Dsurface(TostSurface,asp=1,contour=FALSE,col=terrain.colors(40),zlim=Dlim)
plot.Dsurface(TostSurface.nonU,asp=1,contour=FALSE,col=terrain.colors(40),zlim=Dlim)
# Hmmm, that is pretty unclear; try on log scale:
Dhat=TostSurface; covariates(Dhat)$D.0=log(covariates(Dhat)$D.0)
Dhat.nonU=TostSurface.nonU; covariates(Dhat.nonU)$D.0=log(covariates(Dhat.nonU)$D.0)
Dlim=range(covariates(Dhat.nonU)$D.0,covariates(Dhat)$D.0)
plot(Dhat,asp=1,contour=FALSE,col=terrain.colors(40),zlim=Dlim)
plot(Dhat.nonU,asp=1,contour=FALSE,col=terrain.colors(40),zlim=Dlim)

# Plot measure of probability of detection by given trap, from all points on mask:
lambda0=exp(coef(Tost.hhn.DHab.nonU)["lambda",1]) # on the real scale
sigma=exp(coef(Tost.hhn.DHab.nonU)["sigma",1]) # on the real scale
alpha=coef(Tost.hhn.DHab.nonU)["noneuc.stdGC",1] # on the beta scale
covariates(TostMask1)$noneuc=exp(-alpha*covariates(TostMask1)$stdGC)
# Do the plotting
par(mfrow=c(2,1))
j=6 # trap number
predicted.pj=pdot(TostMask1,Tost.cams[j,],noccasions=1,detectfn='HHN',
                  detectpar=list(lambda0=lambda0,sigma=sigma),userdist=userdfn1)
covariates(TostMask1)$predicted.pj=predicted.pj
plotcovariate(TostMask1,covariate="predicted.pj",contour=FALSE,asp=1)
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,col="white")) 
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,cex=0.25,col="black"))
# and plot stdGC below it:
plotcovariate(TostSurface.nonU,covariate="stdGC",asp=1,contour=FALSE)
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,col="white")) 
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,cex=0.25,col="black"))

j=16 # also try 32, 25, 39
predicted.pj=pdot(TostMask1,Tost.cams[j,],noccasions=1,detectfn='HHN',
                  detectpar=list(lambda0=lambda0,sigma=sigma),userdist=userdfn1)
covariates(TostMask1)$predicted.pj=predicted.pj
plotcovariate(TostMask1,covariate="predicted.pj",contour=FALSE,asp=1)
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,col="white")) 
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,cex=0.25,col="black"))
# and plot stdGC below it:
plotcovariate(TostSurface.nonU,covariate="stdGC",asp=1,contour=FALSE)
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,col="white")) 
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,cex=0.25,col="black"))

# Loop through a bunch by repeatedly running code below:
j=j+1;j
predicted.pj=pdot(TostMask1,Tost.cams[j,],noccasions=1,detectfn='HHN',
                  detectpar=list(lambda0=lambda0,sigma=sigma),userdist=userdfn1)
covariates(TostMask1)$predicted.pj=predicted.pj
plotcovariate(TostMask1,covariate="predicted.pj",contour=FALSE,asp=1)
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,col="white")) 
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,cex=0.25,col="black"))
# and plot stdGC below it:
plotcovariate(TostSurface.nonU,covariate="stdGC",asp=1,contour=FALSE)
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,col="white")) 
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,cex=0.25,col="black"))


# Look at constant-D non-Euclidian model:
# --------------------------------------
# Plot measure of probability of detection by given trap, from all points on mask:
lambda0=exp(coef(Tost.hhn.D.nonU)["lambda",1]) # on the real scale
sigma=exp(coef(Tost.hhn.D.nonU)["sigma",1]) # on the real scale
alpha=coef(Tost.hhn.D.nonU)["noneuc.stdGC",1] # on the beta scale
covariates(TostMask1)$noneuc=exp(-alpha*covariates(TostMask1)$stdGC)
# Do the plotting
par(mfrow=c(2,1))
j=1 # trap number
predicted.pj=pdot(TostMask1,Tost.cams[j,],noccasions=1,detectfn='HHN',
                  detectpar=list(lambda0=lambda0,sigma=sigma),userdist=userdfn1)
covariates(TostMask1)$predicted.pj=predicted.pj
plotcovariate(TostMask1,covariate="predicted.pj",contour=FALSE,asp=1)
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,col="white")) 
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,cex=0.25,col="black"))
# and plot stdGC below it:
plotcovariate(TostSurface.nonU,covariate="stdGC",asp=1,contour=FALSE)
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,col="white")) 
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,cex=0.25,col="black"))



# Test with different distance metric function userdfn2:
Tost.hhn.D.nonU2<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                          model=list(D~1, lambda0~1, sigma~1, 
                                     noneuc ~ stdGC -1), 
                          details = list(userdist = userdfn2),
                          start = list(noneuc = 1))

coefficients(Tost.hhn.D.nonU2)
TostSurface.D.nonU2<-predictDsurface(Tost.hhn.D.nonU2, se.D=TRUE, cl.D=TRUE)
plot(TostSurface.D.nonU2,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Tost.cams,add=TRUE)
plotcovariate(TostSurface.D.nonU2,covariate="stdGC",asp=1,contour=FALSE)
plot(Tost.cams,add=TRUE)

Nhat1D1.nonU2<-region.N(Tost.hhn.D.nonU2)
Nhat1D1.nonU2


# Look at constant-D non-Euclidian model with userdfn2:
# ----------------------------------------------------
# Plot measure of probability of detection by given trap, from all points on mask:
lambda0=exp(coef(Tost.hhn.D.nonU2)["lambda",1]) # on the real scale
sigma=exp(coef(Tost.hhn.D.nonU2)["sigma",1]) # on the real scale
alpha=coef(Tost.hhn.D.nonU2)["noneuc.stdGC",1] # on the beta scale
covariates(TostMask1)$noneuc=exp(alpha*covariates(TostMask1)$stdGC)
covariates(TostMask1)$noneuc=exp(covariates(TostMask1)$stdGC)
# Do the plotting
par(mfrow=c(2,1))
j=39 # trap number
predicted.pj=pdot(TostMask1,Tost.cams[j,],noccasions=1,detectfn='HHN',
                  detectpar=list(lambda0=lambda0,sigma=sigma),userdist=userdfn2)
covariates(TostMask1)$predicted.pj=predicted.pj
plotcovariate(TostMask1,covariate="predicted.pj",contour=FALSE,asp=1)
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,col="white")) 
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,cex=0.25,col="black"))
# and plot stdGC below it:
plotcovariate(TostMask1,covariate="stdGC",asp=1,contour=FALSE)
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,col="white")) 
plot(Tost.cams[j,], add = TRUE,detpar=list(pch=19,cex=0.25,col="black"))

