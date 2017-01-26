library(secr)
library(fields)
library(maptools)
source("scrplotting.r")


#Running SECR for Tost 2012

# Read capture file and boundary
all.data.Tost<-read.capthist(captfile = "./Tost/Tost_capthist2012.csv", trapfile = "./Tost/Tost_cams_rugged2012.csv", detector="count", fmt = "trapID", trapcovnames = c("Effort",	"Topo",	"Altidute",	"Rgd"))
boundaryTost=readShapeSpatial("./Tost//Habitat/TostStudy_Area.shp")
# and plot it
plot(boundaryTost)
plot(x=all.data.Tost, add=TRUE)

# Make mask:
TostMask=make.mask(traps(all.data.Tost), spacing=500, buffer = 25000, type="trapbuffer", poly=boundaryTost)

# Read ruggedness covariate and put it into mask covariate GRIDCODE
SLCost.Tost<-readShapePoly("./Tost//Habitat/Tost_Rgd500m.shp")  #ruggedness pixels averaged over 500m radius
TostMask1<-addCovariates(TostMask, SLCost.Tost)

# Read binary habitat suitability code into mask covariate ...
SLCostBINARY.Tost<-readShapePoly("./Tost//Habitat/tost_sl.shp")  #ruggedness pixels averaged over 500m radius
TostMask1<-addCovariates(TostMask1, SLCostBINARY.Tost)
names(covariates(TostMask1))[3:4] = c("binaryID","BINCODE")
summary(covariates(TostMask1))
# make NAs in BINCODE zeros:
covariates(TostMask1)$BINCODE[is.na(covariates(TostMask1)$BINCODE)] = 0
summary(covariates(TostMask1))


# Standarize Rgd on traps (this makes fits a bit more stable)
# -----------------------------------------------------------
summary(covariates(traps(all.data.Tost)))
covariates(traps(all.data.Tost))$stdRgd = scale(covariates(traps(all.data.Tost))$Rgd)
summary(covariates(traps(all.data.Tost)))
head(covariates(traps(all.data.Tost)))

# Standarize GRIDCODE (in stdGC) and BINCODE (in stdBC) on mask
# ------------------------------------------------------------------------
summary(covariates(TostMask1))
covariates(TostMask1)$stdGC = scale(covariates(TostMask1)$GRIDCODE)
covariates(TostMask1)$stdBC = scale(covariates(TostMask1)$BINCODE)
summary(covariates(TostMask1))
names(covariates(TostMask1))


head(covariates(TostMask1))
summary(covariates(TostMask1))
summary(covariates(traps(all.data.Tost)))

plot(TostMask1, covariate="stdGC", contour=FALSE, col=terrain.colors(10), legend = FALSE)
head(covariates(traps(all.data.Tost)))
head(covariates(TostMask1))


Tost.cams=traps(all.data.Tost)

Tost.hhn<-secr.fit(all.data.Tost, model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", mask=TostMask1)
#Tost.hhn.detrug<-secr.fit(all.data.Tost, model=list(D~1, lambda0~Rgd, sigma~Rgd), detectfn="HHN", mask=TostMask1)
#Tost.hhn.Dxy<-secr.fit(all.data.Tost, model=list(D~x+y, lambda0~1, sigma~1), detectfn="HHN", mask=TostM)
Tost.hhn.DHab<-secr.fit(all.data.Tost, model=list(D~stdGC, lambda0~1, sigma~1), detectfn="HHN", mask=TostMask1)
#Tost.hhn.Dx<-secr.fit(all.data.Tost, model=list(D~x, lambda0~1, sigma~1), detectfn="HHN", mask=TostMask1)

AIC(Tost.hhn, Tost.hhn.DHab)
coefficients(Tost.hhn.DHab)
TostSurface<-predictDsurface(Tost.hhn.DHab, se.D=TRUE, cl.D=TRUE)
plot(TostSurface,asp=1,contour=FALSE)
plotcovariate(TostSurface,covariate="stdGC",asp=1,contour=FALSE)
plot(Tost.cams,add=TRUE)

Nhat1<-region.N(Tost.hhn.DHab) #Estimates the population N of the animals within the region defined by mask
Nhat1

Nhat2<-region.N(Tost.hhn)
Nhat2

# Non-Euclidian fits
# ==================
# This taken straight from secr vignette:
userdfn1 <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') #When function is called, more of a jargon. Tells that it is a non-euclidean function
  require(gdistance) #to load transition and geoCorrection functions
  Sraster <- raster(mask, 'noneuc') #Creates a raster from a set of coordinates and attributes and turn that into a raster. noneuc needs to be made in advance in the mask that is being used in the analysis
  ## conductance is inverse of friction
  trans <- transition(Sraster, transitionFunction = function(x) 1/mean(x),  directions = 16)
  trans <- geoCorrection(trans) #takes care of earth's curvature and also the distance differences between square and diagonally neighbouring cells
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

# Model with stdGC in noneuc:
# ---------------------------
Tost.hhn.DHab.nonU<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                             model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                             details = list(userdist = userdfn1),
                             start = list(noneuc = 1)) #-1 gets rid of the intercept
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

# Model with stdGC and stdBC in noneuc:
# -------------------------------------
Tost.hhn.DHab.nonU.GBGC<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                             model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC + stdBC-1), 
                             details = list(userdist = userdfn1),
                             start = list(noneuc = 1))
coefficients(Tost.hhn.DHab.nonU.GBGC)
TostSurface.nonU<-predictDsurface(Tost.hhn.DHab.nonU.GBGC, se.D=TRUE, cl.D=TRUE)
plot(TostSurface.nonU.GBGC,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Tost.cams,add=TRUE)
plotcovariate(TostSurface.nonU.GBGC,covariate="stdGC",asp=1,contour=FALSE)
plot(Tost.cams,add=TRUE)

Nhat1.nonU.GBGC<-region.N(Tost.hhn.DHab.nonU.GBGC)
Nhat1.nonU.GBGC

# Model with stdBC in noneuc:
# -------------------------------------
Tost.hhn.DHab.nonU.GB<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~stdBC-1), 
                                details = list(userdist = userdfn1),
                                start = list(noneuc = 1))
coefficients(Tost.hhn.DHab.nonU.GB)
TostSurface.nonU<-predictDsurface(Tost.hhn.DHab.nonU.GB, se.D=TRUE, cl.D=TRUE)
plot(TostSurface.nonU.GB,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Tost.cams,add=TRUE)
plotcovariate(TostSurface.nonU.GB,covariate="stdGC",asp=1,contour=FALSE)
plot(Tost.cams,add=TRUE)

Nhat1.nonU.GB<-region.N(Tost.hhn.DHab.nonU.GB)
Nhat1.nonU.GB


# Compare with and without non-Euclidian distance:
# -----------------------------------------------

# First save objects so don't have to refit:
TostSurface<-predictDsurface(Tost.hhn.DHab, se.D=TRUE, cl.D=TRUE)
TostSurface.nonU<-predictDsurface(Tost.hhn.DHab.nonU, se.D=TRUE, cl.D=TRUE)
TostSurface.D.nonU<-predictDsurface(Tost.hhn.D.nonU, se.D=TRUE, cl.D=TRUE)
TostSurface.DHab.nonU.GB<-predictDsurface(Tost.hhn.DHab.nonU.GB, se.D=TRUE, cl.D=TRUE)
TostSurface.DHab.nonU.GBGC<-predictDsurface(Tost.hhn.DHab.nonU.GBGC, se.D=TRUE, cl.D=TRUE)
save(Tost.cams,TostMask1,TostSurface,TostSurface.nonU,TostSurface.D.nonU,Tost.hhn.D.nonU,
     Tost.hhn.DHab.nonU,Tost.hhn.DHab.nonU.GB,Tost.hhn.DHab.nonU.GBGC,Tost.hhn.DHab,file="./Tost/Tost-nonEuc-fits2.RData")
# load fitted objects:
load("./Tost/Tost-nonEuc-fits2.RData")

# Compare AICs:
AIC(Tost.hhn.D.nonU,Tost.hhn.DHab.nonU,Tost.hhn.DHab,Tost.hhn.DHab.nonU.GB,Tost.hhn.DHab.nonU.GBGC)


# Best non-Euclidian
Tost.Nhatbest.nonU<-region.N(Tost.hhn.DHab.nonU)
# Best Euclidian
Tost.Nhatbest.U<-region.N(Tost.hhn.DHab)
# Compare them:
Tost.Nhatbest.nonU
Tost.Nhatbest.U

# get density range so plot on same scale
Dlim=range(covariates(TostSurface.nonU)$D.0,covariates(TostSurface)$D.0)

# Do plots next to each other:
par(mfrow=c(2,1))
plot.Dsurface(TostSurface,asp=1,contour=FALSE,col=terrain.colors(40),zlim=Dlim)
plot(TostSurface.nonU,asp=1,contour=FALSE,col=terrain.colors(40),zlim=Dlim)
# Hmmm, that is pretty unclear; try on log scale:
Dhat=TostSurface; covariates(Dhat)$D.0=log(covariates(Dhat)$D.0)
Dhat.nonU=TostSurface.nonU; covariates(Dhat.nonU)$D.0=log(covariates(Dhat.nonU)$D.0)
Dlim=range(covariates(Dhat.nonU)$D.0,covariates(Dhat)$D.0)
plot(Dhat,asp=1,contour=FALSE,col=terrain.colors(40),zlim=Dlim)
plot(Dhat.nonU,asp=1,contour=FALSE,col=terrain.colors(40),zlim=Dlim)

# Plot non-euclidian distance map, selecting origin with cursor.
dmap <- function (traps, mask, userd, i = 1, ...) {
  if (is.na(i)) i <- nearesttrap(unlist(locator(1)), traps) 
  covariates(mask) <- data.frame(d = userd[i,]) 
  covariates(mask)$d[!is.finite(covariates(mask)$d)] <- NA 
  plotcovariate(mask, covariate = 'd',  ...)
  points(traps[i,], pch = 3, col = 'red')
}
# plot gridcode
windows(h=5,w=10)
plotcovariate(TostSurface.D.nonU,covariate="stdGC",asp=1,contour=FALSE,col=terrain.colors(40))
text(Tost.cams,labels=as.character(1:40),cex=0.75)
# plot distance from given trap
trap=7 # choose trap
dmask=TostMask1
covariates(dmask)$noneuc = log(covariates(TostMask1)$noneuc +1) # log for better plot
im=prep4image(data.frame(x=dmask$x,y=dmask$y,z=covariates(dmask)$stdGC),plot=FALSE)
quartz(h=5,w=10)
crfun=function(x) sqrt(prod(x))
userd = nedist(Tost.cams,dmask,dmask,transitionFunction=crfun)
dmap(Tost.cams, dmask, userd, i=trap, contour=FALSE)
dmap(Tost.cams, dmask, userd, i=NA, contour=FALSE)
points(Tost.cams[trap,],pch=19,cex=2,col="white")
points(Tost.cams[trap,],pch=19,cex=1.5)
contour(im,add=TRUE,col="gray",lwd=0.5,nlevels=6)
text(Tost.cams,labels=as.character(1:40),cex=0.75,col="white")


# Plot measure of probability of detection by given trap, from all points on mask:
lambda0=exp(coef(Tost.hhn.DHab.nonU)["lambda",1]) # on the real scale
sigma=exp(coef(Tost.hhn.DHab.nonU)["sigma",1]) # on the real scale
alpha=coef(Tost.hhn.DHab.nonU)["noneuc.stdGC",1] # on the beta scale
covariates(TostMask1)$noneuc=exp(-alpha*covariates(TostMask1)$stdGC)
# Do the plotting
par(mfrow=c(2,1))
j=7 # trap number
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

j=32 # also try 16, 25, 39
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

