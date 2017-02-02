library(secr)
library(fields)
library(maptools)
source("scrplotting.r")

#Running SECR for Noyon 2013

# Read capture file and boundary
all.data.Noyon<-read.capthist(captfile = "./Noyon2013/Noyon_capthist2013secr.csv", trapfile = "./Noyon2013/Noyon_trap2013secr.csv", detector="count", fmt = "trapID", trapcovnames = c("Effort",	"Topo",	"Substrate",	"Brokenness", "Rgd", "Water"))
boundaryNoyon=readShapeSpatial("./Noyon2013/Habitat/NoyonStudy_Area.shp")
# and plot it
plot(boundaryNoyon)
plot(x=all.data.Noyon, add=TRUE)

# Make mask:
NoyonMask=make.mask(traps(all.data.Noyon), spacing=500, buffer = 25000, type="trapbuffer", poly=boundaryNoyon)

# Read ruggedness covariate and put it into mask covariate GRIDCODE
SLCost.Noyon<-readShapePoly("./Noyon2013//Habitat/Noyon_Rgd500m.shp")  #ruggedness pixels averaged over 500m radius
NoyonMask1<-addCovariates(NoyonMask, SLCost.Noyon)
head(covariates(NoyonMask1))
# Read binary habitat suitability code into mask covariate ...
SLCostBINARY.Noyon<-readShapePoly("./Tost//Habitat/tost_sl.shp")  #Logistic binary SL habitat created using telemetry data
plot(SLCostBINARY.Noyon, add=TRUE)

NoyonMask1<-addCovariates(NoyonMask1, SLCostBINARY.Noyon)
head(covariates(NoyonMask1))
names(covariates(NoyonMask1))[3:4] = c("binaryID","BINCODE") #Rename headers
summary(covariates(NoyonMask1))
# make NAs in BINCODE zeros:
covariates(NoyonMask1)$BINCODE[is.na(covariates(NoyonMask1)$BINCODE)] = 0
summary(covariates(NoyonMask1))


# Standarize Rgd on traps (this makes fits a bit more stable)
# -----------------------------------------------------------
summary(covariates(traps(all.data.Noyon)))
covariates(traps(all.data.Noyon))$stdRgd = scale(covariates(traps(all.data.Noyon))$Rgd)
summary(covariates(traps(all.data.Noyon)))
head(covariates(traps(all.data.Noyon)))

# Standarize GRIDCODE (in stdGC) and BINCODE (in stdBC) on mask
# ------------------------------------------------------------------------
summary(covariates(NoyonMask1))
covariates(NoyonMask1)$stdGC = scale(covariates(NoyonMask1)$GRIDCODE)
covariates(NoyonMask1)$stdBC = scale(covariates(NoyonMask1)$BINCODE)
summary(covariates(NoyonMask1))
names(covariates(NoyonMask1))


head(covariates(NoyonMask1))
summary(covariates(NoyonMask1))
summary(covariates(traps(all.data.Noyon)))

plot(NoyonMask1, covariate="stdGC", contour=FALSE, col=terrain.colors(10), legend = FALSE)
head(covariates(traps(all.data.Noyon)))
head(covariates(NoyonMask1))


Noyon.cams=traps(all.data.Noyon)

Noyon.hhn<-secr.fit(all.data.Noyon, model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", mask=NoyonMask1)
Noyon.hhn.detrgd<-secr.fit(all.data.Noyon, model=list(D~1, lambda0~stdRgd, sigma~stdRgd), detectfn="HHN", mask=NoyonMask1)
Noyon.hhn.DHab<-secr.fit(all.data.Noyon, model=list(D~stdGC, lambda0~1, sigma~1), detectfn="HHN", mask=NoyonMask1)
Noyon.hhn.DHab.DetRgd01<-secr.fit(all.data.Noyon, model=list(D~stdGC, lambda0~1, sigma~stdRgd), detectfn="HHN",mask=NoyonMask1)
Noyon.hhn.detWater<-secr.fit(all.data.Noyon, model=list(D~1, lambda0~Water, sigma~1), detectfn="HHN", mask=NoyonMask1)
Noyon.hhn.detTopo10<-secr.fit(all.data.Noyon, model=list(D~1, lambda0~Topo, sigma~1), detectfn="HHN", mask=NoyonMask1)
Noyon.hhn.detTopoWater<-secr.fit(all.data.Noyon, model=list(D~1, lambda0~Topo+Water, sigma~1), detectfn="HHN", mask=NoyonMask1)

AIC(Noyon.hhn, Noyon.hhn.detrgd, Noyon.hhn.DHab, Noyon.hhn.DHab.DetRgd01, Noyon.hhn.detWater, 
    Noyon.hhn.detTopo10, Noyon.hhn.detTopoWater)
coefficients(Noyon.hhn.DHab.DetRgd01)

NoyonSurface<-predictDsurface(Noyon.hhn.DHab.DetRgd01, se.D=TRUE, cl.D=TRUE)
plot(NoyonSurface,asp=1,contour=FALSE)
plotcovariate(NoyonSurface,covariate="stdGC",asp=1,contour=FALSE)
plot(Noyon.cams,add=TRUE)

Nhat1<-region.N(Noyon.hhn.DHab.DetRgd01) #Estimates the population N of the animals within the region defined by mask
Nhat1

Nhat2<-region.N(Noyon.hhn)
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

# Model with stdGC in noneuc:
# ---------------------------
Noyon.hhn.DHab.nonU<-secr.fit(all.data.Noyon, detectfn="HHN", mask=NoyonMask1,
                             model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                             details = list(userdist = userdfn1),
                             start = list(noneuc = 1)) #-1 gets rid of the intercept
AIC(Noyon.hhn,Noyon.hhn.DHab.nonU)
coefficients(Noyon.hhn.DHab.nonU)
NoyonSurface.nonU<-predictDsurface(Noyon.hhn.DHab.nonU, se.D=TRUE, cl.D=TRUE)
windows()
plot(NoyonSurface.nonU,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Noyon.cams,add=TRUE)
plotcovariate(NoyonSurface.nonU,covariate="stdGC",asp=1,contour=FALSE)
plot(Noyon.cams,add=TRUE)

Nhat1.nonU<-region.N(Noyon.hhn.DHab.nonU)
Nhat1.nonU

# Model with Flat density and non-Euclidian:
#-------------------------------------------
Noyon.hhn.D.nonU<-secr.fit(all.data.Noyon, detectfn="HHN", mask=NoyonMask1,
                          model=list(D~1, lambda0~1, sigma~1, 
                                     noneuc ~ stdGC -1), 
                          details = list(userdist = userdfn1),
                          start = list(noneuc = 1))

coefficients(Noyon.hhn.D.nonU)
NoyonSurface.D.nonU<-predictDsurface(Noyon.hhn.D.nonU, se.D=TRUE, cl.D=TRUE)
plot(NoyonSurface.D.nonU,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Noyon.cams,add=TRUE)
plotcovariate(NoyonSurface.D.nonU,covariate="stdGC",asp=1,contour=FALSE)
plot(Noyon.cams,add=TRUE)

Nhat1D1.nonU<-region.N(Noyon.hhn.D.nonU)
Nhat1D1.nonU
AIC(Noyon.hhn,Noyon.hhn.D.nonU,Noyon.hhn.DHab.nonU)

# Model with stdGC and stdBC in noneuc:
# -------------------------------------
Noyon.hhn.DHab.nonU.GBGC<-secr.fit(all.data.Noyon, detectfn="HHN", mask=NoyonMask1,
                                  model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC + stdBC-1), 
                                  details = list(userdist = userdfn1),
                                  start = list(noneuc = 1))
coefficients(Noyon.hhn.DHab.nonU.GBGC)

NoyonSurface.nonU<-predictDsurface(Noyon.hhn.DHab.nonU.GBGC, se.D=TRUE, cl.D=TRUE)
plot(NoyonSurface.nonU.GBGC,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Noyon.cams,add=TRUE)
plotcovariate(NoyonSurface.nonU.GBGC,covariate="stdGC",asp=1,contour=FALSE)
plot(Noyon.cams,add=TRUE)

Nhat1.nonU<-region.N(Noyon.hhn.DHab.nonU.GBGC)
Nhat1.nonU
AIC(Noyon.hhn,Noyon.hhn.D.nonU,Noyon.hhn.DHab.nonU, Noyon.hhn.DHab.nonU.GBGC)

# Model with stdBC only in noneuc:
# -------------------------------------
Noyon.hhn.DHab.nonU.GB<-secr.fit(all.data.Noyon, detectfn="HHN", mask=NoyonMask1,
                                model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~stdBC-1), 
                                details = list(userdist = userdfn1),
                                start = list(noneuc = 1))
coefficients(Noyon.hhn.DHab.nonU.GB)
NoyonSurface.nonU.GB<-predictDsurface(Noyon.hhn.DHab.nonU.GB, se.D=TRUE, cl.D=TRUE)
plot(NoyonSurface.nonU.GB,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Noyon.cams,add=TRUE)
plotcovariate(NoyonSurface.nonU.GB,covariate="stdGC",asp=1,contour=FALSE)
plot(Noyon.cams,add=TRUE)

Nhat1.nonU<-region.N(Noyon.hhn.DHab.nonU.GB)
Nhat1.nonU

# Model with stdGC in noneuc Topo in Detection:
# ---------------------------
Noyon.hhn.DHab.Topo10.nonU<-secr.fit(all.data.Noyon, detectfn="HHN", mask=NoyonMask1,
                              model=list(D~stdGC, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                              details = list(userdist = userdfn1),
                              start = list(noneuc = 1)) #-1 gets rid of the intercept

# Model with stdGC in noneuc Water in Detection:
# ---------------------------
Noyon.hhn.DHab.DetW.nonU<-secr.fit(all.data.Noyon, detectfn="HHN", mask=NoyonMask1,
                                     model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                                     details = list(userdist = userdfn1),
                                     start = list(noneuc = 1)) #-1 gets rid of the intercept

# Model with stdGC in noneuc Topo & Water in Detection:
# ---------------------------
Noyon.hhn.DHab.Topo10W.nonU<-secr.fit(all.data.Noyon, detectfn="HHN", mask=NoyonMask1,
                                     model=list(D~stdGC, lambda0~Water+Topo, sigma~1, noneuc ~ stdGC -1), 
                                     details = list(userdist = userdfn1),
                                     start = list(noneuc = 1)) #-1 gets rid of the intercept

AIC(Noyon.hhn,Noyon.hhn.D.nonU,Noyon.hhn.DHab.nonU, Noyon.hhn.DHab.Topo10.nonU, Noyon.hhn.DHab.DetW.nonU,
    Noyon.hhn.DHab.Topo10W.nonU)

# Compare models with and without non-Euclidian distance:
# -----------------------------------------------

# Save objects:
NoyonSurface<-predictDsurface(Noyon.hhn.DHab, se.D=TRUE, cl.D=TRUE)
NoyonSurface.nonU<-predictDsurface(Noyon.hhn.DHab.nonU, se.D=TRUE, cl.D=TRUE)
NoyonSurface.D.nonU<-predictDsurface(Noyon.hhn.D.nonU, se.D=TRUE, cl.D=TRUE)
NoyonSurface.DHab.nonU.GB<-predictDsurface(Noyon.hhn.DHab.nonU.GB, se.D=TRUE, cl.D=TRUE)
NoyonSurface.DHab.nonU.GBGC<-predictDsurface(Noyon.hhn.DHab.nonU.GBGC, se.D=TRUE, cl.D=TRUE)
save(Noyon.cams,NoyonMask1,NoyonSurface,NoyonSurface.nonU,NoyonSurface.D.nonU,Noyon.hhn.DHab,
      Noyon.hhn,Noyon.hhn.detrgd,Noyon.hhn.D.nonU,Noyon.hhn.DHab.nonU,Noyon.hhn.DHab.nonU.GB,
      Noyon.hhn.DHab.nonU.GBGC,Noyon.hhn.DHab, Noyon.hhn.DHab.Topo10.nonU, Noyon.hhn.DHab.DetW.nonU,
     Noyon.hhn.DHab.Topo10W.nonU, file="./Noyon2013/Noyon-nonEuc-fits3.RData")
# load fitted objects:
load("./Noyon2013/Noyon-nonEuc-fits.RData")

load("./Noyon2013/Noyon-nonEuc-fits.RData") #Second round analysis (with water and topography)

# Compare all model AICs:
AIC(Noyon.hhn, Noyon.hhn.detrgd,Noyon.hhn.DHab,Noyon.hhn.DHab.nonU,Noyon.hhn.D.nonU, 
    Noyon.hhn.DHab.nonU.GBGC, Noyon.hhn.DHab.nonU.GB)

# get density range so plot on same scale
Dlim=range(covariates(NoyonSurface.nonU)$D.0,covariates(NoyonSurface)$D.0)

# Do plots next to each other:
par(mfrow=c(2,1))
plot.Dsurface(NoyonSurface,asp=1,contour=FALSE,col=terrain.colors(40),zlim=Dlim)
plot(NoyonSurface.nonU,asp=1,contour=FALSE,col=terrain.colors(40),zlim=Dlim)
# Using log scale to plot clearly:
Dhat=NoyonSurface; covariates(Dhat)$D.0=log(covariates(Dhat)$D.0)
Dhat.nonU=NoyonSurface.nonU; covariates(Dhat.nonU)$D.0=log(covariates(Dhat.nonU)$D.0)
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
windows(h=5,w=10) #DID NOT FIND THIS FUNCTION!!!
plotcovariate(NoyonSurface.D.nonU,covariate="stdGC",asp=1,contour=FALSE,col=terrain.colors(40))
text(Noyon.cams,labels=as.character(1:40),cex=0.75)

# plot distance from given trap
##for some reason, error stating replacement has 0 rows, data has 9927 at row 215!
trap=7 # choose trap
dmask=NoyonMask1
head(covariates(NoyonMask1))
names(covariates(NoyonMask1))
covariates(NoyonMask1)$noneuc=exp(-alpha*covariates(NoyonMask1)$stdGC)
covariates(dmask)$noneuc = log(covariates(NoyonMask1)$noneuc +1) # log for better plot
im=prep4image(data.frame(x=dmask$x,y=dmask$y,z=covariates(dmask)$stdGC),plot=FALSE)
windows(h=5,w=10)
crfun=function(x) sqrt(prod(x))
userd = nedist(Tost.cams,dmask,dmask,transitionFunction=crfun)
dmap(Tost.cams, dmask, userd, i=trap, contour=FALSE)
dmap(Tost.cams, dmask, userd, i=NA, contour=FALSE)
points(Tost.cams[trap,],pch=19,cex=2,col="white")
points(Tost.cams[trap,],pch=19,cex=1.5)
contour(im,add=TRUE,col="gray",lwd=0.5,nlevels=6)
text(Tost.cams,labels=as.character(1:40),cex=0.75,col="white")

