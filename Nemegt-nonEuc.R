library(secr)
library(fields)
library(maptools)
source("scrplotting.r")

#Running SECR for Nemegt 2013

# Read capture file and boundary
all.data.Nemegt<-read.capthist(captfile = "./Nemegt/Nemegt2013_Capture.csv", trapfile = "./Nemegt/Nemegt2013_Cams.csv", detector="count", fmt = "trapID", trapcovnames = c("Effort",	"Topo",	"Brokenness",	"Grass", "Rgd"))
boundaryNemegt=readShapeSpatial("./Nemegt//Habitat/Nemegt_StudyArea.shp")
# and plot it
plot(boundaryNemegt)
plot(x=all.data.Nemegt, add=TRUE)

# Make mask:
NemegtMask=make.mask(traps(all.data.Nemegt), spacing=500, buffer = 25000, type="trapbuffer", poly=boundaryNemegt)

# Read ruggedness covariate and put it into mask covariate GRIDCODE
SLCost.Nemegt<-readShapePoly("./Nemegt//Habitat/Nemegt_Rgd500m.shp")  #ruggedness pixels averaged over 500m radius
NemegtMask1<-addCovariates(NemegtMask, SLCost.Nemegt)

# Read binary habitat suitability code into mask covariate ...
SLCostBINARY.Nemegt<-readShapePoly("./Tost//Habitat/tost_sl.shp")  #Logistic binary SL habitat created using telemetry data
plot(SLCostBINARY.Nemegt, add=TRUE)

NemegtMask1<-addCovariates(NemegtMask1, SLCostBINARY.Nemegt)
head(covariates(NemegtMask1))
names(covariates(NemegtMask1))[3:4] = c("binaryID","BINCODE") #Rename headers
summary(covariates(NemegtMask1))
# make NAs in BINCODE zeros:
covariates(NemegtMask1)$BINCODE[is.na(covariates(NemegtMask1)$BINCODE)] = 0
summary(covariates(NemegtMask1))


# Standarize Rgd on traps (this makes fits a bit more stable)
# -----------------------------------------------------------
summary(covariates(traps(all.data.Nemegt)))
covariates(traps(all.data.Nemegt))$stdRgd = scale(covariates(traps(all.data.Nemegt))$Rgd)
summary(covariates(traps(all.data.Nemegt)))
head(covariates(traps(all.data.Nemegt)))

# Standarize GRIDCODE (in stdGC) and BINCODE (in stdBC) on mask
# ------------------------------------------------------------------------
summary(covariates(NemegtMask1))
covariates(NemegtMask1)$stdGC = scale(covariates(NemegtMask1)$GRIDCODE)
covariates(NemegtMask1)$stdBC = scale(covariates(NemegtMask1)$BINCODE)
summary(covariates(NemegtMask1))
names(covariates(NemegtMask1))


head(covariates(NemegtMask1))
summary(covariates(NemegtMask1))
summary(covariates(traps(all.data.Nemegt)))

plot(NemegtMask1, covariate="stdGC", contour=FALSE, col=terrain.colors(10), legend = FALSE)
head(covariates(traps(all.data.Nemegt)))
head(covariates(NemegtMask1))


Nemegt.cams=traps(all.data.Nemegt)

Nemegt.hhn<-secr.fit(all.data.Nemegt, model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", mask=NemegtMask1)
Nemegt.hhn.detrgd<-secr.fit(all.data.Nemegt, model=list(D~1, lambda0~stdRgd, sigma~stdRgd), detectfn="HHN", mask=NemegtMask1)
Nemegt.hhn.DHab<-secr.fit(all.data.Nemegt, model=list(D~stdGC, lambda0~1, sigma~1), detectfn="HHN", mask=NemegtMask1)

#AIC(Noyon.hhn, Noyon.hhn.detrgd, Noyon.hhn.DHab)
coefficients(Nemegt.hhn.DHab)
NemegtSurface<-predictDsurface(Nemegt.hhn.DHab, se.D=TRUE, cl.D=TRUE)
plot(NemegtSurface,asp=1,contour=FALSE)
plotcovariate(NemegtSurface,covariate="stdGC",asp=1,contour=FALSE)
plot(Nemegt.cams,add=TRUE)

Nhat1<-region.N(Nemegt.hhn.DHab) #Estimates the population N of the animals within the region defined by mask
Nhat1

Nhat2<-region.N(Nemegt.hhn)
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
Nemegt.hhn.DHab.nonU<-secr.fit(all.data.Nemegt, detectfn="HHN", mask=NemegtMask1,
                              model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                              details = list(userdist = userdfn1),
                              start = list(noneuc = 1)) #-1 gets rid of the intercept
AIC(Nemegt.hhn,Nemegt.hhn.DHab.nonU)
coefficients(Nemegt.hhn.DHab.nonU)
NemegtSurface.nonU<-predictDsurface(Nemegt.hhn.DHab.nonU, se.D=TRUE, cl.D=TRUE)
windows()
plot(NemegtSurface.nonU,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Nemegt.cams,add=TRUE)
plotcovariate(NemegtSurface.nonU,covariate="stdGC",asp=1,contour=FALSE)
plot(Nemegt.cams,add=TRUE)

Nhat1.nonU<-region.N(Nemegt.hhn.DHab.nonU)
Nhat1.nonU

# Model with Flat density and non-Euclidian:
#-------------------------------------------
Nemegt.hhn.D.nonU<-secr.fit(all.data.Nemegt, detectfn="HHN", mask=NemegtMask1,
                           model=list(D~1, lambda0~1, sigma~1, 
                                      noneuc ~ stdGC -1), 
                           details = list(userdist = userdfn1),
                           start = list(noneuc = 1))

coefficients(Nemegt.hhn.D.nonU)
NemegtSurface.D.nonU<-predictDsurface(Nemegt.hhn.D.nonU, se.D=TRUE, cl.D=TRUE)
plot(NemegtSurface.D.nonU,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Nemegt.cams,add=TRUE)
plotcovariate(NemegtSurface.D.nonU,covariate="stdGC",asp=1,contour=FALSE)
plot(Nemegt.cams,add=TRUE)

Nhat1D1.nonU<-region.N(Nemegt.hhn.D.nonU)
Nhat1D1.nonU
#AIC(Nemegt.hhn,Nemegt.hhn.D.nonU,Nemegt.hhn.DHab.nonU)

# Model with stdGC and stdBC in noneuc:
# -------------------------------------
Nemegt.hhn.DHab.nonU.GBGC<-secr.fit(all.data.Nemegt, detectfn="HHN", mask=NemegtMask1,
                                   model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC + stdBC-1), 
                                   details = list(userdist = userdfn1),
                                   start = list(noneuc = 1))
coefficients(Nemegt.hhn.DHab.nonU.GBGC)

NemegtSurface.nonU.GBGC<-predictDsurface(Nemegt.hhn.DHab.nonU.GBGC, se.D=TRUE, cl.D=TRUE)
plot(NemegtSurface.nonU.GBGC,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Nemegt.cams,add=TRUE)
plotcovariate(NemegtSurface.nonU.GBGC,covariate="stdGC",asp=1,contour=FALSE)
plot(Nemegt.cams,add=TRUE)

Nhat1.nonU<-region.N(Nemegt.hhn.DHab.nonU.GBGC)
Nhat1.nonU
AIC(Nemegt.hhn,Nemegt.hhn.D.nonU,Nemegt.hhn.DHab.nonU, Nemegt.hhn.DHab.nonU.GBGC)

# Model with stdBC only in noneuc:
# -------------------------------------
Nemegt.hhn.DHab.nonU.GB<-secr.fit(all.data.Nemegt, detectfn="HHN", mask=NemegtMask1,
                                 model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~stdBC-1), 
                                 details = list(userdist = userdfn1),
                                 start = list(noneuc = 1))
coefficients(Nemegt.hhn.DHab.nonU.GB)
NemegtSurface.nonU.GB<-predictDsurface(Nemegt.hhn.DHab.nonU.GB, se.D=TRUE, cl.D=TRUE)
plot(NemegtSurface.nonU.GB,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Nemegt.cams,add=TRUE)
plotcovariate(NemegtSurface.nonU.GB,covariate="stdGC",asp=1,contour=FALSE)
plot(Nemegt.cams,add=TRUE)

Nhat1.nonU<-region.N(Nemegt.hhn.DHab.nonU.GB)
Nhat1.nonU

# Compare models with and without non-Euclidian distance:
# -----------------------------------------------

# Save objects:
NemegtSurface<-predictDsurface(Nemegt.hhn.DHab, se.D=TRUE, cl.D=TRUE)
NemegtSurface.nonU<-predictDsurface(Nemegt.hhn.DHab.nonU, se.D=TRUE, cl.D=TRUE)
NemegtSurface.D.nonU<-predictDsurface(Nemegt.hhn.D.nonU, se.D=TRUE, cl.D=TRUE)
NemegtSurface.DHab.nonU.GB<-predictDsurface(Nemegt.hhn.DHab.nonU.GB, se.D=TRUE, cl.D=TRUE)
NemegtSurface.DHab.nonU.GBGC<-predictDsurface(Nemegt.hhn.DHab.nonU.GBGC, se.D=TRUE, cl.D=TRUE)
save(Nemegt.cams,NemegtMask1,NemegtSurface,NemegtSurface.nonU,NemegtSurface.D.nonU,Nemegt.hhn.DHab,
     Nemegt.hhn,Nemegt.hhn.detrgd,Nemegt.hhn.D.nonU,Nemegt.hhn.DHab.nonU,Nemegt.hhn.DHab.nonU.GB,
     Nemegt.hhn.DHab.nonU.GBGC,Nemegt.hhn.DHab,file="./Nemegt/Nemegt-nonEuc-fits.RData")
# load fitted objects:
load("./Nemegt/Nemegt-nonEuc-fits.RData")

# Compare all model AICs:
AIC(Nemegt.hhn, Nemegt.hhn.detrgd,Nemegt.hhn.DHab,Nemegt.hhn.DHab.nonU,Nemegt.hhn.D.nonU, 
    Nemegt.hhn.DHab.nonU.GBGC, Nemegt.hhn.DHab.nonU.GB)

# get density range so plot on same scale
Dlim=range(covariates(NemegtSurface.nonU)$D.0,covariates(NemegtSurface)$D.0)

# Do plots next to each other:
par(mfrow=c(2,1))
plot.Dsurface(NemegtSurface,asp=1,contour=FALSE,col=terrain.colors(40),zlim=Dlim)
plot(NemegtSurface.nonU,asp=1,contour=FALSE,col=terrain.colors(40),zlim=Dlim)
# Using log scale to plot clearly:
Dhat=NemegtSurface; covariates(Dhat)$D.0=log(covariates(Dhat)$D.0)
Dhat.nonU=NemegtSurface.nonU; covariates(Dhat.nonU)$D.0=log(covariates(Dhat.nonU)$D.0)
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
quartz(h=5,w=10) #DID NOT FIND THIS FUNCTION!!!
plotcovariate(NemegtSurface.D.nonU,covariate="stdGC",asp=1,contour=FALSE,col=terrain.colors(40))
text(Nemegt.cams,labels=as.character(1:40),cex=0.75)