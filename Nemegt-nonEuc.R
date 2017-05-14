library(secr)
library(fields)
library(maptools)
source("scrplotting.r")

#Running SECR for Nemegt 2013

# Read capture file and boundary
all.data.Nemegt<-read.capthist(captfile = "./Nemegt/Nemegt2013_Capture.csv", 
                               trapfile = "./Nemegt/Nemegt2013_Cams.csv", 
                               detector="count", fmt = "trapID", 
                               trapcovnames = c("Topo",	"Brokenness",	"Grass", "Rgd", "Water"),
                               binary.usage=FALSE)
#Read capture file-reduced
all.data.Nemegt_R<-read.capthist(captfile = "./Nemegt/Nemegt2013_Capture_reduced.csv", 
                               trapfile = "./Nemegt/Nemegt2013_Cams_add.csv", 
                               detector="count", fmt = "trapID", 
                               trapcovnames = c("Topo",	"Brokenness",	"Grass", "Rgd", "Water"),
                               binary.usage=FALSE)

#Read capture file-clubbed
all.data.Nemegt_C<-read.capthist(captfile = "./Nemegt/Nemegt2013_Capture_BC_AS_clubbed.csv", 
                                 trapfile = "./Nemegt/Nemegt2013_Cams_add.csv", 
                                 detector="count", fmt = "trapID", 
                                 trapcovnames = c("Topo",	"Brokenness",	"Grass", "Rgd", "Water"),
                                 binary.usage=FALSE)

boundaryNemegt=readShapeSpatial("./Nemegt//Habitat/Nemegt_StudyArea.shp")
# and plot it
summary(all.data.Nemegt)

summary(all.data.Nemegt_R)

summary(all.data.Nemegt_C)
#traps(all.data.Nemegt)<-addCovariates(traps(all.data.Nemegt), #Add a covariate of waterholes

plot(boundaryNemegt)
plot(x=all.data.Nemegt, add=TRUE)
text(traps(all.data.Nemegt),labels=as.character(1:40),cex=0.75)

# Make mask:
NemegtMask=make.mask(traps(all.data.Nemegt), spacing=500, buffer = 25000, type="trapbuffer", poly=boundaryNemegt)

# Read ruggedness covariate and put it into mask covariate GRIDCODE
SLCost.Nemegt<-readShapePoly("./Nemegt//Habitat/Nemegt_Rgd500m.shp")  #ruggedness pixels averaged over 500m radius
NemegtMask1<-addCovariates(NemegtMask, SLCost.Nemegt)

# Read binary habitat suitability code into mask covariate ...
SLCostBINARY.Nemegt<-readShapePoly("./Tost//Habitat/tost_sl.shp")  #Logistic binary SL habitat created using telemetry data
#plot(SLCostBINARY.Nemegt, add=TRUE)

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

# Standarize Rgd on traps for REDUCED Dataset (this makes fits a bit more stable)
# -----------------------------------------------------------
summary(covariates(traps(all.data.Nemegt_R)))
covariates(traps(all.data.Nemegt_R))$stdRgd = scale(covariates(traps(all.data.Nemegt_R))$Rgd)
summary(covariates(traps(all.data.Nemegt_R)))
head(covariates(traps(all.data.Nemegt_R)))
tail(covariates(traps(all.data.Nemegt_R)))

# Standarize Rgd on traps for CLUBBED Dataset (this makes fits a bit more stable)
# -----------------------------------------------------------
summary(covariates(traps(all.data.Nemegt_C)))
covariates(traps(all.data.Nemegt_C))$stdRgd = scale(covariates(traps(all.data.Nemegt_C))$Rgd)
summary(covariates(traps(all.data.Nemegt_C)))
head(covariates(traps(all.data.Nemegt_C)))
tail(covariates(traps(all.data.Nemegt_C)))

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
summary(covariates(traps(all.data.Nemegt_R)))
summary(covariates(traps(all.data.Nemegt_C)))

plot(NemegtMask1, covariate="stdGC", contour=FALSE, col=terrain.colors(10), legend = FALSE)
head(covariates(traps(all.data.Nemegt)))
head(covariates(NemegtMask1))


Nemegt.cams=traps(all.data.Nemegt)
Nemegt.camsR=traps(all.data.Nemegt_R) #Reduced Dataset captures
Nemegt.camsC=traps(all.data.Nemegt_C) #Clubbed dataset captures

### Now Jump to line 

Nemegt.hhn<-secr.fit(all.data.Nemegt, model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", mask=NemegtMask1)
Nemegt.hhn.detrgd<-secr.fit(all.data.Nemegt, model=list(D~1, lambda0~stdRgd, sigma~stdRgd), detectfn="HHN", mask=NemegtMask1)
Nemegt.hhn.DHab<-secr.fit(all.data.Nemegt, model=list(D~stdGC, lambda0~1, sigma~1), detectfn="HHN", mask=NemegtMask1)
Nemegt.hhn.DHab.detrgd10<-secr.fit(all.data.Nemegt, model=list(D~stdGC, lambda0~stdRgd, sigma~1), detectfn="HHN", mask=NemegtMask1)
Nemegt.hhn.DHab.detrgd01<-secr.fit(all.data.Nemegt, model=list(D~stdGC, lambda0~1, sigma~stdRgd), detectfn="HHN", mask=NemegtMask1)
Nemegt.hhnx<-Nemegt.hhn
Nemegt.hhn.detrgdx<-Nemegt.hhn.detrgd
Nemegt.hhn.DHabx<-Nemegt.hhn.DHab
Nemegt.hhn.DHab.detrgd10x<-Nemegt.hhn.DHab.detrgd10
Nemegt.hhn.DHab.detrgd01x<-Nemegt.hhn.DHab.detrgd01

AIC(Nemegt.hhn, Nemegt.hhn.detrgd,Nemegt.hhn.DHab, Nemegt.hhn.DHab.detrgd10, Nemegt.hhn.DHab.detrgd01)

coefficients(Nemegt.hhn.DHab)
coefficients(Nemegt.hhn.detrgd)
NemegtSurface<-predictDsurface(Nemegt.hhn.DHab, se.D=TRUE, cl.D=TRUE)

windows()
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
Nemegt.hhn.DHab.nonUx<-secr.fit(all.data.Nemegt, detectfn="HHN", mask=NemegtMask1,
                              model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                              details = list(userdist = userdfn1),
                              start = list(noneuc = 1)) #-1 gets rid of the intercept
AIC(Nemegt.hhn,Nemegt.hhn.DHab.nonU)
coefficients(Nemegt.hhn.DHab.nonUx)
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
Nemegt.hhn.D.nonUx<-secr.fit(all.data.Nemegt, detectfn="HHN", mask=NemegtMask1,
                           model=list(D~1, lambda0~1, sigma~1, 
                                      noneuc ~ stdGC -1), 
                           details = list(userdist = userdfn1),
                           start = list(noneuc = 1))

coefficients(Nemegt.hhn.D.nonU)
AIC(Nemegt.hhn,Nemegt.hhn.DHab.nonU, Nemegt.hhn.D.nonU)

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
Nemegt.hhn.DHab.nonU.GBGCx<-secr.fit(all.data.Nemegt, detectfn="HHN", mask=NemegtMask1,
                                   model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC + stdBC-1), 
                                   details = list(userdist = userdfn1),
                                   start = list(noneuc = 1))
coefficients(Nemegt.hhn.DHab.nonU.GBGC)
AIC(Nemegt.hhn,Nemegt.hhn.DHab.nonU, Nemegt.hhn.DHab.nonU.GBGC, Nemegt.hhn.D.nonU)

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
Nemegt.hhn.DHab.nonU.GBx<-secr.fit(all.data.Nemegt, detectfn="HHN", mask=NemegtMask1,
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

# Model with stdGC in noneuc, topography for lambda:
# ---------------------------
Nemegt.hhn.DHab.nonU.LamTopox<-secr.fit(all.data.Nemegt, detectfn="HHN", mask=NemegtMask1,
                               model=list(D~stdGC, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
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

# Model with stdGC in noneuc, Water for lambda:
# ---------------------------
Nemegt.hhn.DHab.nonU.LamWx<-secr.fit(all.data.Nemegt, detectfn="HHN", mask=NemegtMask1,
                                       model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                                       details = list(userdist = userdfn1),
                                       start = list(noneuc = 1)) #-1 gets rid of the intercept

coefficients(Nemegt.hhn.DHab.nonU.LamWx)
# Model with stdGC in noneuc, Topo+Water for lambda:
# ---------------------------
Nemegt.hhn.DHab.nonU.LamTopoWx<-secr.fit(all.data.Nemegt, detectfn="HHN", mask=NemegtMask1,
                                    model=list(D~stdGC, lambda0~Topo+Water, sigma~1, noneuc ~ stdGC -1), 
                                    details = list(userdist = userdfn1),
                                    start = list(noneuc = 1)) #-1 gets rid of the intercept
coefficients(Nemegt.hhn.DHab.nonU.LamTopoW)
coefficients(Nemegt.hhn.D.nonUx)


NemegtAIC=AIC(Nemegt.hhnx, Nemegt.hhn.detrgdx, Nemegt.hhn.DHabx, Nemegt.hhn.DHab.detrgd10x, 
              Nemegt.hhn.DHab.detrgd01x, Nemegt.hhn.DHab.nonUx, Nemegt.hhn.D.nonUx, Nemegt.hhn.DHab.nonU.GBGCx, 
              Nemegt.hhn.DHab.nonU.GBx, Nemegt.hhn.DHab.nonU.LamTopox, Nemegt.hhn.DHab.nonU.LamWx,
              Nemegt.hhn.DHab.nonU.LamTopoWx)
NemegtAIC

coefficients(Nemegt.hhn.DHab.nonU.LamWx)
coefficients(Nemegt.hhn.DHab.nonUx)

write.csv(NemegtAIC, file = "AICNemegtx.csv")

NhatNem.Topmodelx<-region.N(Nemegt.hhn.DHab.nonU.LamWx)
NhatNem.Nullx<-region.N(Nemegt.hhnx)
windows()
NemegtSurfaceX<-predictDsurface(Nemegt.hhn.DHab.nonU.LamWx, se.D=TRUE, cl.D=TRUE)
plot(NemegtSurfaceX,asp=1,contour=FALSE, add=TRUE)
plot(x=all.data.Nemegt, col=lwd())
coefficients(Nemegt.hhn.DHab.nonU.LamWx)

NemegtSurface<-predictDsurface(Nemegt.hhn.DHab.detW, se.D=TRUE, cl.D=TRUE)
plot(NemegtSurface,asp=1,contour=FALSE,col=terrain.colors(40), add = TRUE)
head(covariates(NemegtSurface)$D.0)
head(NemegtSurface)


NemegtSurfaceNU<-predictDsurface(Nemegt.hhn.DHab.nonU.LamW, se.D=TRUE, cl.D=TRUE)
plot(NemegtSurfaceNU,asp=1,contour=FALSE,col=terrain.colors(40))


NemegtSurfaceNU1<-predictDsurface(Nemegt.hhn.DHab.nonU, se.D=TRUE, cl.D=TRUE)
plot(NemegtSurfaceNU1,asp=1,contour=FALSE,col=terrain.colors(40))
# Compare models with and without non-Euclidian distance:
# -----------------------------------------------

# Save objects:
NemegtSurface<-predictDsurface(Nemegt.hhn.DHab.detW, se.D=TRUE, cl.D=TRUE)
NemegtSurface.nonU<-predictDsurface(Nemegt.hhn.DHab.nonU, se.D=TRUE, cl.D=TRUE)
NemegtSurface.D.nonU<-predictDsurface(Nemegt.hhn.D.nonU, se.D=TRUE, cl.D=TRUE)
NemegtSurface.DHab.nonU.GB<-predictDsurface(Nemegt.hhn.DHab.nonU.GB, se.D=TRUE, cl.D=TRUE)
NemegtSurface.DHab.nonU.GBGC<-predictDsurface(Nemegt.hhn.DHab.nonU.GBGC, se.D=TRUE, cl.D=TRUE)
save(Nemegt.cams,NemegtMask1,NemegtSurface,NemegtSurface.nonU,NemegtSurface.D.nonU,Nemegt.hhn.DHab,
     Nemegt.hhn,Nemegt.hhn.detrgd,Nemegt.hhn.D.nonU,Nemegt.hhn.DHab.nonU,Nemegt.hhn.DHab.nonU.GB,
     Nemegt.hhn.DHab.nonU.GBGC,Nemegt.hhn.DHab, Nemegt.hhn.DHab.detTopo10, Nemegt.hhn.DHab.detTopo01, 
     Nemegt.hhn.DHab.detTopoRgd,file="./Nemegt/Nemegt-nonEuc-fits.RData")
# load fitted objects:
load("./Nemegt/Nemegt-nonEuc-fits.RData")

# Compare all model AICs:
AIC(Nemegt.hhn, Nemegt.hhn.detrgd,Nemegt.hhn.DHab,Nemegt.hhn.DHab.nonU,Nemegt.hhn.D.nonU, 
    Nemegt.hhn.DHab.nonU.GBGC, Nemegt.hhn.DHab.nonU.GB)


# Best non-Euclidian
Nemegt.Nhatbest.nonU<-region.N(Nemegt.hhn.D.nonU)
# Best Euclidian
Nemegt.Nhatbest.U<-region.N(Nemegt.hhn)
# Compare them
Nemegt.Nhatbest.nonU
Nemegt.Nhatbest.U

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

########################
#########################
########################

#Running SECR for Nemegt 2013 with Additional Camera Points!!!
#########################
##########################
##########################

# Read capture file and boundary
all.data.Nemegt2x<-read.capthist(captfile = "./Nemegt/Nemegt2013_Capture.csv", 
                               trapfile = "./Nemegt/Nemegt2013_Cams_add.csv", 
                               detector="count", fmt = "trapID", 
                               trapcovnames = c("Topo",	"Brokenness",	"Grass", "Rgd", "Water"))

boundaryNemegt2x=readShapeSpatial("./Nemegt//Habitat/Nemegt_StudyArea.shp")
# and plot it
summary(all.data.Nemegt2x)

#traps(all.data.Nemegt)<-addCovariates(traps(all.data.Nemegt), #Add a covariate of waterholes

plot(boundaryNemegt2x)
plot(x=all.data.Nemegt2x, add=TRUE)
text(traps(all.data.Nemegt2x),labels=as.character(1:40),cex=0.75)

# Make mask:
NemegtMask2x=make.mask(traps(all.data.Nemegt2x), spacing=500, buffer = 25000, type="trapbuffer", poly=boundaryNemegt2x)

# Read ruggedness covariate and put it into mask covariate GRIDCODE
SLCost.Nemegt2x<-readShapePoly("./Nemegt//Habitat/Nemegt_Rgd500m.shp")  #ruggedness pixels averaged over 500m radius
NemegtMask12x<-addCovariates(NemegtMask2x, SLCost.Nemegt2x)

# Read binary habitat suitability code into mask covariate ...
SLCostBINARY.Nemegt2x<-readShapePoly("./Tost//Habitat/tost_sl.shp")  #Logistic binary SL habitat created using telemetry data
#plot(SLCostBINARY.Nemegt, add=TRUE)

NemegtMask12x<-addCovariates(NemegtMask12x, SLCostBINARY.Nemegt2x)
head(covariates(NemegtMask12x))
names(covariates(NemegtMask12x))[3:4] = c("binaryID","BINCODE") #Rename headers
summary(covariates(NemegtMask12x))
# make NAs in BINCODE zeros:
covariates(NemegtMask12x)$BINCODE[is.na(covariates(NemegtMask12x)$BINCODE)] = 0
summary(covariates(NemegtMask12x))


# Standarize Rgd on traps (this makes fits a bit more stable)
# -----------------------------------------------------------
summary(covariates(traps(all.data.Nemegt2x)))
covariates(traps(all.data.Nemegt2x))$stdRgd = scale(covariates(traps(all.data.Nemegt2x))$Rgd)
summary(covariates(traps(all.data.Nemegt2x)))
head(covariates(traps(all.data.Nemegt2x)))

# Standarize GRIDCODE (in stdGC) and BINCODE (in stdBC) on mask
# ------------------------------------------------------------------------
summary(covariates(NemegtMask12x))
covariates(NemegtMask12x)$stdGC = scale(covariates(NemegtMask12x)$GRIDCODE)
covariates(NemegtMask12x)$stdBC = scale(covariates(NemegtMask12x)$BINCODE)
summary(covariates(NemegtMask12x))
names(covariates(NemegtMask12x))


head(covariates(NemegtMask12x))
summary(covariates(NemegtMask12x))
summary(covariates(traps(all.data.Nemegt2x)))

plot(NemegtMask12x, covariate="stdGC", contour=FALSE, col=terrain.colors(10), legend = FALSE)
head(covariates(traps(all.data.Nemegt2x)))
head(covariates(NemegtMask12x))


Nemegt.cams=traps(all.data.Nemegt2x)

Nemegt.hhn2x<-secr.fit(all.data.Nemegt2x, model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", mask=NemegtMask12x)
Nemegt.hhn.detrgd2x<-secr.fit(all.data.Nemegt2x, model=list(D~1, lambda0~stdRgd, sigma~stdRgd), detectfn="HHN", mask=NemegtMask12x)
Nemegt.hhn.DHab2x<-secr.fit(all.data.Nemegt2x, model=list(D~stdGC, lambda0~1, sigma~1), detectfn="HHN", mask=NemegtMask12x)
Nemegt.hhn.DHab.detrgd102x<-secr.fit(all.data.Nemegt2x, model=list(D~stdGC, lambda0~stdRgd, sigma~1), detectfn="HHN", mask=NemegtMask12x)
Nemegt.hhn.DHab.detrgd012x<-secr.fit(all.data.Nemegt2x, model=list(D~stdGC, lambda0~1, sigma~stdRgd), detectfn="HHN", mask=NemegtMask12x)

# Same analysis as above, with reduced Dataset
Nemegt.hhn2xR<-secr.fit(all.data.Nemegt_R, model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", mask=NemegtMask12x)
Nemegt.hhn.detrgd2xR<-secr.fit(all.data.Nemegt_R, model=list(D~1, lambda0~stdRgd, sigma~stdRgd), detectfn="HHN", mask=NemegtMask12x)
Nemegt.hhn.DHab2xR<-secr.fit(all.data.Nemegt_R, model=list(D~stdGC, lambda0~1, sigma~1), detectfn="HHN", mask=NemegtMask12x)
Nemegt.hhn.DHab.detrgd102xR<-secr.fit(all.data.Nemegt_R, model=list(D~stdGC, lambda0~stdRgd, sigma~1), detectfn="HHN", mask=NemegtMask12x)
Nemegt.hhn.DHab.detrgd012xR<-secr.fit(all.data.Nemegt_R, model=list(D~stdGC, lambda0~1, sigma~stdRgd), detectfn="HHN", mask=NemegtMask12x)

# Same analysis as above, with clubbed Dataset
Nemegt.hhn2xC<-secr.fit(all.data.Nemegt_C, model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", mask=NemegtMask12x)
Nemegt.hhn.detrgd2xC<-secr.fit(all.data.Nemegt_C, model=list(D~1, lambda0~stdRgd, sigma~stdRgd), detectfn="HHN", mask=NemegtMask12x)
Nemegt.hhn.DHab2xC<-secr.fit(all.data.Nemegt_C, model=list(D~stdGC, lambda0~1, sigma~1), detectfn="HHN", mask=NemegtMask12x)
Nemegt.hhn.DHab.detrgd102xC<-secr.fit(all.data.Nemegt_C, model=list(D~stdGC, lambda0~stdRgd, sigma~1), detectfn="HHN", mask=NemegtMask12x)
Nemegt.hhn.DHab.detrgd012xC<-secr.fit(all.data.Nemegt_C, model=list(D~stdGC, lambda0~1, sigma~stdRgd), detectfn="HHN", mask=NemegtMask12x)


AIC(Nemegt.hhn2xR, Nemegt.hhn.detrgd2xR,Nemegt.hhn.DHab2xR, Nemegt.hhn.DHab.detrgd102xR, Nemegt.hhn.DHab.detrgd012xR)
AIC(Nemegt.hhn2xR, Nemegt.hhn.DHab.nonU2xR)

coefficients(Nemegt.hhn.DHab2x)
coefficients(Nemegt.hhn.detrgd2x)
NemegtSurface<-predictDsurface(Nemegt.hhn.DHab2x, se.D=TRUE, cl.D=TRUE)

windows()
plot(NemegtSurface,asp=1,contour=FALSE)
plotcovariate(NemegtSurface,covariate="stdGC",asp=1,contour=FALSE)
plot(Nemegt.cams,add=TRUE)

Nhat1<-region.N(Nemegt.hhn.DHab2x) #Estimates the population N of the animals within the region defined by mask
Nhat1

Nhat2<-region.N(Nemegt.hhn2x)
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
Nemegt.hhn.DHab.nonU2x<-secr.fit(all.data.Nemegt2x, detectfn="HHN", mask=NemegtMask12x,
                                model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                                details = list(userdist = userdfn1),
                                start = list(noneuc = 1)) #-1 gets rid of the intercept

Nemegt.hhn.DHab.nonU2xR<-secr.fit(all.data.Nemegt_R, detectfn="HHN", mask=NemegtMask12x,
                                 model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                                 details = list(userdist = userdfn1),
                                 start = list(noneuc = 1)) #-1 gets rid of the intercept

Nemegt.hhn.DHab.nonU2xC<-secr.fit(all.data.Nemegt_C, detectfn="HHN", mask=NemegtMask12x,
                                  model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                                  details = list(userdist = userdfn1),
                                  start = list(noneuc = 1)) #-1 gets rid of the intercept

AIC(Nemegt.hhn,Nemegt.hhn.DHab.nonU2x)
coefficients(Nemegt.hhn.DHab.nonU2x)
NemegtSurface.nonU<-predictDsurface(Nemegt.hhn.DHab.nonU2x, se.D=TRUE, cl.D=TRUE)
windows()
plot(NemegtSurface.nonU,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Nemegt.cams,add=TRUE)
plotcovariate(NemegtSurface.nonU,covariate="stdGC",asp=1,contour=FALSE)
plot(Nemegt.cams,add=TRUE)

Nhat1.nonU<-region.N(Nemegt.hhn.DHab.nonU2x)
Nhat1.nonU

# Model with Flat density and non-Euclidian:
#-------------------------------------------
Nemegt.hhn.D.nonU2x<-secr.fit(all.data.Nemegt2x, detectfn="HHN", mask=NemegtMask12x,
                             model=list(D~1, lambda0~1, sigma~1, 
                                        noneuc ~ stdGC -1), 
                             details = list(userdist = userdfn1),
                             start = list(noneuc = 1))

Nemegt.hhn.D.nonU2xR<-secr.fit(all.data.Nemegt_R, detectfn="HHN", mask=NemegtMask12x,
                              model=list(D~1, lambda0~1, sigma~1, 
                                         noneuc ~ stdGC -1), 
                              details = list(userdist = userdfn1),
                              start = list(noneuc = 1))

Nemegt.hhn.D.nonU2xC<-secr.fit(all.data.Nemegt_C, detectfn="HHN", mask=NemegtMask12x,
                               model=list(D~1, lambda0~1, sigma~1, 
                                          noneuc ~ stdGC -1), 
                               details = list(userdist = userdfn1),
                               start = list(noneuc = 1))

coefficients(Nemegt.hhn.D.nonU2x)
AIC(Nemegt.hhn2x,Nemegt.hhn.DHab.nonU2x, Nemegt.hhn.D.nonU2x)

NemegtSurface.D.nonU2x<-predictDsurface(Nemegt.hhn.D.nonU2x, se.D=TRUE, cl.D=TRUE)
plot(NemegtSurface.D.nonU2x,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Nemegt.cams,add=TRUE)
plotcovariate(NemegtSurface.D.nonU2x,covariate="stdGC",asp=1,contour=FALSE)
plot(Nemegt.cams,add=TRUE)

Nhat1D1.nonU<-region.N(Nemegt.hhn.D.nonU2x)
Nhat1D1.nonU
#AIC(Nemegt.hhn,Nemegt.hhn.D.nonU,Nemegt.hhn.DHab.nonU)

# Model with stdGC and stdBC in noneuc:
# ------------------------------------- NOT TO RUN
#Nemegt.hhn.DHab.nonU.GBGC2x<-secr.fit(all.data.Nemegt2x, detectfn="HHN", mask=NemegtMask12x,
#                                     model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC + stdBC-1), 
#                                     details = list(userdist = userdfn1),
#                                     start = list(noneuc = 1))

#Nemegt.hhn.DHab.nonU.GBGC2xR<-secr.fit(all.data.Nemegt_R, detectfn="HHN", mask=NemegtMask12x,
#                                      model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC + stdBC-1), 
#                                      details = list(userdist = userdfn1),
#                                      start = list(noneuc = 1))

#Nemegt.hhn.DHab.nonU.GBGC2xC<-secr.fit(all.data.Nemegt_C, detectfn="HHN", mask=NemegtMask12x,
#                                       model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC + stdBC-1), 
#                                       details = list(userdist = userdfn1),
#                                       start = list(noneuc = 1))


coefficients(Nemegt.hhn.DHab.nonU.GBGC2x)
AIC(Nemegt.hhn2x,Nemegt.hhn.DHab.nonU2x, Nemegt.hhn.DHab.nonU.GBGC2x, Nemegt.hhn.D.nonU2x)

NemegtSurface.nonU.GBGC2x<-predictDsurface(Nemegt.hhn.DHab.nonU.GBGC2x, se.D=TRUE, cl.D=TRUE)
plot(NemegtSurface.nonU.GBGC2x,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Nemegt.cams,add=TRUE)
plotcovariate(NemegtSurface.nonU.GBGC2x,covariate="stdGC",asp=1,contour=FALSE)
plot(Nemegt.cams,add=TRUE)

Nhat1.nonU<-region.N(Nemegt.hhn.DHab.nonU.GBGC2x)
Nhat1.nonU
AIC(Nemegt.hhn,Nemegt.hhn.D.nonU2x,Nemegt.hhn.DHab.nonU2x, Nemegt.hhn.DHab.nonU.GBGC2x)

# Model with stdBC only in noneuc:
# -------------------------------------
Nemegt.hhn.DHab.nonU.GB2x<-secr.fit(all.data.Nemegt, detectfn="HHN", mask=NemegtMask12x,
                                     model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~stdBC-1), 
                                     details = list(userdist = userdfn1),
                                     start = list(noneuc = 1))

Nemegt.hhn.DHab.nonU.GB2xR<-secr.fit(all.data.Nemegt_R, detectfn="HHN", mask=NemegtMask12x,
                                   model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~stdBC-1), 
                                   details = list(userdist = userdfn1),
                                   start = list(noneuc = 1))

Nemegt.hhn.DHab.nonU.GB2xC<-secr.fit(all.data.Nemegt_C, detectfn="HHN", mask=NemegtMask12x,
                                     model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~stdBC-1), 
                                     details = list(userdist = userdfn1),
                                     start = list(noneuc = 1))


coefficients(Nemegt.hhn.DHab.nonU.GB2x)
NemegtSurface.nonU.GB2x<-predictDsurface(Nemegt.hhn.DHab.nonU.GB2x, se.D=TRUE, cl.D=TRUE)
plot(NemegtSurface.nonU.GB2x,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Nemegt.cams,add=TRUE)
plotcovariate(NemegtSurface.nonU.GB2x,covariate="stdGC",asp=1,contour=FALSE)
plot(Nemegt.cams,add=TRUE)

Nhat1.nonU<-region.N(Nemegt.hhn.DHab.nonU.GB2x)
Nhat1.nonU

# Model with stdGC in noneuc, topography for lambda:
# ---------------------------
Nemegt.hhn.DHab.nonU.LamTopo2x<-secr.fit(all.data.Nemegt2x, detectfn="HHN", mask=NemegtMask12x,
                                        model=list(D~stdGC, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                                        details = list(userdist = userdfn1),
                                        start = list(noneuc = 1)) #-1 gets rid of the intercept

Nemegt.hhn.DHab.nonU.LamTopoR<-secr.fit(all.data.Nemegt_R, detectfn="HHN", mask=NemegtMask12x,
                                         model=list(D~stdGC, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                                         details = list(userdist = userdfn1),
                                         start = list(noneuc = 1)) #-1 gets rid of the intercept

Nemegt.hhn.DHab.nonU.LamTopoC<-secr.fit(all.data.Nemegt_C, detectfn="HHN", mask=NemegtMask12x,
                                        model=list(D~stdGC, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                                        details = list(userdist = userdfn1),
                                        start = list(noneuc = 1)) #-1 gets rid of the intercept



AIC(Nemegt.hhn,Nemegt.hhn.DHab.nonU2x)
coefficients(Nemegt.hhn.DHab.nonU2x)
NemegtSurface.nonU<-predictDsurface(Nemegt.hhn.DHab.nonU2x, se.D=TRUE, cl.D=TRUE)
windows()
plot(NemegtSurface.nonU,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Nemegt.cams,add=TRUE)
plotcovariate(NemegtSurface.nonU,covariate="stdGC",asp=1,contour=FALSE)
plot(Nemegt.cams,add=TRUE)

Nhat1.nonU<-region.N(Nemegt.hhn.DHab.nonU2x)
Nhat1.nonU

# Model with stdGC in noneuc, Water for lambda:
# ---------------------------
Nemegt.hhn.DHab.nonU.LamW2x<-secr.fit(all.data.Nemegt2x, detectfn="HHN", mask=NemegtMask12x,
                                     model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                                     details = list(userdist = userdfn1),
                                     start = list(noneuc = 1)) #-1 gets rid of the intercept

Nemegt.hhn.DHab.nonU.LamW2xR<-secr.fit(all.data.Nemegt_R, detectfn="HHN", mask=NemegtMask12x,
                                      model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                                      details = list(userdist = userdfn1),
                                      start = list(noneuc = 1)) #-1 gets rid of the intercept

Nemegt.hhn.DHab.nonU.GB.LamW2xR<-secr.fit(all.data.Nemegt_R, detectfn="HHN", mask=NemegtMask12x,
                                       model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdBC -1), 
                                       details = list(userdist = userdfn1),
                                       start = list(noneuc = 1)) #-1 gets rid of the intercept

Nemegt.hhn.DHab.nonU.LamW2xC<-secr.fit(all.data.Nemegt_C, detectfn="HHN", mask=NemegtMask12x,
                                       model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                                       details = list(userdist = userdfn1),
                                       start = list(noneuc = 1)) #-1 gets rid of the intercept

coefficients(Nemegt.hhn.DHab.nonU.LamW2x)
# Model with stdGC in noneuc, Topo+Water for lambda:
# ---------------------------
Nemegt.hhn.DHab.nonU.LamTopoW2x<-secr.fit(all.data.Nemegt2x, detectfn="HHN", mask=NemegtMask12x,
                                         model=list(D~stdGC, lambda0~Topo+Water, sigma~1, noneuc ~ stdGC -1), 
                                         details = list(userdist = userdfn1),
                                         start = list(noneuc = 1)) #-1 gets rid of the intercept

Nemegt.hhn.DHab.nonU.LamTopoW2R<-secr.fit(all.data.Nemegt_R, detectfn="HHN", mask=NemegtMask12x,
                                          model=list(D~stdGC, lambda0~Topo+Water, sigma~1, noneuc ~ stdGC -1), 
                                          details = list(userdist = userdfn1),
                                          start = list(noneuc = 1)) #-1 gets rid of the intercept

Nemegt.hhn.DHab.nonU.LamTopoW2C<-secr.fit(all.data.Nemegt_C, detectfn="HHN", mask=NemegtMask12x,
                                          model=list(D~stdGC, lambda0~Topo+Water, sigma~1, noneuc ~ stdGC -1), 
                                          details = list(userdist = userdfn1),
                                          start = list(noneuc = 1)) #-1 gets rid of the intercept

# Warning messages: 1: In autoini(ch, msk, binomN = details$binomN, ignoreusage = details$ignoreusage) :
#  'autoini' failed to find g0; setting initial g0 = 0.1
# 2: In secr.fit(all.data.Nemegt_R, detectfn = "HHN", mask = NemegtMask12x,  :
#                 at least one variance calculation failed 

coefficients(Nemegt.hhn.DHab.nonU.LamTopoW2x)
coefficients(Nemegt.hhn.D.nonU2x)
S1<-predictDsurface(Nemegt.hhn.DHab.nonU.LamTopoW2x)
S2<-predictDsurface(Nemegt.hhn.D.nonU2x)

NemegtAIC2x=AIC(Nemegt.hhn2x, Nemegt.hhn.detrgd2x, Nemegt.hhn.DHab2x, Nemegt.hhn.DHab.detrgd102x, 
                 Nemegt.hhn.DHab.detrgd012x, Nemegt.hhn.DHab.nonU2x, Nemegt.hhn.D.nonU2x,  
                 Nemegt.hhn.DHab.nonU.GB2x, Nemegt.hhn.DHab.nonU.LamTopo2x, Nemegt.hhn.DHab.nonU.LamW2x,
                 Nemegt.hhn.DHab.nonU.LamTopoW2x)
NemegtAIC2x

NemegtAIC2xR=AIC(Nemegt.hhn2xR, Nemegt.hhn.detrgd2xR, Nemegt.hhn.DHab2xR, Nemegt.hhn.DHab.detrgd102xR, 
              Nemegt.hhn.DHab.detrgd012xR, Nemegt.hhn.DHab.nonU2xR, Nemegt.hhn.D.nonU2xR,  
              Nemegt.hhn.DHab.nonU.GB2xR, Nemegt.hhn.DHab.nonU.LamTopoR, Nemegt.hhn.DHab.nonU.LamW2xR,
              Nemegt.hhn.DHab.nonU.LamTopoW2R)
NemegtAIC2xR
coefficients(Nemegt.hhn.DHab.nonU.LamW2xR)
region.N(Nemegt.hhn.DHab.nonU.LamW2xR)

NemegtAIC2xC=AIC(Nemegt.hhn2xC, Nemegt.hhn.detrgd2xC, Nemegt.hhn.DHab2xC, Nemegt.hhn.DHab.detrgd102xC, 
                 Nemegt.hhn.DHab.detrgd012xC, Nemegt.hhn.DHab.nonU2xC, Nemegt.hhn.D.nonU2xC,  
                 Nemegt.hhn.DHab.nonU.GB2xC, Nemegt.hhn.DHab.nonU.LamTopoC, Nemegt.hhn.DHab.nonU.LamW2xC,
                 Nemegt.hhn.DHab.nonU.LamTopoW2C)
NemegtAIC2xC
coefficients(Nemegt.hhn.DHab.nonU.LamW2xR)
coefficients(Nemegt.hhn.DHab.nonU.LamW2xC)
coefficients(Nemegt.hhn.DHab.nonU.LamW2xR)

region.N(Nemegt.hhn.DHab.nonU.LamW2xR)
region.N(Nemegt.hhn.DHab.nonU.LamW2xC)
region.N(Nemegt.hhn.DHab.nonU.LamW2x)

write.csv(NemegtAIC2x, file = "AICNemegt2x.csv")
NemegtSurface2x<-predictDsurface(Nemegt.hhn.DHab.nonU.LamW2x)
plot(NemegtSurface2x, asp=1,contour=FALSE)

coefficients(Nemegt.hhn.DHab.nonU.LamW2xR)
coefficients(Nemegt.hhn.DHab.nonU2xR)
coefficients(Nemegt.hhn.DHab.nonU.GB2xR)

region.N(Nemegt.hhn.DHab.nonU.LamW2xR)
region.N(Nemegt.hhn2x)


region.N(Nemegt.hhn.DHab.nonU.LamW2x)
region.N(Nemegt.hhn2x)

save(Nemegt.hhn2xR, Nemegt.hhn.detrgd2xR, Nemegt.hhn.DHab2xR, Nemegt.hhn.DHab.detrgd102xR,
     Nemegt.hhn.DHab.detrgd012xR, Nemegt.hhn.DHab.nonU2xR, Nemegt.hhn.D.nonU2xR,
     Nemegt.hhn.DHab.nonU.GB2xR, Nemegt.hhn.DHab.nonU.LamTopoR, Nemegt.hhn.DHab.nonU.LamW2xR,
     Nemegt.hhn.DHab.nonU.LamTopoW2R, file="./Nemegt/Nemegt-nonEuc-fit2xR.RData")
load("./Nemegt/Nemegt-nonEuc-fit2xR.RData")

save(NemegtMask12x,Nemegt.hhn2x, Nemegt.hhn.detrgd2x, Nemegt.hhn.DHab2x, Nemegt.hhn.DHab.detrgd102x, 
     Nemegt.hhn.DHab.detrgd012x, Nemegt.hhn.DHab.nonU2x, Nemegt.hhn.D.nonU2x, Nemegt.hhn.DHab.nonU.GBGC2x, 
     Nemegt.hhn.DHab.nonU.GB2x, Nemegt.hhn.DHab.nonU.LamTopo2x, Nemegt.hhn.DHab.nonU.LamW2x,
     Nemegt.hhn.DHab.nonU.LamTopoW2x,file="./Nemegt/Nemegt-nonEuc-fit2x.RData")
# load fitted objects:
load("./Nemegt/Nemegt-nonEuc-fit2x.RData")


