library(secr)
library(fields)
library(maptools)
source("scrplotting.r")
library(gdistance)


#Running SECR for Tost 2012

# Read capture file and boundary
setwd("C:/Users/koust/Dropbox (Snow Leopard Trust)/CREEM/Analyses/snowleopards")
all.data.Tost<-read.capthist(captfile = "./Tost/Tost_capthist2012.csv", binary.usage = FALSE,  
trapfile = "./Tost/Tost_cams_rugged2012.csv", detector="count", 
fmt = "trapID", trapcovnames = c("Rgd", "Topo",	"Altidute",	"Water"))

#all.data.Tost<-read.capthist(captfile = "./Tost/Tost_capthist2012.csv", trapfile = "./Tost/Tost_cams_rugged2012.csv", detector="count", fmt = "trapID", trapcovnames = c("Effort",	"Topo",	"Altidute",	"Rgd", "Water"))
plot(all.data.Tost)
plot(traps(all.data.Tost))
getwd()
boundaryTost=readShapeSpatial("./Tost/Habitat/TostStudy_Area.shp")
# and plot it
plot(boundaryTost, add = TRUE)
plot(x=all.data.Tost, add=TRUE)
summary(all.data.Tost)

# Make mask:
TostMask=make.mask(traps(all.data.Tost), spacing=500, buffer = 25000, type="trapbuffer", poly=boundaryTost)
plot(TostMask)
# Read ruggedness covariate and put it into mask covariate GRIDCODE
SLCost.Tost<-readShapePoly("./Tost//Habitat/Tost_Rgd500m.shp")  #ruggedness pixels averaged over 500m radius
TostMask1<-addCovariates(TostMask, SLCost.Tost)
summary(TostMask1)
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

Tost.hhnx<-secr.fit(all.data.Tost, model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", mask=TostMask1)
Tost.hhn.detTopo10x<-secr.fit(all.data.Tost, model=list(D~1, lambda0~Topo, sigma~1), detectfn="HHN", mask=TostMask1)
Tost.hhn.detWaterx<-secr.fit(all.data.Tost, model=list(D~1, lambda0~Water, sigma~1), detectfn="HHN", mask=TostMask1)
Tost.hhn.DHabx<-secr.fit(all.data.Tost, model=list(D~stdGC, lambda0~1, sigma~1), detectfn="HHN", mask=TostMask1)
Tost.hhn.DHab.LamTopox<-secr.fit(all.data.Tost, model=list(D~stdGC, lambda0~Topo, sigma~1), detectfn="HHN", mask=TostMask1)
Tost.hhn.DHabx.LamWatx<-secr.fit(all.data.Tost, model=list(D~stdGC, lambda0~Water, sigma~1), detectfn="HHN", mask=TostMask1)

AIC(Tost.hhnx, Tost.hhn.detTopo10x, Tost.hhn.detWaterx, Tost.hhn.DHabx, 
    Tost.hhn.DHab.LamTopox, Tost.hhn.DHabx.LamWatx)
Null.Abundance = region.N(Tost.hhnx) 
#AIC(Tost.hhn,Tost.hhn1)
#Tost.hhn.Dx<-secr.fit(all.data.Tost, model=list(D~x, lambda0~1, sigma~1), detectfn="HHN", mask=TostMask1)

coefficients(Tost.hhn.DHabx)
TostSurface<-predictDsurface(Tost.hhn.DHab, se.D=TRUE, cl.D=TRUE)
plot(TostSurface,asp=1,contour=FALSE)
plotcovariate(TostSurface,covariate="stdGC",asp=1,contour=FALSE)
plot(Tost.cams,add=TRUE)

Nhat1<-region.N(Tost.hhn.DHabx) #Estimates the population N of the animals within the region defined by mask
Nhat1

Nhat2<-region.N(Tost.hhn)
Nhat2

# define functions for non-Euclidian distance calculation:
# -------------------------------------------------------

# Create a user-defined tranistion function:
myConductanceFun = function(x) exp(mean(log(x)))

# Create a user-defined function to calculate least-cost distances, from the 
# mask covariate "noneuc" with the user-defined function "myConductanceFun"
mydistFun = function (xy1, xy2, mask) {
  if (missing(xy1)) return("noneuc") # required by secr.fit
  if (!require(gdistance))
    stop ("install package gdistance to use this function")
  # Make raster from mask, because transition() requires raster
  # The mask must contain a covariate called "noneuc", this
  # is added invisibly by secr.fit:
  Sraster <- raster(mask, "noneuc")
  # Calculate the conductances between all points in xy1 and all in xy2, 
  # using a user-defined function called myConductanceFun.
  # (See help for "transition" for more on transitionFunctions.)
  tr <- transition(Sraster, transitionFunction=myConductanceFun, directions=16) 
  # Correction to get more accurate distances:
  # (See help for "geoCorrection" for more on this function.)
  tr <- geoCorrection(tr, type = "c", multpl = FALSE)
  # Pass the object containing the conductances, to costDistance, 
  # which uses their inverse to calculate least-cost distances 
  # by means of Dijkstra's Algorithm.
  costDistance(tr, as.matrix(xy1), as.matrix(xy2))
}



# Non-Euclidian fits
#Clunky, according to David but still works... will be interesting to compare with conductance function above 
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
Tost.hhn.DHab.nonUx<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                             model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                             details = list(userdist = userdfn1),
                             start = list(noneuc = 1)) #-1 gets rid of the intercept

Tost.hhn.DHab.nonUz<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                              model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                              details = list(userdist = mydistFun),
                              start = list(noneuc = 1)) #-1 gets rid of the intercept
Tost.hhn.DHab.a0.nonUz<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                   model=list(D~stdGC, a0~1, sigma~1, noneuc ~ stdGC -1), 
                                   details = list(userdist = mydistFun),
                                   start = list(noneuc = 1))

coefficients(Tost.hhn.DHab.nonUx)
coefficients(Tost.hhn.DHab.nonUz)

AIC(Tost.hhn, Tost.hhn.DHab, Tost.hhn.detTopo10, Tost.hhn.detWater, Tost.hhn.DHab.nonU)

# Model with stdGC in noneuc & Detection Topo:
# ---------------------------
Tost.hhn.DHab.nonU.Topo10x<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                             model=list(D~stdGC, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                             details = list(userdist = userdfn1),
                             start = list(noneuc = 1)) #-1 gets rid of the intercept

Tost.hhn.DHab.nonU.Topo10z<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                     model=list(D~stdGC, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                                     details = list(userdist = mydistFun),
                                     start = list(noneuc = 1)) #-1 gets rid of the intercept
Tost.hhn.DHab.a0Topo.nonU<-secr.fit(all.data.Tost,detectfn="HHN", mask=TostMask1,
                                    model=list(D~stdGC, a0~Topo, sigma~1, noneuc~stdGC -1),
                                    details = list(userdist=mydistFun), start = list(noneuc=1))
Tost.DHab.a0.STopo.nonU<-secr.fit(all.data.Tost,detectfn="HHN", mask=TostMask1,
                                    model=list(D~stdGC, a0~1, sigma~Topo, noneuc~stdGC -1),
                                    details = list(userdist=mydistFun), start = list(noneuc=1))

coefficients(Tost.hhn.DHab.nonU.Topo10x)
coefficients(Tost.hhn.DHab.nonU.Topo10z)

# Model with stdGC in noneuc & Detection Water:
# ---------------------------
Tost.hhn.DHab.nonU.Wx<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                    model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                                    details = list(userdist = userdfn1),
                                    start = list(noneuc = 1)) #-1 gets rid of the intercept

Tost.hhn.DHab.nonU.Wz<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                                details = list(userdist = mydistFun),
                                start = list(noneuc = 1)) #-1 gets rid of the intercept

Tost.hhn.DHab.a0W.nonUz<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                model=list(D~stdGC, a0~Water, sigma~1, noneuc ~ stdGC -1), 
                                details = list(userdist = mydistFun),
                                start = list(noneuc = 1)) #-1 gets rid of the intercept

Tost.hhn.DHab.a0SW.nonUz<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                  model=list(D~stdGC, a0~1, sigma~Water, noneuc ~ stdGC -1), 
                                  details = list(userdist = mydistFun),
                                  start = list(noneuc = 1)) #-1 gets rid of the intercept

coefficients(Tost.hhn.DHab.nonU.Wx)
coefficients(Tost.hhn.DHab.nonU.Wz)

# Model with stdGC in noneuc & Detection Topo water:
# ---------------------------
Tost.hhn.DHab.nonU.T01Wx<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                               model=list(D~stdGC, lambda0~Topo+Water, sigma~1, noneuc ~ stdGC -1), 
                               details = list(userdist = userdfn1),
                               start = list(noneuc = 1)) #-1 gets rid of the intercept

Tost.hhn.DHab.nonU.T01Wz<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                   model=list(D~stdGC, lambda0~Topo+Water, sigma~1, noneuc ~ stdGC -1), 
                                   details = list(userdist = mydistFun),
                                   start = list(noneuc = 1)) #-1 gets rid of the intercept

Tost.hhn.DHab.a0TW.nonUz<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                   model=list(D~stdGC, a0~Topo+Water, sigma~1, noneuc ~ stdGC -1), 
                                   details = list(userdist = mydistFun),
                                   start = list(noneuc = 1))

Tost.hhn.DHab.a0ST_W.nonUz<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                   model=list(D~stdGC, a0~1, sigma~Topo+Water, noneuc ~ stdGC -1), 
                                   details = list(userdist = mydistFun),
                                   start = list(noneuc = 1))


TostSurface.nonU<-predictDsurface(Tost.hhn.DHab.nonU, se.D=TRUE, cl.D=TRUE)
plot(TostSurface.nonU,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Tost.cams,add=TRUE)
plotcovariate(TostSurface.nonU,covariate="stdGC",asp=1,contour=FALSE)
plot(Tost.cams,add=TRUE)

Nhat1.nonU<-region.N(Tost.hhn.DHab.nonU)
Nhat1.nonU

# Try with flat density and non-Euclidian:
Tost.hhn.D.nonUx<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                             model=list(D~1, lambda0~1, sigma~1, 
                                        noneuc ~ stdGC -1), 
                             details = list(userdist = userdfn1),
                             start = list(noneuc = 1))

Tost.hhn.D.nonUz<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                           model=list(D~1, lambda0~1, sigma~1, 
                                      noneuc ~ stdGC -1), 
                           details = list(userdist = mydistFun),
                           start = list(noneuc = 1))

Tost.hhn.D.a0.nonUz<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                           model=list(D~1, a0~1, sigma~1, 
                                      noneuc ~ stdGC -1), 
                           details = list(userdist = mydistFun),
                           start = list(noneuc = 1))

# Try with flat density, water and non-Euclidian:
Tost.hhn.D.DetW.nonUx<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                           model=list(D~1, lambda0~Water, sigma~1, 
                                      noneuc ~ stdGC -1), 
                           details = list(userdist = userdfn1),
                           start = list(noneuc = 1))

Tost.hhn.D.DetW.nonUz<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                model=list(D~1, lambda0~Water, sigma~1, 
                                           noneuc ~ stdGC -1), 
                                details = list(userdist = mydistFun),
                                start = list(noneuc = 1))

Tost.hhn.a0W.nonUz<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                model=list(D~1, a0~Water, sigma~1, 
                                           noneuc ~ stdGC -1), 
                                details = list(userdist = mydistFun),
                                start = list(noneuc = 1))

Tost.hhn.a0SW.nonUz<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                            model=list(D~1, a0~1, sigma~Water, 
                                       noneuc ~ stdGC -1), 
                            details = list(userdist = mydistFun),
                            start = list(noneuc = 1))


# Try with flat density, Topo and non-Euclidian:
Tost.hhn.D.DetT10.nonUx<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                model=list(D~1, lambda0~Topo, sigma~1, 
                                           noneuc ~ stdGC -1), 
                                details = list(userdist = userdfn1),
                                start = list(noneuc = 1))

Tost.hhn.D.DetT10.nonUz<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                  model=list(D~1, lambda0~Topo, sigma~1, 
                                             noneuc ~ stdGC -1), 
                                  details = list(userdist = mydistFun),
                                  start = list(noneuc = 1))

Tost.hhn.D.a0Topo.nonUz<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                  model=list(D~1, a0~Topo, sigma~1, noneuc~stdGC -1),
                                  details = list(userdist = mydistFun), start = list(noneuc=1))

coefficients(Tost.hhn.D.a0Topo.nonUz)
coefficients(Tost.hhn.D.nonU)
TostSurface.D.nonU<-predictDsurface(Tost.hhn.D.nonU, se.D=TRUE, cl.D=TRUE)
plot(TostSurface.D.nonU,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Tost.cams,add=TRUE)
plotcovariate(TostSurface.D.nonU,covariate="stdGC",asp=1,contour=FALSE)
plot(Tost.cams,add=TRUE)

Nhat1D1.nonU<-region.N(Tost.hhn.D.nonU)
Nhat1D1.nonU

# Model with stdGC in D and stdBC in noneuc:
# -------------------------------------
Tost.hhn.DHab.nonU.GBx<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                 model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdBC-1), 
                                 details = list(userdist = userdfn1),
                                 start = list(noneuc = 1),trace=0)

Tost.hhn.DHab.nonU.GBy<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                 model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdBC-1), 
                                 details = list(userdist = userdfn2),
                                 start = list(noneuc = 1),trace=0)

Tost.hhn.DHab.nonU.GBz<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                 model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdBC-1), 
                                 details = list(userdist = mydistFun),
                                 start = list(noneuc = 1),trace=0)

coefficients(Tost.hhn.DHab.nonU.GBz)
coefficients(Tost.hhn.DHab.nonU.GBy)
coefficients(Tost.hhn.DHab.nonU.GBx)

# Look at parameter correlations (look at that between sigma and noneuc.stdBC)
cov2cor(Tost.hhn.DHab.nonU.GBz$beta.vcv)
cov2cor(Tost.hhn.DHab.nonU.GBy$beta.vcv)
cov2cor(Tost.hhn.DHab.nonU.GBx$beta.vcv)

# AICs:
AIC(Tost.hhn.DHab.nonU.GBx,Tost.hhn.DHab.nonU.GBy,Tost.hhn.DHab.nonU.GBz)

# Abundances
Nest.x = region.N(Tost.hhn.DHab.nonU.GBx)
Nest.y = region.N(Tost.hhn.DHab.nonU.GBy)
Nest.z = region.N(Tost.hhn.DHab.nonU.GBz)
Nest.x;Nest.y;Nest.z

# Double-check abundances manually
cell.area = 25
Dbeta0.x = coefficients(Tost.hhn.DHab.nonU.GBx)["D","beta"]
DbetastdGC.x = coefficients(Tost.hhn.DHab.nonU.GBx)["D.stdGC","beta"]
Nhat.x = sum(exp(Dbeta0.x + DbetastdGC.x*covariates(TostMask1)$stdGC))*cell.area
Dbeta0.y = coefficients(Tost.hhn.DHab.nonU.GBy)["D","beta"]
DbetastdGC.y = coefficients(Tost.hhn.DHab.nonU.GBy)["D.stdGC","beta"]
Nhat.y = sum(exp(Dbeta0.y + DbetastdGC.y*covariates(TostMask1)$stdGC))*cell.area
Dbeta0.z = coefficients(Tost.hhn.DHab.nonU.GBz)["D","beta"]
DbetastdGC.z = coefficients(Tost.hhn.DHab.nonU.GBz)["D.stdGC","beta"]
Nhat.z = sum(exp(Dbeta0.z + DbetastdGC.z*covariates(TostMask1)$stdGC))*cell.area
Nhat.x;Nhat.y;Nhat.z

# Look at conductances in space (at a point - not quite what is used)
beta.x = coefficients(Tost.hhn.DHab.nonU.GBx)["noneuc.stdBC","beta"]
noneuc.x = exp(beta.x * covariates(TostMask1)$stdBC)
beta.y = coefficients(Tost.hhn.DHab.nonU.GBy)["noneuc.stdBC","beta"]
noneuc.y = exp(beta.y * covariates(TostMask1)$stdBC)
beta.z = coefficients(Tost.hhn.DHab.nonU.GBz)["noneuc.stdBC","beta"]
noneuc.z = exp(beta.z * covariates(TostMask1)$stdBC)
predmask = TostMask1
covariates(predmask)$conductance.x = 1/noneuc.x
covariates(predmask)$conductance.y = noneuc.x
covariates(predmask)$conductance.z = noneuc.z # this is only approximate (`cause exponent outside mean)
quartz(h=9,w=4)
par(mfrow=c(3,1))
plotcovariate(predmask,covariate="conductance.x",contour=FALSE,col=c("green","red"))
plotcovariate(predmask,covariate="conductance.y",contour=FALSE,col=c("green","red"))
plotcovariate(predmask,covariate="conductance.z",contour=FALSE,col=c("green","red"))
quartz(h=9,w=4)
par(mfrow=c(3,1))
hist(covariates(predmask)$conductance.x)
hist(covariates(predmask)$conductance.y)
hist(covariates(predmask)$conductance.z)

# Plot density surfaces
Dsurf.x = predictDsurface(Tost.hhn.DHab.nonU.GBx, se.D=TRUE, cl.D=TRUE)
Dsurf.y = predictDsurface(Tost.hhn.DHab.nonU.GBy, se.D=TRUE, cl.D=TRUE)
Dsurf.z = predictDsurface(Tost.hhn.DHab.nonU.GBz, se.D=TRUE, cl.D=TRUE)
plotcovariate(Dsurf.x,covariate="D.0",contour=FALSE,asp=1)
plot(Tost.cams,add=TRUE,detpar=list(pch=19,col="white"));plot(Tost.cams,add=TRUE)
plotcovariate(Dsurf.y,covariate="D.0",contour=FALSE,asp=1)
plot(Tost.cams,add=TRUE,detpar=list(pch=19,col="white"));plot(Tost.cams,add=TRUE)
plotcovariate(Dsurf.z,covariate="D.0",contour=FALSE,asp=1)
plot(Tost.cams,add=TRUE,detpar=list(pch=19,col="white"));plot(Tost.cams,add=TRUE)

plotcovariate(TostMask1,covariate="stdBC",contour=FALSE,col=c("green","red"),asp=1)
plot(Tost.cams,add=TRUE,detpar=list(pch=19,col="white"));plot(Tost.cams,add=TRUE)


# Plot where the density accounting for the top X% of abudance is
pc = 25 # Top percentage of abundance to plot
plotNpc(Dsurf.x,pc,contour=TRUE,boundary=boundaryTost)
plotNpc(Dsurf.y,pc,contour=TRUE,boundary=boundaryTost)
plotNpc(Dsurf.z,pc,contour=TRUE,boundary=boundaryTost)



TostSurface.nonU<-predictDsurface(Tost.hhn.DHab.nonU.GBGC, se.D=TRUE, cl.D=TRUE)
plot(TostSurface.nonU.GBGC,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Tost.cams,add=TRUE)
plotcovariate(TostSurface.nonU.GBGC,covariate="stdGC",asp=1,contour=FALSE)
plot(Tost.cams,add=TRUE)

Nhat1.nonU.GBGC<-region.N(Tost.hhn.DHab.nonU.GBGC)
Nhat1.nonU.GBGC

# Model with stdGC in noneuc and stdBC in D:
# -------------------------------------
Tost.hhn.DGB.nonU.GCx<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                model=list(D~stdBC, lambda0~1, sigma~1, noneuc ~stdGC-1), 
                                details = list(userdist = userdfn1),
                                start = list(noneuc = 1))

Tost.hhn.DGB.nonU.GCz<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                 model=list(D~stdBC, lambda0~1, sigma~1, noneuc ~stdGC-1), 
                                 details = list(userdist = mydistFun),
                                 start = list(noneuc = 1))

# Model with stdGC in nonEuc and Sigma affected by Topography
Tost.hhn.SigTopo.DGC.nonU.GCz<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                model=list(D~stdGC, lambda0~1, sigma~Topo, noneuc ~stdGC-1), 
                                details = list(userdist = mydistFun),
                                start = list(noneuc = 1))

Tost.a0.SigTopo.DGC.nonU.GCz<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                        model=list(D~stdGC, a0~1, sigma~Topo, noneuc ~stdGC-1), 
                                        details = list(userdist = mydistFun),
                                        start = list(noneuc = 1))

Tost.a0Topo.DGC.nonU.GCz<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
                                       model=list(D~stdGC, a0~Topo, sigma~1, noneuc ~stdGC-1), 
                                       details = list(userdist = mydistFun),
                                       start = list(noneuc = 1))

AIC(Tost.a0Topo.DGC.nonU.GCz, Tost.a0.SigTopo.DGC.nonU.GCz, Tost.hhn.SigTopo.DGC.nonU.GCz)
coefficients(Tost.a0.SigTopo.DGC.nonU.GCz)

coefficients(Tost.hhn.DGB.nonU.GCx)
coefficients(Tost.hhn.DGB.nonU.GCz)
TostSurface.nonU<-predictDsurface(Tost.hhn.DHab.nonU.GBGCx, se.D=TRUE, cl.D=TRUE)
plot(TostSurface.nonU,asp=1,contour=FALSE,col=terrain.colors(40))
plot(Tost.cams,add=TRUE)
plotcovariate(TostSurface.nonU.GB,covariate="stdGC",asp=1,contour=FALSE)
plot(Tost.cams,add=TRUE)

Nhat1.nonU.GB<-region.N(Tost.hhn.DHab.nonU.GBGCx)
Nhat1.nonU.GB
region.N(Tost.hhnx)

#Ignoring the models with binary habitat covariate GB for ocmparison
AICTost=AIC(Tost.hhnx, Tost.hhn.detTopo10x, Tost.hhn.detWaterx, Tost.hhn.DHabx, Tost.hhn.DHab.nonUx, 
            Tost.hhn.DHab.nonU.Topo10x,Tost.hhn.DHab.nonU.Wx, Tost.hhn.DHab.nonU.T01Wx, 
            Tost.hhn.D.nonUx, Tost.hhn.DHab.nonU.GBx, Tost.hhn.DHab.nonU.GBGCx,
            Tost.hhn.DHab.LamTopox, Tost.hhn.DHabx.LamWatx)

AICTostz=AIC(Tost.hhnx, Tost.hhn.detTopo10x, Tost.hhn.detWaterx, Tost.hhn.DHabx, Tost.hhn.DHab.nonUz, 
              Tost.hhn.DHab.nonU.Topo10z, Tost.hhn.DHab.nonU.Wz, Tost.hhn.DHab.nonU.T01Wz, Tost.hhn.D.nonUz, 
              Tost.hhn.D.DetW.nonUz, Tost.hhn.D.DetT10.nonUz, Tost.hhn.SigTopo.DGC.nonU.GCz,
             Tost.a0Topo.DGC.nonU.GCz, Tost.a0.SigTopo.DGC.nonU.GCz, Tost.hhn.D.a0Topo.nonUz,
             Tost.hhn.DHab.a0.nonUz, criterion = "AICc")

AIC(Tost.hhnx, Tost.hhn.detTopo10x, Tost.hhn.detWaterx, Tost.hhn.DHabx, criterion = "AICc")
AICTostz
AICTost
AICTostcoefficients(Tost.hhn.DHab.nonU.GBx)
coefficients(Tost.hhn.DHab.nonU.Topo10z)

save(Tost.hhnx, Tost.hhn.detTopo10x, Tost.hhn.detWaterx, Tost.hhn.DHabx, Tost.hhn.DHab.nonUz, 
  Tost.hhn.DHab.nonU.Topo10z, Tost.hhn.DHab.nonU.Wz, Tost.hhn.DHab.nonU.T01Wz, Tost.hhn.D.nonUz, 
  Tost.hhn.D.DetW.nonUz, Tost.hhn.D.DetT10.nonUz, Tost.hhn.SigTopo.DGC.nonU.GCz,
  Tost.a0Topo.DGC.nonU.GCz, Tost.a0.SigTopo.DGC.nonU.GCz, Tost.hhn.D.a0Topo.nonUz,
  Tost.hhn.DHab.a0.nonUz, file="./Tost/Tost-nonEuc-fitsz.RData")

AICTost

coefficients(Tost.hhn.D.DetT10.nonUz)
coefficients(Tost.a0.SigTopo.DGC.nonU.GCz)
region.N(Tost.hhn.D.DetT10.nonUz)
region.N(Tost.hhn.DHab.nonU.GBz)

write.csv(AICTost, file = "AICTostx.csv")
AICTostx = read.csv("AICTostx.csv")
AICTostx

save(Tost.hhnx, Tost.hhn.detTopo10x, Tost.hhn.detWaterx, Tost.hhn.DHabx, Tost.hhn.DHab.nonUx, 
      Tost.hhn.DHab.nonU.Topo10x,Tost.hhn.DHab.nonU.Wx, Tost.hhn.DHab.nonU.T01Wx, Tost.hhn.D.nonUx, 
      Tost.hhn.DHab.nonU.GBGCx, Tost.hhn.DHab.nonU.GBx, Tost.hhn.DHab.LamTopox, Tost.hhn.DHabx.LamWatx,
     file="./Tost/Tost-nonEuc-fitsx.RData")
load("./Tost/Tost-nonEuc-fitsx.RData")

Nhat1.nonU.Topo10<-region.N(Tost.hhn.DHab.nonU.Topo10)
Nhat1null<-region.N(Tost.hhn)

TopSurface=predictDsurface(Tost.hhn.DHab.nonUx)
plot(TopSurface,asp=1,contour=FALSE,col=terrain.colors(40))
coefficients(Tost.hhn.DHab.nonUx)

Nhat.TopModelx<-region.N(Tost.hhn.DHab.nonUx)
Nhat.null<-region.N(Tost.hhnx)

TostSurfaceX<-predictDsurface(Tost.hhn.DHab.nonU.Topo10x, se.D=TRUE, cl.D=TRUE)
plot(TostSurfaceX,asp=1,contour=FALSE)


# Compare with and without non-Euclidian distance:
# -----------------------------------------------

# First save objects so don't have to refit:
TostSurface<-predictDsurface(Tost.hhn.DHab, se.D=TRUE, cl.D=TRUE)
TostSurface.nonU<-predictDsurface(Tost.hhn.DHab.nonU, se.D=TRUE, cl.D=TRUE)
TostSurface.D.nonU<-predictDsurface(Tost.hhn.D.nonU, se.D=TRUE, cl.D=TRUE)
TostSurface.DHab.nonU.GB<-predictDsurface(Tost.hhn.DHab.nonU.GB, se.D=TRUE, cl.D=TRUE)
TostSurface.DHab.nonU.GBGC<-predictDsurface(Tost.hhn.DHab.nonU.GBGC, se.D=TRUE, cl.D=TRUE)
save(Tost.cams,TostMask1,TostSurface,TostSurface.nonU,TostSurface.D.nonU,Tost.hhn.D.nonU,
     Tost.hhn.DHab.nonU,Tost.hhn.DHab.nonU.GB,Tost.hhn.DHab.nonU.GBGC,Tost.hhn.DHab,
     file="./Tost/Tost-nonEuc-fits2.RData")
# load fitted objects:
load("./Tost/Tost-nonEuc-fits2.RData") #first round analysis

save(Tost.cams,TostMask1,TostSurface,TostSurface.nonU,TostSurface.D.nonU,Tost.hhn.DHab, Tost.hhn.detWater, 
     Tost.hhn.detTopo10, Tost.hhn, Tost.hhn.DHab.nonU.W, Tost.hhn.DHab.nonU.T01W, Tost.hhn.DHab.nonU.Topo10, 
     Tost.hhn.DHab.nonU, file="./Tost/Tost-nonEuc-fits3.RData")
load("./Tost/Tost-nonEuc-fits3.RData") #second round analysis, using nonU, water, topo as det covs

esa(derived(Tost.hhn.DHab.nonU))
esa.plot(Tost.hhn.DHab, xlim = c(0,500), ylim = c(0,10))
abline(v = 4 * 25.6, col = "red", lty = 2)
(region.N(Tost.hhn.DHab.nonU)*100/211650)
DRange<-predictDsurface(Tost.hhn.DHab, cl.D=TRUE, alpha=0.05)
covariates(DRange)[which(covariates(DRange)$D.0==(max(covariates(DRange)$D.0))),]
covariates(DRange)[which(covariates(DRange)$D.0==(min(covariates(DRange)$D.0))),]

max(covariates(DRange)$lcl)
Dlim=range(covariates(DRange)$D.0)



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

