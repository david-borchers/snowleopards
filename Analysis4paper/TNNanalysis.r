library(secr)
library(fields)
library(maptools)
source("scrplotting.r")

#----------------------- Load the RData objects made by Make_TNN_RData.r --------------------------
load("./Analysis4paper/TNN_boundaries.RData") # Tostboundary,Noyonboundary,Nemegtboundary
load("./Analysis4paper/TNN_masks.RData") # TostMask,NoyonMask,NemegtMask
load("./Analysis4paper/TNN_caphists.RData") # Tost_ch,Noyon_ch,Nemegt_chTNN_ch
#----------------------- ----------------------------------------------- --------------------------

#----------------------- Do some plots to check data seems OK --------------------------

# Find extent of bounding boxes of boundarys:
bbox.Nemegt = bbox(Nemegtboundary)
bbox.Noyon = bbox(Noyonboundary)
bbox.Tost = bbox(Tostboundary)
bbxlim = range(bbox.Nemegt["x",],bbox.Noyon["x",],bbox.Tost["x",]) #Set plot limit for all 3 areas together
bbylim = range(bbox.Nemegt["y",],bbox.Noyon["y",],bbox.Tost["y",])

dxlim = diff(bbxlim)
xlim = c(bbxlim[1],bbxlim[2]+0.15*dxlim)
ylim = bbylim

quartz(w=8,h=4)

# Plot boundaries
plot(xlim,ylim,xlim=xlim,ylim=ylim,type="n",asp=1,bty="n",xlab="Easting",ylab="Northing") 
plot(Tostboundary,add=TRUE)
plot(Noyonboundary,add=TRUE)
plot(Nemegtboundary,add=TRUE)

# Plot GC means
zlim = range(covariates(NemegtMask)$rmeanGC,
             covariates(NoyonMask)$rmeanGC,
             covariates(TostMask)$rmeanGC)
plot(xlim,ylim,xlim=xlim,ylim=ylim,type="n",asp=1,bty="n",xlab="Easting",ylab="Northing") 
plotcovariate(NemegtMask,covariate="rmeanGC",add=TRUE,zlim=zlim,contour=FALSE)
plotcovariate(NoyonMask,covariate="rmeanGC",add=TRUE,zlim=zlim,contour=FALSE)
plotcovariate(TostMask,covariate="rmeanGC",add=TRUE,zlim=zlim,contour=FALSE)

# Plot deviations from GC mean
zlim = range(covariates(NemegtMask)$rmeanGCdev,
             covariates(NoyonMask)$rmeanGCdev,
             covariates(TostMask)$rmeanGCdev)
plot(xlim,ylim,xlim=xlim,ylim=ylim,type="n",asp=1,bty="n",xlab="Easting",ylab="Northing") 
plotcovariate(NemegtMask,covariate="rmeanGCdev",add=TRUE,zlim=zlim,contour=FALSE)
plotcovariate(NoyonMask,covariate="rmeanGCdev",add=TRUE,zlim=zlim,contour=FALSE)
plotcovariate(TostMask,covariate="rmeanGCdev",add=TRUE,zlim=zlim,contour=FALSE)


zlim = range(covariates(NemegtMask)$GC,
             covariates(NoyonMask)$GC,
             covariates(TostMask)$GC)
#pdf("Allregions.pdf",h=5,w=10)
plot(bbxlim,bbylim,xlim=xlim,ylim=ylim,xlab="",ylab="",bty="n",type="n",xaxt="n",yaxt="n",asp=1) 
# Plot the terrain
plotcovariate(NoyonMask, covariate="GC", contour = FALSE, col = terrain.colors(16), zlim=zlim, add = TRUE)
plotcovariate(NemegtMask, covariate="GC", contour=FALSE, col=terrain.colors(16), zlim=zlim, add = TRUE)
plotcovariate(TostMask, covariate="GC", contour=FALSE, col=terrain.colors(16), zlim=zlim, add=TRUE)
# Add the traps
plot(traps(all.data.TNN)[[1]],add=TRUE,detpar=list(col="black",pch="+"))
plot(traps(all.data.TNN)[[2]],add=TRUE,detpar=list(col="black",pch="+"))
plot(traps(all.data.TNN)[[3]],add=TRUE,detpar=list(col="black",pch="+"))
# and the borders over the top
plot(Nemegtboundary,add=TRUE)#,border=3)
plot(Noyonboundary,add=TRUE)#,border=2)
plot(Tostboundary,add=TRUE)#,border=1)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
#dev.off()

# Plot Nemegt data with ADDITOINAL traps added:
plotcovariate(NemegtMask, covariate="GC", contour = FALSE, col = terrain.colors(16), zlim=zlim,asp=1)
NemegtADDITIONAL_ch<-read.capthist(captfile = "./Analysis4paper/Data/Nemegt_capthist.csv", 
                            trapfile = "./Analysis4paper/Data/Nemegt_ADDITIONALcams.csv", 
                            detector="count", binary.usage=FALSE, fmt = "trapID", 
                            trapcovnames = c("Rgd", "Topo", "Water", "Winter"))
plot(traps(NemegtADDITIONAL_ch),add=TRUE)
plot(NemegtADDITIONAL_ch,tracks=TRUE,add=TRUE)


# -------------------------- End of Plotting ---------------------------------



# ----------------------- Base-case simple models for individual regions ---------------------------
Nemegt.basic<-secr.fit(Nemegt_ch, detectfn="HHN", mask=NemegtMask,
                         model=list(D~1, lambda0~Water, sigma~1))
Noyon.basic<-secr.fit(Noyon_ch, detectfn="HHN", mask=NoyonMask,
                     model=list(D~1, lambda0~1, sigma~1))
Tost.basic<-secr.fit(Tost_ch, detectfn="HHN", mask=TostMask,
                   model=list(D~1, lambda0~1, sigma~1))



# ----------------------- Individual region models ---------------------------
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

Nemegt.NUstdGC<-secr.fit(Nemegt_ch, detectfn="HHN", mask=NemegtMask,
                         model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                         details = list(userdist = userdfn1),
                         start = list(noneuc = 1))
# Equivalent model but with rmeanGCdev (should give same density estimates):
Nemegt.NUrmeanGCdev<-secr.fit(Nemegt_ch, detectfn="HHN", mask=NemegtMask,
                         model=list(D~rmeanGCdev, lambda0~Water, sigma~1, noneuc ~ rmeanGCdev -1), 
                         details = list(userdist = userdfn1),
                         start = list(noneuc = 1))
# Do again, starting where left off:
Nemegt.NUrmeanGCdev<-secr.fit(Nemegt_ch, detectfn="HHN", mask=NemegtMask,
                              model=list(D~rmeanGCdev, lambda0~Water, sigma~1, noneuc ~ rmeanGCdev -1), 
                              details = list(userdist = userdfn1),
                              start = Nemegt.NUrmeanGCdev)
coefficients(Nemegt.NUstdGC)
coefficients(Nemegt.NUrmeanGCdev)
Nest.Nemegt.NUstdGC = region.N(Nemegt.NUstdGC)
Nest.Nemegt.NUrmeanGCdev = region.N(Nemegt.NUrmeanGCdev)
Nest.Nemegt.NUstdGC;Nest.Nemegt.NUrmeanGCdev # Slightly different SEs & CIs: must be lack of asymptotics
AIC(Nemegt.NUstdGC,Nemegt.NUrmeanGCdev)

# replace noneuc(stdGC) with noneuc(BC)
Nemegt.NUBC<-secr.fit(Nemegt_ch, detectfn="HHN", mask=NemegtMask,
                      model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ BC -1), 
                      details = list(userdist = userdfn1),
                      start = list(noneuc = 1))
AIC(Nemegt.NUBC,Nemegt.NUstdGC)

# Try without lambda0~Water
Nemegt.NUstdGC.W<-secr.fit(Nemegt_ch, detectfn="HHN", mask=NemegtMask,
                         model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                         details = list(userdist = userdfn1),
                         start = list(noneuc = 1))
AIC(Nemegt.NUBC,Nemegt.NUstdGC,Nemegt.NUstdGC.W)


# Try Nemegt with best model but with ADDITIONAL traps in Steppe:
NemegtADD.NUstdGC<-secr.fit(NemegtADDITIONAL_ch, detectfn="HHN", mask=NemegtMask,
                         model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                         details = list(userdist = userdfn1),
                         start = list(noneuc = 1))



# Noyon: D(GC); noneuc(BC):
# ---------------------------
Noyon.NUBC<-secr.fit(Noyon_ch, detectfn="HHN", mask=NoyonMask,
                     model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ BC -1), 
                     details = list(userdist = userdfn1),
                     start = list(noneuc = 1))
# Cannot get Tost.NUBC to converge!
Noyon.NUBC<-secr.fit(Noyon_ch, detectfn="HHN", mask=NoyonMask,
                     model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ BC -1), 
                     details = list(userdist = userdfn1),
                     start = Noyon.NUBC)
# replace noneuc(BC) with noneuc(stdGC)
Noyon.NUstdGC<-secr.fit(Noyon_ch, detectfn="HHN", mask=NoyonMask,
                        model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                        details = list(userdist = userdfn1),
                        start = list(noneuc = 1))
Noyon.NUstdGC<-secr.fit(Noyon_ch, detectfn="HHN", mask=NoyonMask,
                        model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                        details = list(userdist = userdfn1),
                        start = Noyon.NUstdGC)
# With water in lambda0
Noyon.NUstdGC.W<-secr.fit(Noyon_ch, detectfn="HHN", mask=NoyonMask,
                        model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                        details = list(userdist = userdfn1),
                        start = Noyon.NUstdGC)
AIC(Noyon.NUstdGC,Noyon.NUstdGC.W)

# Tost: D(GC); noneuc(BC):
# ---------------------------
Tost.NUBC<-secr.fit(Tost_ch, detectfn="HHN", mask=TostMask,
                    model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ BC -1), 
                    details = list(userdist = userdfn1),
                    start = list(noneuc = 1))
Tost.NUBC<-secr.fit(Tost_ch, detectfn="HHN", mask=TostMask,
                    model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ BC -1), 
                    details = list(userdist = userdfn1),
                    start = Tost.NUBC)
# replace noneuc(BC) with noneuc(stdGC)
Tost.NUstdGC<-secr.fit(Tost_ch, detectfn="HHN", mask=TostMask,
                       model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                       details = list(userdist = userdfn1),
                       start = list(noneuc = 1))
Tost.NUstdGC.W<-secr.fit(Tost_ch, detectfn="HHN", mask=TostMask,
                       model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                       details = list(userdist = userdfn1),
                       start = list(noneuc = 1))
AIC(Tost.NUBC,Tost.NUstdGC,Tost.NUstdGC.W)

# Conclusion so far: model: D~stdGC in all regions
#                           lambda0~1 for Tost & Noyon; lambda0~Water for Nemegt
#                           sigma~1 # for all regions
#                           noneuc ~ stdGC -1 # for  all regions


# ----------------------- Combined region models ---------------------------
# Fit (almost) the same model simultaneously to all, but with log-linear relationship between the 
# density model intercepts across the regions:
TNN.GCmean_dev.WW <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                           model=list(D~rmeanGC+rmeanGCdev, lambda0~Water:Winter, sigma~1, noneuc ~ stdGC -1), 
                           details = list(userdist = userdfn1),
                           start = list(noneuc = 1))
region.aics = AIC(Nemegt.NUstdGC,Noyon.NUstdGC,Tost.NUstdGC)
sum(region.aics$AIC)
AIC(TNN.GCmean_dev.WW)$AIC
AIC(TNN.GCmean_dev.WW)$AIC - sum(region.aics$AIC)
coefficients(Nemegt.NUstdGC)
coefficients(Noyon.NUstdGC)
coefficients(Tost.NUstdGC)
coefficients(TNN.GCmean_dev.WW)

# Conclusion so far:
# 1. Nemegt has stdGC REDUCING density and INCREASING conductance; also lambda0~Water
# 2. Tost and Noyon have the opposite effects, and no lambda0 effect
# 3. Combined model, which requires stdGC effects to be the same in all regions, is worse by 39 AIC units
# 
save(Tost.NUBC,Tost.NUstdGC,Tost.NUstdGC.W,
     Noyon.NUstdGC,Noyon.NUstdGC.W,
     Nemegt.NUBC,Nemegt.NUstdGC,Nemegt.NUstdGC.W,
     TNN.GCmean_dev.WW,
     file="./Analysis4paper/Fits1.RData")

# Try with region interaction with stdGC in Density and noneuc:
TNN.GCmean_dev.WW.session <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                              model=list(D~rmeanGC+rmeanGCdev:session, lambda0~Water:Winter, sigma~1, noneuc ~ stdGC:session -1), 
                              details = list(userdist = userdfn1),
                              start = list(noneuc = 1))
TNN.GCmean_dev.WW.session <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                      model=list(D~rmeanGC+rmeanGCdev:session, lambda0~Water:Winter, sigma~1, noneuc ~ stdGC:session -1), 
                                      details = list(userdist = userdfn1),
                                      start = TNN.GCmean_dev.WW.session)

# Try with region interaction with stdGC in Density and noneuc and sigma:
TNN.GCmean_dev.WW.session.sigma <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                      model=list(D~rmeanGC+rmeanGCdev:session, lambda0~Water:Winter, sigma~session, noneuc ~ stdGC:session -1), 
                                      details = list(userdist = userdfn1),
                                      start = list(noneuc = 1))
save(TNN.GCmean_dev.WW.session,file="./Analysis4paper/tempTNN.GCmean_dev.WW.session.sigma.RData")
TNN.GCmean_dev.WW.session.sigma <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                            model=list(D~rmeanGC+rmeanGCdev:session, lambda0~Water:Winter, sigma~session, noneuc ~ stdGC:session -1), 
                                            details = list(userdist = userdfn1),
                                            start = TNN.GCmean_dev.WW.session.sigma)
save(TNN.GCmean_dev.WW.session.sigma,
     TNN.GCmean_dev.WW.session,
     TNN.GCmean_dev.WW,
     file="./Analysis4paper/TNNFits1.RData")

AIC(TNN.GCmean_dev.WW.session.sigma,TNN.GCmean_dev.WW.session,TNN.GCmean_dev.WW)
AIC(TNN.GCmean_dev.WW)$AIC - sum(aics$AIC)
AIC(TNN.GCmean_dev.WW.session)$AIC - sum(aics$AIC)
AIC(TNN.GCmean_dev.WW.session.sigma)$AIC - sum(aics$AIC)
coefficients(Nemegt.NUstdGC)
coefficients(Noyon.NUstdGC)
coefficients(Tost.NUstdGC)
coefficients(TNN.GCmean_dev.WW)
coefficients(TNN.GCmean_dev.WW.session)
coefficients(TNN.GCmean_dev.WW.session.sigma)

region.N(TNN.GCmean_dev.WW.session.sigma)
region.N(TNN.GCmean_dev.WW.session)

# Best individual area spatial model estiamtes:
Nhat.Nemegt = region.N(Nemegt.NUstdGC)
Nhat.Noyon = region.N(Noyon.NUstdGC)
Nhat.Tost = region.N(Tost.NUstdGC)
Nhats = data.frame(
  Nhat=c(Nhat.Tost["R.N","estimate"],
         Nhat.Noyon["R.N","estimate"],
         Nhat.Nemegt["R.N","estimate"]),
  lower=c(Nhat.Tost["R.N","lcl"],
         Nhat.Noyon["R.N","lcl"],
         Nhat.Nemegt["R.N","lcl"]),
  upper=c(Nhat.Tost["R.N","ucl"],
         Nhat.Noyon["R.N","ucl"],
         Nhat.Nemegt["R.N","ucl"])
)
row.names(Nhats)=c("Tost","Noyon","Nemegt")
logNhats = log(Nhats)


# Basic model estiamtes:
Basic.Nemegt = region.N(Nemegt.basic)
Basic.Noyon = region.N(Noyon.basic)
Basic.Tost = region.N(Tost.basic)
Basics = data.frame(
  Nhat=c(Basic.Tost["R.N","estimate"],
         Basic.Noyon["R.N","estimate"],
         Basic.Nemegt["R.N","estimate"]),
  lower=c(Basic.Tost["R.N","lcl"],
          Basic.Noyon["R.N","lcl"],
          Basic.Nemegt["R.N","lcl"]),
  upper=c(Basic.Tost["R.N","ucl"],
          Basic.Noyon["R.N","ucl"],
          Basic.Nemegt["R.N","ucl"])
)
row.names(Basics)=c("Tost","Noyon","Nemegt")
logBasics = log(Basics)

# Do some plots
meanGC = c(unique(covariates(TostMask)$rmeanGC),
           unique(covariates(NoyonMask)$rmeanGC),
           unique(covariates(NemegtMask)$rmeanGC))
names(meanGC) = c("Tost","Noyon","Nemegt")
# Estimates on natrual scale:
ylim = range(Nhats,Basics)
plot(meanGC,Nhats$Nhat,pch=19,ylim=ylim)
lines(meanGC,Nhats$Nhat,lty=2)
segments(meanGC,Nhats$lower,meanGC,Nhats$upper)
offset=0.1
points(meanGC+offset,Basics$Nhat,pch=19,ylim=ylim,col="gray")
lines(meanGC+offset,Basics$Nhat,lty=2,col="gray")
segments(meanGC+offset,Basics$lower,meanGC+offset,Basics$upper,col="gray")
# Estimates on log scale:
ylim = range(logNhats,logBasics)
plot(meanGC,logNhats$Nhat,pch=19,ylim=ylim)
lines(meanGC,logNhats$Nhat,lty=2)
segments(meanGC,logNhats$lower,meanGC,logNhats$upper)
offset=0.1
points(meanGC+offset,logBasics$Nhat,pch=19,ylim=ylim,col="gray")
lines(meanGC+offset,logBasics$Nhat,lty=2,col="gray")
segments(meanGC+offset,logBasics$lower,meanGC+offset,logBasics$upper,col="gray")























# Test a few models
# model with region mean effect only:
GC.TNN.hhn<-secr.fit(all.data.TNN, 
                  model=list(D~rmeanGC, lambda0~1, sigma~1), detectfn="HHN", 
                  mask=list(TostMask, NoyonMask, NemegtMask))
#save(GC.TNN.hhn,file="./Tost_Noyon_Nemegt/GC.TNN.hhn.RData")

# model with region and GC effect: DOES NOT FIT
#GC.reg.TNN.hhn<-secr.fit(all.data.TNN, 
#                         model=list(D~session+GC, lambda0~1, sigma~1), detectfn="HHN", 
#                         mask=list(TostMask, NoyonMask, NemegtMask))
#save(GC.reg.TNN.hhn,file="./Tost_Noyon_Nemegt/GC.reg.TNN.hhn.RData")

# model with region, region mean and GC deviation effect:
GC.dev.TNN.hhn<-secr.fit(all.data.TNN, 
                     model=list(D~rmeanGC+rmeanGCdev, lambda0~1, sigma~1), detectfn="HHN", 
                     mask=list(TostMask, NoyonMask, NemegtMask))
#save(GC.dev.TNN.hhn,file="./Tost_Noyon_Nemegt/GC.dev.TNN.hhn.RData")

GC.dev_region.TNN.hhn<-secr.fit(all.data.TNN, 
                         model=list(D~rmeanGC+rmeanGCdev+rmeanGCdev:session, lambda0~1, sigma~1), detectfn="HHN", 
                         mask=list(TostMask, NoyonMask, NemegtMask))
#save(GC.dev_region.TNN.hhn,file="./Tost_Noyon_Nemegt/GC.dev_region.TNN.hhn.RData")

region.GC.dev_region.TNN.hhn<-secr.fit(all.data.TNN, 
                                       model=list(D~session+rmeanGC+rmeanGCdev+rmeanGCdev:session, lambda0~1, sigma~1), detectfn="HHN", 
                                       mask=list(TostMask, NoyonMask, NemegtMask))
#save(region.GC.dev_region.TNN.hhn,file="./Tost_Noyon_Nemegt/region.GC.dev_region.TNN.hhn.RData")

region.GC.dev.region.TNN.hhn<-secr.fit(all.data.TNN, 
                                model=list(D~session+rmeanGC+rmeanGCdev, lambda0~1, sigma~1), detectfn="HHN", 
                                mask=list(TostMask, NoyonMask, NemegtMask))
#save(region.GC.dev.region.TNN.hhn,file="./Tost_Noyon_Nemegt/region.GC.dev.region.TNN.hhn")

load("./Tost_Noyon_Nemegt/GC.TNN.hhn.RData")
#load("./Tost_Noyon_Nemegt/GC.reg.TNN.hhn.RData")
load("./Tost_Noyon_Nemegt/GC.dev.TNN.hhn.RData")
load("./Tost_Noyon_Nemegt/GC.dev_region.TNN.hhn.RData")
load("./Tost_Noyon_Nemegt/region.GC.dev_region.TNN.hhn.RData")
load("./Tost_Noyon_Nemegt/region.GC.dev.region.TNN.hhn")

AIC(
  GC.TNN.hhn,
#  GC.reg.TNN.hhn,
  GC.dev.TNN.hhn,
  GC.dev_region.TNN.hhn,
  region.GC.dev.region.TNN.hhn,
  region.GC.dev_region.TNN.hhn
  )


# --------- Some NonEuc fits -----------
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


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# IGNORE THE INDIVIDUAL AREA FITS FOR NOW - UNTIL RESOLVE NUMBER OF TRAPS IN NEMEGT
# (NEMEGT DATA HAS 6 MORE TRAPS THAN IN COMBINED DATA)
# Nemegt: D(GC); noneuc(GC)l lambda0(Water:Winter):
# ---------------------------
Nemegt.NUstdGC<-secr.fit(all.data.Nemegt, detectfn="HHN", mask=NemegtMask,
                                       model=list(D~stdGC, lambda0~Water:Winter, sigma~1, noneuc ~ stdGC -1), 
                                       details = list(userdist = userdfn1),
                                       start = list(noneuc = 1))
# replace noneuc(stdGC) with noneuc(BC)
Nemegt.NUBC<-secr.fit(all.data.Nemegt, detectfn="HHN", mask=NemegtMask,
                      model=list(D~stdGC, lambda0~Water:Winter, sigma~1, noneuc ~ BC -1), 
                      details = list(userdist = userdfn1),
                      start = list(noneuc = 1))
AIC(Nemegt.NUBC,Nemegt.NUstdGC)

# Noyon: D(GC); noneuc(BC):
# ---------------------------
Noyon.NUBC<-secr.fit(all.data.Noyon, detectfn="HHN", mask=NoyonMask,
                     model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ BC -1), 
                     details = list(userdist = userdfn1),
                     start = list(noneuc = 1))
# replace noneuc(BC) with noneuc(stdGC)
Noyon.NUstdGC<-secr.fit(all.data.Noyon, detectfn="HHN", mask=NoyonMask,
                     model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                     details = list(userdist = userdfn1),
                     start = list(noneuc = 1))
AIC(Noyon.NUBC,Noyon.NUstdGC)

# Tost: D(GC); noneuc(BC):
# ---------------------------
Tost.NUBC<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask,
                     model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ BC -1), 
                     details = list(userdist = userdfn1),
                     start = list(noneuc = 1))
# replace noneuc(BC) with noneuc(stdGC)
Tost.NUstdGC<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask,
                    model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                    details = list(userdist = userdfn1),
                    start = list(noneuc = 1))
AIC(Tost.NUBC,Tost.NUstdGC)

#save(Nemegt.NUBC,Nemegt.NUstdGC,Noyon.NUBC,Noyon.NUstdGC,Tost.NUBC,Tost.NUstdGC,
#     file="NUBC_NUstdGC_fits.RData")
load("NUBC_NUstdGC_fits.RData")

# Go with the NUstdGC models
# --------------------------
# Look at parameter estimates
coefficients(Nemegt.NUstdGC) # 
coefficients(Noyon.NUstdGC) # 
coefficients(Tost.NUstdGC) # 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -





# Try joint model wth the region-specific mean
# Do not put Water in lambda0's for now
TNN.GCmean_dev <- secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                       model=list(D~rmeanGC+rmeanGCdev, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                       details = list(userdist = userdfn1),
                       start = list(noneuc = 1))
save(TNN.GCmean_dev,file="TNN.GCmean_dev.RData")
aics = AIC(Nemegt.NUstdGC,Noyon.NUstdGC,Tost.NUstdGC)
sum(aics$AIC)
AIC(TNN.GCmean_dev)$AIC

# Need to use session instead of region as factor
TNN.GCmean_devregion <- secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                 model=list(D~rmeanGC+rmeanGCdev:session, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                                 details = list(userdist = userdfn1),
                                 start = list(noneuc = 1))
AIC(TNN.GCmean_dev,TNN.GCmean_devregion)

TNN.dev_GCmean_devregion <- secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                 model=list(D~session+rmeanGC+rmeanGCdev:session, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                                 details = list(userdist = userdfn1),
                                 start = list(noneuc = 1))
AIC(TNN.GCmean_dev,TNN.GCmean_devregion,TNN.dev_GCmean_devregion)
save(TNN.GCmean_dev,TNN.GCmean_devregion,TNN.dev_GCmean_devregion,file="TNN.noWfits.RData")


# Same, but with Water in lambda0:
TNN.GCmean_dev.W <- secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                           model=list(D~rmeanGC+rmeanGCdev, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                           details = list(userdist = userdfn1),
                           start = list(noneuc = 1))

TNN.GCmean_devregion.W <- secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                 model=list(D~rmeanGC+rmeanGCdev:session, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                                 details = list(userdist = userdfn1),
                                 start = list(noneuc = 1))

TNN.dev_GCmean_devregion.W <- secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                     model=list(D~session+rmeanGC+rmeanGCdev:session, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                                     details = list(userdist = userdfn1),
                                     start = list(noneuc = 1))


# Same, but with Water:Winter in lambda0:
TNN.GCmean_dev.WW <- secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                             model=list(D~rmeanGC+rmeanGCdev, lambda0~Water:Winter, sigma~1, noneuc ~ stdGC -1), 
                             details = list(userdist = userdfn1),
                             start = list(noneuc = 1))

TNN.GCmean_devregion.WW <- secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                   model=list(D~rmeanGC+rmeanGCdev:session, lambda0~Water:Winter, sigma~1, noneuc ~ stdGC -1), 
                                   details = list(userdist = userdfn1),
                                   start = list(noneuc = 1))

TNN.dev_GCmean_devregion.WW <- secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                       model=list(D~session+rmeanGC+rmeanGCdev:session, lambda0~Water:Winter, sigma~1, noneuc ~ stdGC -1), 
                                       details = list(userdist = userdfn1),
                                       start = list(noneuc = 1))

AIC(TNN.GCmean_dev,TNN.GCmean_devregion,TNN.dev_GCmean_devregion,
    TNN.GCmean_dev.W,TNN.GCmean_devregion.W,TNN.dev_GCmean_devregion.W,
    TNN.GCmean_dev.WW,TNN.GCmean_devregion.WW,TNN.dev_GCmean_devregion.WW)



# Repeat all above, with all.data.TNN_R instead of all.data.TNN:
# -------------------------------------------------------------
# Try joint model wth the region-specific mean
# Do not put Water in lambda0's for now
TNN_R.GCmean_dev <- secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                           model=list(D~rmeanGC+rmeanGCdev, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                           details = list(userdist = userdfn1),
                           start = list(noneuc = 1))
save(TNN_R.GCmean_dev,file="TNN_R.GCmean_dev.RData")
aics = AIC(Nemegt.NUstdGC,Noyon.NUstdGC,Tost.NUstdGC)
sum(aics$AIC)
AIC(TNN_R.GCmean_dev)$AIC

# Need to use session instead of region as factor
TNN_R.GCmean_devregion <- secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                 model=list(D~rmeanGC+rmeanGCdev:session, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                                 details = list(userdist = userdfn1),
                                 start = list(noneuc = 1))
AIC(TNN_R.GCmean_dev,TNN_R.GCmean_devregion)

TNN_R.dev_GCmean_devregion <- secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                     model=list(D~session+rmeanGC+rmeanGCdev:session, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                                     details = list(userdist = userdfn1),
                                     start = list(noneuc = 1))
AIC(TNN_R.GCmean_dev,TNN_R.GCmean_devregion,TNN_R.dev_GCmean_devregion)
save(TNN_R.GCmean_dev,TNN_R.GCmean_devregion,TNN_R.dev_GCmean_devregion,file="TNN_R.noWfits.RData")


# Same, but with Water in lambda0:
TNN_R.GCmean_dev.W <- secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                             model=list(D~rmeanGC+rmeanGCdev, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                             details = list(userdist = userdfn1),
                             start = list(noneuc = 1))

TNN_R.GCmean_devregion.W <- secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                   model=list(D~rmeanGC+rmeanGCdev:session, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                                   details = list(userdist = userdfn1),
                                   start = list(noneuc = 1))

TNN_R.dev_GCmean_devregion.W <- secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                       model=list(D~session+rmeanGC+rmeanGCdev:session, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                                       details = list(userdist = userdfn1),
                                       start = list(noneuc = 1))


# Same, but with Water:Winter in lambda0:
TNN_R.GCmean_dev.WW <- secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                              model=list(D~rmeanGC+rmeanGCdev, lambda0~Water:Winter, sigma~1, noneuc ~ stdGC -1), 
                              details = list(userdist = userdfn1),
                              start = list(noneuc = 1))

TNN_R.GCmean_devregion.WW <- secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                    model=list(D~rmeanGC+rmeanGCdev:session, lambda0~Water:Winter, sigma~1, noneuc ~ stdGC -1), 
                                    details = list(userdist = userdfn1),
                                    start = list(noneuc = 1))

TNN_R.dev_GCmean_devregion.WW <- secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                        model=list(D~session+rmeanGC+rmeanGCdev:session, lambda0~Water:Winter, sigma~1, noneuc ~ stdGC -1), 
                                        details = list(userdist = userdfn1),
                                        start = list(noneuc = 1))

AIC(TNN_R.GCmean_dev,TNN_R.GCmean_devregion,TNN_R.dev_GCmean_devregion,
    TNN_R.GCmean_dev.W,TNN_R.GCmean_devregion.W,TNN_R.dev_GCmean_devregion.W,
    TNN_R.GCmean_dev.WW,TNN_R.GCmean_devregion.WW,TNN_R.dev_GCmean_devregion.WW)






















# ------------------------------ old code -------------------------------
# Standarize GRIDCODE (in stdGC) and BINCODE (in stdBC) on mask
# ------------------------------------------------------------------------
summary(covariates(NemegtMask))
covariates(NemegtMask)$stdGC = scale(covariates(NemegtMask)$GRIDCODE)
covariates(NemegtMask)$stdBC = scale(covariates(NemegtMask)$BINCODE)
summary(covariates(NemegtMask))
names(covariates(NemegtMask))

summary(covariates(NoyonMask))
covariates(NoyonMask)$stdGC = scale(covariates(NoyonMask)$GRIDCODE)
covariates(NoyonMask)$stdBC = scale(covariates(NoyonMask)$BINCODE)
summary(covariates(NoyonMask))
names(covariates(NoyonMask))

summary(covariates(TostMask))
covariates(TostMask)$stdGC = scale(covariates(TostMask)$GRIDCODE)
covariates(TostMask)$stdBC = scale(covariates(TostMask)$BINCODE)
summary(covariates(TostMask))
names(covariates(TostMask))

# Plot sessions together:
pdf("Allregions.pdf",h=5,w=10)
plot(bbxlim,bbylim,xlim=xlim,ylim=ylim,xlab="",ylab="",bty="n",type="n",xaxt="n",yaxt="n",asp=1) #Generates Error!
# Plot the terrain
plot(NoyonMask, covariate="stdGC", contour = FALSE, col = terrain.colors(16), legend = FALSE, add = TRUE)
plot(NemegtMask, covariate="stdGC", contour=FALSE, col=terrain.colors(16), legend = FALSE, add = TRUE)
plot(TostMask, covariate="stdGC", contour=FALSE, col=terrain.colors(16), legend = FALSE, add=TRUE)
# Add the traps
plot(traps(all.data.TNN)[[1]],add=TRUE,detpar=list(col="black",pch="+"))
plot(traps(all.data.TNN)[[2]],add=TRUE,detpar=list(col="black",pch="+"))
plot(traps(all.data.TNN)[[3]],add=TRUE,detpar=list(col="black",pch="+"))
# and the borders over the top
plot(Nemegtboundary,add=TRUE)#,border=3)
plot(Noyonboundary,add=TRUE)#,border=2)
plot(Tostboundary,add=TRUE)#,border=1)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
dev.off()

TNN.hhn<-secr.fit(all.data.TNN, 
                  model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", 
                  mask=list(TostMask, NoyonMask, NemegtMask))

TNN.hhnR<-secr.fit(all.data.TNN_R, 
                  model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", 
                  mask=list(TostMask, NoyonMask, NemegtMask))

coefficients(TNN.hhn)
coefficients(TNN.hhnR)


####Is the code below correct? Can't seem to add different covariates for different sessions (here areas)  
TNN.hhn.DRgd<-secr.fit(all.data.TNN, model = list(D~stdGC, lambda0~1, sigma~1), detectfn="HHN",
                       mask = list(TostMask, NoyonMask, NemegtMask))

TNN.hhn.DRgd.DetTopo10<-secr.fit(all.data.TNN, model = list(D~stdGC, lambda0~Topo, sigma~1), detectfn="HHN",
                       mask = list(TostMask, NoyonMask, NemegtMask))


predict(TNN.hhn.DRgd)
TNNSurface.DRgd<-predictDsurface(TNN.hhn.DRgd)
# When you have sessions, you have to plot by session:
windows()
plot(TNNSurface.DRgd[[1]],asp=1,contour=FALSE,asp=1) 
plot(TNNSurface.DRgd[[2]],asp=1,contour=FALSE,asp=1) 
plot(TNNSurface.DRgd[[3]],asp=1,contour=FALSE,asp=1) 

# Try with session factor
# -----------------------
sess = as.factor(1:3)
TNN.hhn.DRgd.sess<-secr.fit(all.data.TNN, model = list(D~stdGC+sfac, lambda0~1, sigma~1), detectfn="HHN",
                       mask = list(TostMask, NoyonMask, NemegtMask),sessioncov=data.frame(sfac=sess))

TNN.hhn.DRgd.sessR<-secr.fit(all.data.TNN_R, model = list(D~stdGC+sfac, lambda0~1, sigma~1), detectfn="HHN",
                            mask = list(TostMask, NoyonMask, NemegtMask),sessioncov=data.frame(sfac=sess))

TNN.hhn.DRgd.sess.DetW<-secr.fit(all.data.TNN, model = list(D~stdGC+sfac, lambda0~sfac*Water, sigma~1), detectfn="HHN",
                            mask = list(TostMask, NoyonMask, NemegtMask),sessioncov=data.frame(sfac=sess))

TNN.hhn.DRgd.sess.DetWR<-secr.fit(all.data.TNN_R, model = list(D~stdGC+sfac, lambda0~sfac*Water, sigma~1), detectfn="HHN",
                                 mask = list(TostMask, NoyonMask, NemegtMask),sessioncov=data.frame(sfac=sess))

TNN.hhn.DRgd.DetTopo10W<-secr.fit(all.data.TNN, model = list(D~stdGC, lambda0~Topo+Water, sigma~1), detectfn="HHN",
                                mask = list(TostMask, NoyonMask, NemegtMask))
TNN.hhn.DRgd.DetTopo10WR<-secr.fit(all.data.TNN_R, model = list(D~stdGC, lambda0~Topo+Water, sigma~1), detectfn="HHN",
                                  mask = list(TostMask, NoyonMask, NemegtMask))

TNN.hhn.DRgd.sess_interact<-secr.fit(all.data.TNN, model = list(D~stdGC*sfac, lambda0~1, sigma~1), detectfn="HHN",
                            mask = list(TostMask, NoyonMask, NemegtMask),sessioncov=data.frame(sfac=sess))


AIC(TNN.hhn.DRgd.DetTopo10WR, TNN.hhn.DRgd.sess.DetWR, TNN.hhn.DRgd.sessR, TNN.hhnR)

coefficients(TNN.hhn.DRgd.DetTopo10WR)

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
TNN.hhn.DHab.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                            model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                            details = list(userdist = userdfn1),
                               start = list(noneuc = 1)) #-1 gets rid of the intercept
TNN.hhn.DHab.nonUR<-secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                            model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                            details = list(userdist = userdfn1),
                            start = list(noneuc = 1)) #-1 gets rid of the intercept

# Model with stdGC in noneuc, Topo in Detect:
# ---------------------------
TNN.hhn.DHab.DetTopo10.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                            model=list(D~stdGC, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                            details = list(userdist = userdfn1),
                            start = list(noneuc = 1)) #-1 gets rid of the intercept

TNN.hhn.DHab.DetTopo10.nonU1<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                                      model=list(D~stdGC, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                                      details = list(userdist = userdfn1),
                                      start = list(noneuc = 1)) #-1 gets rid of the intercept


TNN.hhn.DHab.DetTopo10.nonUR<-secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                                      model=list(D~stdGC, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                                      details = list(userdist = userdfn1),
                                      start = list(noneuc = 1)) #-1 gets rid of the intercept

# Model with BC in noneuc, Topo in Detect, GC in density:
# ---------------------------

TNN.hhn.DGC.DetTopo10.nonUGBR<-secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                                       model=list(D~stdGC, lambda0~Topo, sigma~1, noneuc ~ stdBC -1), 
                                       details = list(userdist = userdfn1),
                                       start = list(noneuc = 1)) #-1 gets rid of the intercept

# Model with GC in noneuc, Topo in Detect, GB in density:
# ---------------------------

TNN.hhn.DGB.DetTopo10.nonUR<-secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                                        model=list(D~stdBC, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                                        details = list(userdist = userdfn1),
                                        start = list(noneuc = 1)) #-1 gets rid of the intercept


coefficients(TNN.hhn.DHab.DetTopo10.nonU)
coefficients(TNN.hhn.DHab.DetTopo10.nonUR)
coefficients(TNN.hhn.DGC.DetTopo10.nonUGBR)
coefficients(TNN.hhn.DGB.DetTopo10.nonUR)

summary(all.data.TNN)
summary(all.data.TNN_R)
# Same model as above, with reduced dataset


# Model with stdGC in noneuc, Topo in Detect per session:
# ---------------------------
TNN.hhn.DHab.DetToposess.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                                      model=list(D~stdGC, lambda0~Topo*sfac, sigma~1, noneuc ~ stdGC -1), 
                                      details = list(userdist = userdfn1),sessioncov=data.frame(sfac=sess),
                                      start = list(noneuc = 1)) #-1 gets rid of the intercept
TNN.hhn.DHab.DetToposess.nonUR<-secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                                        model=list(D~stdGC, lambda0~Topo*sfac, sigma~1, noneuc ~ stdGC -1), 
                                        details = list(userdist = userdfn1),sessioncov=data.frame(sfac=sess),
                                        start = list(noneuc = 1)) #-1 gets rid of the intercept


coefficients(TNN.hhn.DHab.DetToposess.nonU)


# Model with D->Rgd session interact, noneuc->rgd, Topo in lam0:
# ---------------------------
TNN.hhn.DHab.S.DetTopo10.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                                      model=list(D~stdGC*sfac, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                                      details = list(userdist = userdfn1),sessioncov=data.frame(sfac=sess),
                                      start = list(noneuc = 1)) #-1 gets rid of the intercept

TNN.hhn.DHab.S.DetTopo10.nonUR<-secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                                        model=list(D~stdGC*sfac, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                                        details = list(userdist = userdfn1),sessioncov=data.frame(sfac=sess),
                                        start = list(noneuc = 1)) #-1 gets rid of the intercept

# Model with D->Rgd session add, noneuc->rgd, Topo in lam0:
# ---------------------------
TNN.hhn.DHab_S.DetTopo10.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                                        model=list(D~stdGC+sfac, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                                        details = list(userdist = userdfn1),sessioncov=data.frame(sfac=sess),
                                        start = list(noneuc = 1)) #-1 gets rid of the intercept

TNN.hhn.DHab_S.DetTopo10.nonUR<-secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                                        model=list(D~stdGC+sfac, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                                        details = list(userdist = userdfn1),sessioncov=data.frame(sfac=sess),
                                        start = list(noneuc = 1)) #-1 gets rid of the intercept

coefficients(TNN.hhn.DHab.DetTopo10.nonU)
# Model with stdGC in noneuc, Topo+Water in Detect:
# ---------------------------
TNN.hhn.DHab.LamTopoWat.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                                      model=list(D~stdGC, lambda0~Topo+Water, sigma~1, noneuc ~ stdGC -1), 
                                      details = list(userdist = userdfn1),
                                      start = list(noneuc = 1)) #-1 gets rid of the intercept

TNN.hhn.DHab.LamTopoWat.nonUR<-secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                                       model=list(D~stdGC, lambda0~Topo+Water, sigma~1, noneuc ~ stdGC -1), 
                                       details = list(userdist = userdfn1),
                                       start = list(noneuc = 1)) #-1 gets rid of the intercept

coefficients(TNN.hhn.DHab.LamTopoWat.nonU)

# Model with constt D in noneuc:
# ---------------------------
TNN.hhn.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                            model=list(D~1, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                            details = list(userdist = userdfn1), start = list(noneuc = 1)) #-1 gets rid of the intercept

TNN.hhn.nonUR<-secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                       model=list(D~1, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                       details = list(userdist = userdfn1), start = list(noneuc = 1)) #-1 gets rid of the intercept

TNNAIC2x<-AIC(TNN.hhn.nonU, TNN.hhn.DHab.LamTopoWat.nonU, TNN.hhn.DHab.DetTopo10.nonU,TNN.hhn.DRgd, 
    TNN.hhn.DHab.nonU,TNN.hhn.DRgd.sess.DetW,TNN.hhn.DRgd.DetTopo10W, TNN.hhn.DHab_S.DetTopo10.nonU,
    TNN.hhn.DHab.S.DetTopo10.nonU, TNN.hhn.DHab.DetToposess.nonU)

TNNAIC2xR<- AIC(TNN.hhn.DRgd.sessR, TNN.hhn.DRgd.sess.DetWR, TNN.hhn.DRgd.DetTopo10WR, TNN.hhn.DHab.nonUR, 
                TNN.hhn.DHab.DetTopo10.nonUR, TNN.hhn.DHab.DetToposess.nonUR, 
                TNN.hhn.DHab.S.DetTopo10.nonUR, TNN.hhn.DHab_S.DetTopo10.nonUR, 
                TNN.hhn.DHab.LamTopoWat.nonUR, TNN.hhn.nonUR,TNN.hhn.DGC.DetTopo10.nonUGBR,
                TNN.hhn.DGB.DetTopo10.nonUR)
TNNAIC2x
TNNAIC2xR

coefficients(TNN.hhn.DGC.DetTopo10.nonUGBR)
coefficients(TNN.hhn.DHab.DetTopo10.nonUR)
coefficients(TNN.hhn.DHab.LamTopoWat.nonUR)
coefficients(TNN.hhn.DHab.DetToposess.nonUR)
coefficients(TNN.hhn.DHab_S.DetTopo10.nonUR)
#Very very long time to model (some models took up to 4 hours, but end results look inverse now!)

write.csv(TNNAIC2x, file = "TNNAIC2x.csv")
coefficients(TNN.hhn.DHab.DetTopo10.nonU)
coefficients(TNN.hhn.DRgd.sess.DetW)

save(TNN.hhn.nonU, TNN.hhn.DHab.LamTopoWat.nonU, TNN.hhn.DHab.DetTopo10.nonU,TNN.hhn.DRgd, 
     TNN.hhn.DHab.nonU,TNN.hhn.DRgd.sess.DetW,TNN.hhn.DRgd.DetTopo10W, TNN.hhn.DHab_S.DetTopo10.nonU,
     TNN.hhn.DHab.S.DetTopo10.nonU, TNN.hhn.DHab.DetToposess.nonU, file="./Tost_Noyon_Nemegt/TNN-nonEuc-fits2x.RData")

save(TNN.hhn.DRgd.sessR, TNN.hhn.DRgd.sess.DetWR, TNN.hhn.DRgd.DetTopo10WR, TNN.hhn.DHab.nonUR, 
         TNN.hhn.DHab.DetTopo10.nonUR, TNN.hhn.DHab.DetToposess.nonUR, TNN.hhn.DHab.S.DetTopo10.nonUR, 
         TNN.hhn.DHab_S.DetTopo10.nonUR, TNN.hhn.DHab.LamTopoWat.nonUR, TNN.hhn.DGB.DetTopo10.nonUR,
        TNN.hhn.nonUR,TNN.hhn.DGC.DetTopo10.nonUGBR, file = "./Tost_Noyon_Nemegt/TNN-nonEuc-fits2xR.RData")

load("./Tost_Noyon_Nemegt/TNN-NonEuc-fits2x.RData")

load("./Tost_Noyon_Nemegt/TNN-NonEuc-fits2xR.RData")

# Model with D dependent on Rgd OR session in noneuc:
# ---------------------------
TNN.hhn.DRgd_sess.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask),
                                 model=list(D~stdGC+sfac, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                                 sessioncov=data.frame(sfac=sess),details = list(userdist = userdfn1), 
                                 start = list(noneuc = 1)) #-1 gets rid of the intercept


# Model with D dependent on Rgd & different for each session in noneuc:
# ---------------------------
TNN.hhn.DRgd.sess.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask),
                                 model=list(D~stdGC*sfac, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                                 sessioncov=data.frame(sfac=sess),details = list(userdist = userdfn1), 
                                 start = list(noneuc = 1)) #-1 gets rid of the intercept

# Model with D dependent on Rgd &  different for each session, lambda on water in noneuc:
# ---------------------------
TNN.hhn.DRgd.sess.D.W_sess.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask),
                                          model=list(D~stdGC*sfac, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                                          sessioncov=data.frame(sfac=sess), details = list(userdist = userdfn1), 
                                          start = list(noneuc = 1)) #-1 gets rid of the intercept

# Model with D dependent on Rgd, lambda on water & BOTH different for each session in noneuc:
# ---------------------------
TNN.hhn.DRgd.sess.D.W.sess.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                            model=list(D~stdGC*sfac, lambda0~Water*sfac, sigma~1, noneuc ~ stdGC -1), 
                            sessioncov=data.frame(sfac=sess),
                            details = list(userdist = userdfn1), start = list(noneuc = 1)) #-1 gets rid of the intercept

# Model with D dependent on Rgd, lambda on water+topo & BOTH different for each session in noneuc:
# ---------------------------
TNN.hhn.DRgd.sess.D_W_sess.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                                          model=list(D~stdGC*sfac, lambda0~(Water+Topo)*sfac, sigma~1, noneuc ~ stdGC -1), 
                                          sessioncov=data.frame(sfac=sess),
                                          details = list(userdist = userdfn1), start = list(noneuc = 1)) #-1 gets rid of the intercept
#Took 1938 iterations to converge!!!

TNNAIC=AIC(TNN.hhn.detTopo10, TNN.hhn.detTopo01, TNN.hhn, TNN.hhn.DRgd, TNN.hhn.DRgd.DetWater, 
    TNN.hhn.DRgd.DetTopo10W, TNN.hhn.DHab.nonU, TNN.hhn.nonU, TNN.hhn.DRgd.sess.D.W.sess.nonU,TNN.hhn.DRgd.sess.D.W_sess.nonU,
    TNN.hhn.DRgd.sess.nonU,TNN.hhn.DRgd_sess.nonU, TNN.hhn.DRgd.sess.D_W_sess.nonU, TNN.hhn.DRgd.DetTopo10)
write.csv(TNNAIC, file = "AICTNNx.csv")
TNNAIC = read.csv("AICTNNx.csv")

getwd()


coefficients(TNN.hhn.DHab.nonU)
coefficients(TNN.hhn.DRgd.sess.D.W.sess.nonU)
coefficients(TNN.hhn.DRgd.sess.D_W_sess.nonU)
predict(TNN.hhn.DRgd.DetWater)
predict(TNN.hhn.DRgd.sess.D.W.sess.nonU)

TNNSurface.DRgdX<-predictDsurface(TNN.hhn.DRgd.DetWater)
windows()
plot(TNNSurface.DRgd,asp=1,contour=FALSE) #This generates an error. Something I am doing wrong here it seems...

save(TNN.hhn, TNN.hhn.DRgd, TNN.hhn.DRgd.sess, TNN.hhn.DRgd.sess.DetW, TNN.hhn.DRgd.DetTopo10W, TNN.hhn.DHab.nonU, 
     TNN.hhn.nonU, TNN.hhn.DRgd_sess.nonU, TNN.hhn.DRgd.sess.nonU, TNN.hhn.DRgd.sess.D.W_sess.nonU, 
     TNN.hhn.DRgd.sess.D.W.sess.nonU, TNN.hhn.DRgd.sess.D_W_sess.nonU, file="./Tost_Noyon_Nemegt/TNN-nonEuc-fitsx.RData")
# load fitted objects:
load("./Tost_Noyon_Nemegt/TNN-nonEuc-fitsx.RData")


# How to get region.N for each of the 3 areas?
# region.N for each of the 3 areas:
region.N(TNN.hhn.DRgd.sess,region=TostMask,session="1")
region.N(TNN.hhn.DRgd.sess,region=NoyonMask,session="2")
region.N(TNN.hhn.DRgd.sess,region=NemegtMask,session="3")

region.N(TNN.hhn.DRgd.sess,region=TostMask,session="1")
region.N(TNN.hhn.DRgd.sess,region=NoyonMask,session="2")
region.N(TNN.hhn.DRgd.sess,region=NemegtMask,session="3")

region.N(TNN.hhn.DRgd.sess.D.W.sess.nonU,region=TostMask,session="1")
region.N(TNN.hhn.DRgd.sess.D.W.sess.nonU,region=NoyonMask,session="2")
region.N(TNN.hhn.DRgd.sess.D.W.sess.nonU, region=NemegtMask,session="3")

coefficients(TNN.hhn.DRgd.sess)
predict(TNN.hhn.DRgd.sess)
TNNSurface.DRgd.sess<-predictDsurface(TNN.hhn.DRgd.sess)
TNNSurface.Top<-predictDsurface(TNN.hhn.DRgd.sess.D.W.sess.nonU)

TNNSurfacex.Top<-predictDsurface(TNN.hhn.DHab.nonU)
plot(TNNSurfacex.Top[[1]], asp=1,contour=FALSE)

# When you have sessions, you have to plot by session:
windows()
plot(TNNSurface.Top[[1]],asp=1,contour=FALSE) 
plot(TNNSurface.Top[[2]],asp=1,contour=FALSE) 
plot(TNNSurface.Top[[3]],asp=1,contour=FALSE) 

plot(TNNSurface.DRgd.sess[[1]],asp=1,contour=FALSE) 
plot(TNNSurface.DRgd.sess[[2]],asp=1,contour=FALSE) 
plot(TNNSurface.DRgd.sess[[3]],asp=1,contour=FALSE) 


plot(TNN.hhn.DRgd.sess.D.W.sess.nonU[[1]], asp=1, contour=FALSE)
plot(TNN.hhn.DRgd.sess.D.W.sess.nonU[[2]], asp=1, contour=FALSE)
plot(TNN.hhn.DRgd.sess.D.W.sess.nonU[[3]], asp=1, contour=FALSE)



# --------------------------- David additions March 2017 ----------------------------------

# David trial Model with stdGC in noneuc and smooth (k=3) stdGC effect on D and lambda0~sfac*Water:
# ---------------------------
TNN.hhn.DHabS3.lambdaSfacWater.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                              model=list(D~s(stdGC,k=3), lambda0~sfac*Water, sigma~1, noneuc ~ stdGC -1), 
                              sessioncov=data.frame(sfac=sess), details = list(userdist = userdfn1),
                              start = list(noneuc = 1)) #-1 gets rid of the intercept
# Some convergence problems with the above model; need to refit starting from this fit, but not yet done that.

# David trial Model with stdGC in noneuc and smooth (k=3) stdGC effect on D:
# ---------------------------
TNN.hhn.DHabS3.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                              model=list(D~s(stdGC,k=3), lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                              details = list(userdist = userdfn1),
                              start = list(noneuc = 1)) #-1 gets rid of the intercept
# Fit again, using estimates from above as starting values:
TNN.hhn.DHabS3a.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask, NoyonMask,NemegtMask), 
                              model=list(D~s(stdGC,k=3), lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                              details = list(userdist = userdfn1),
                              start = TNN.hhn.DHabS3.nonU) 

TNN.hhn.DHabS3.nonU = TNN.hhn.DHabS3a.nonU

save(TNN.hhn, TNN.hhn.DRgd, TNN.hhn.DRgd.sess, TNN.hhn.DRgd.sess.DetW, TNN.hhn.DRgd.DetTopo10W, TNN.hhn.DHab.nonU, 
     TNN.hhn.nonU, TNN.hhn.DRgd_sess.nonU, TNN.hhn.DRgd.sess.nonU, TNN.hhn.DRgd.sess.D.W_sess.nonU, 
     TNN.hhn.DRgd.sess.D.W.sess.nonU, TNN.hhn.DRgd.sess.D_W_sess.nonU,TNN.hhn.DHabS3.nonU, 
     file="./Tost_Noyon_Nemegt/TNN-nonEuc-fitsx1.RData")
# load fitted objects:
load("./Tost_Noyon_Nemegt/TNN-nonEuc-fitsx1.RData")

AIC(TNN.hhn.DRgd, TNN.hhn.DRgd.sess, TNN.hhn.DRgd.sess.DetW, TNN.hhn.DRgd.DetTopo10W, TNN.hhn.DHab.nonU, 
    TNN.hhn.nonU, TNN.hhn.DRgd_sess.nonU, TNN.hhn.DRgd.sess.nonU, TNN.hhn.DRgd.sess.D.W_sess.nonU, 
    TNN.hhn.DRgd.sess.D.W.sess.nonU, TNN.hhn.DRgd.sess.D_W_sess.nonU,TNN.hhn.DHabS3.nonU)

# individual regions:
AICTostx = read.csv("AICTostx.csv")
AICTostx
AICNoyonx = read.csv("AICNoyonx.csv")
AICNoyonx
AICNemegtx = read.csv("AICNemegtx.csv")
AICNemegtx

# Best separate AICs:
aics = c(465.603,436.753,286.066)
aicsum = sum(aics);aicsum # Combined AIC for models fitted separately to each stratum

# Abundance and plots for given model
#fit = TNN.hhn.DHabS3.lambdaSfacWater.nonU
fit = TNN.hhn.DHabS3.nonU

# region.N for each of the 3 areas:
region.N(fit,region=TostMask,session="1")
region.N(fit,region=NoyonMask,session="2")
region.N(fit,region=NemegtMask,session="3")

coefficients(fit)
predict(fit)
fitpred<-predictDsurface(fit)

# When you have sessions, you have to plot by session:
windows()
plot(fitpred[[1]],asp=1,contour=FALSE) 
plot(fitpred[[2]],asp=1,contour=FALSE) 
plot(fitpred[[3]],asp=1,contour=FALSE) 

# difficult to see on natural scale - dominated by massive density spikes at a few points (high gridcode)
# so look on the log density scale:
logfitpred = fitpred
for(i in 1:3) covariates(logfitpred[[i]])$D.0 = log(covariates(fitpred[[i]])$D.0)
quartz(h=10,w=5)
par(mfrow=c(3,1))
plot(logfitpred[[1]],asp=1,contour=FALSE) 
plot(logfitpred[[2]],asp=1,contour=FALSE) 
plot(logfitpred[[3]],asp=1,contour=FALSE) 

# Plot effect of stdGC on density:
masks = list(TostMask, NoyonMask,NemegtMask)
quartz(h=9,w=9)
par(mfrow=c(3,2))
for(i in 1:3) {
  Dhat = covariates(fitpred[[i]])$D.0
  ord = order(Dhat)
  stdGC4plot = covariates(masks[[i]])$stdGC[ord]
  plot(stdGC4plot,sort(Dhat),type="l")
  plot(masks[[i]],covariate="stdGC",contour=FALSE,asp=1)
}

