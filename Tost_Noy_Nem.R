library(secr)
library(fields)
library(maptools)
source("scrplotting.r")
setwd("C:/Users/koust/Dropbox (Snow Leopard Trust)/CREEM/Analyses/snowleopards")

# Testing reading traps as .csv or .txt (the former does not work)
#test = read.csv("./Tost_Noyon_Nemegt/Tost_Traps.csv");head(test)
#test1 = read.traps(file="./Tost_Noyon_Nemegt/Tost_Traps.txt", detector="count", covnames = c("Effort","Rgd","Topo", "Water"))
#test2 = read.traps(file="./Tost_Noyon_Nemegt/Noyon_Traps.txt", detector="count", covnames = c("Effort","Rgd","Topo", "Water"))
#test3 = read.traps(file="./Tost_Noyon_Nemegt/Nemegt_Traps.txt", detector="count", covnames = c("Effort","Rgd","Topo", "Water"))
#summary(test1)
#summary(test2)
#summary(test3)

# Read combined capture file and boundary for Tost, Noyon & Nemegt
TNN.trapfiles = c("./Tost_Noyon_Nemegt/Tost_Traps.txt",
                  "./Tost_Noyon_Nemegt/Noyon_Traps.txt",
                  "./Tost_Noyon_Nemegt/Nemegt_Traps.txt")
all.data.TNN<-read.capthist(captfile = "./Tost_Noyon_Nemegt/TNN_Capture.csv", 
                            binary.usage = FALSE, trapfile = TNN.trapfiles, 
                            detector="count", fmt = "trapID", 
                            trapcovnames = c("Rgd","Topo", "Water"))
summary(all.data.TNN)

head(TNN.trapfiles)

# Reduced Nemegt Data for 2 uncertain ids
all.data.TNN_R<-read.capthist(captfile = "./Tost_Noyon_Nemegt/TNN_Capture_R.csv", 
                            binary.usage = FALSE, trapfile = TNN.trapfiles, 
                            detector="count", fmt = "trapID", 
                            trapcovnames = c("Rgd","Topo", "Water"))
summary(all.data.TNN_R)

#binary.usage=FALSE to be inserted in the trap file read command. ID X Y USage / <covariates>
# old command has all traps in a single file, not by session:
#all.data.TNN<-read.capthist(captfile = "./Tost_Noyon_Nemegt/TNN_Capture.csv", trapfile = "./Tost_Noyon_Nemegt/TNN_Traps.csv", detector="count", fmt = "trapID", trapcovnames = c("Effort","Rgd","Topo", "Water"))
covariates(traps(all.data.TNN))
covariates(traps(all.data.TNN_R))

summary(all.data.TNN)

# Read boundary files
boundaryNemegt=readShapeSpatial("./Nemegt//Habitat/Nemegt_StudyArea.shp")
boundaryNoyon=readShapeSpatial("./Noyon2013/Habitat/NoyonStudy_Area.shp")
boundaryTost=readShapeSpatial("./Tost//Habitat/TostStudy_Area.shp")
# and plot them
plot(boundaryNemegt)
plot(boundaryNoyon)
plot(boundaryTost)

# plot all together:
# Find extent of bounding boxes of boundarys:
bbox.Nemegt = bbox(boundaryNemegt)
bbox.Noyon = bbox(boundaryNoyon)
bbox.Tost = bbox(boundaryTost)
bbxlim = range(bbox.Nemegt["x",],bbox.Noyon["x",],bbox.Tost["x",]) #Set plot limit for all 3 areas together
bbylim = range(bbox.Nemegt["y",],bbox.Noyon["y",],bbox.Tost["y",])

# redundant code for old traps without sessions
#trapsxy = traps(all.data.TNN)[[1]] # (strange that each element of list has ALL traps in it)
#trnames = row.names(trapsxy)
#sessnum = rep(NA,length(trnames))
#sessnum[which(substr(trnames,1,6)=="Nemegt")] = 3
#sessnum[which(substr(trnames,1,5)=="Noyon")] = 2
#sessnum[which(substr(trnames,1,4)=="Tost")] = 1
plot(bbxlim,bbylim,xlim=bbxlim,ylim=bbylim,type="n") 
plot(boundaryNemegt,add=TRUE,border=3)
plot(boundaryNoyon,add=TRUE,border=2)
plot(boundaryTost,add=TRUE,border=1)


# Make 3 masks
NemegtMask=make.mask(traps(all.data.Nemegt), spacing=500, buffer = 25000, type="trapbuffer", poly=boundaryNemegt)

NoyonMask=make.mask(traps(all.data.Noyon), spacing=500, buffer = 25000, type="trapbuffer", poly=boundaryNoyon)

TostMask=make.mask(traps(all.data.Tost), spacing=500, buffer = 25000, type="trapbuffer", poly=boundaryTost)

summary(NemegtMask) #To get the total area of the mask
summary(NoyonMask)
summary(TostMask)

# Read ruggedness covariate and put it into mask covariate GRIDCODE
SLCost.Nemegt<-readShapePoly("./Nemegt/Habitat/Nemegt_Rgd500m.shp")  #ruggedness pixels averaged over 500m radius
NemegtMask1<-addCovariates(NemegtMask, SLCost.Nemegt)

SLCost.Noyon<-readShapePoly("./Noyon2013/Habitat/Noyon_Rgd500m.shp")  #ruggedness pixels averaged over 500m radius
NoyonMask1<-addCovariates(NoyonMask, SLCost.Noyon)

SLCost.Tost<-readShapePoly("./Tost/Habitat/Tost_Rgd500m.shp")  #ruggedness pixels averaged over 500m radius
TostMask1<-addCovariates(TostMask, SLCost.Tost)


# Read binary habitat suitability code into mask covariate ...
SLCostBINARY.Nemegt<-readShapePoly("./Tost//Habitat/tost_sl.shp")  #Logistic binary SL habitat created using telemetry data
NemegtMask1<-addCovariates(NemegtMask1, SLCostBINARY.Nemegt)
names(covariates(NemegtMask1))[3:4] = c("binaryID","BINCODE") #Rename headers


SLCostBINARY.Noyon<-readShapePoly("./Tost//Habitat/tost_sl.shp")  #Logistic binary SL habitat created using telemetry data
NoyonMask1<-addCovariates(NoyonMask1, SLCostBINARY.Noyon)
names(covariates(NoyonMask1))[3:4] = c("binaryID","BINCODE") #Rename headers

SLCostBINARY.Tost<-readShapePoly("./Tost//Habitat/tost_sl.shp")  #ruggedness pixels averaged over 500m radius
TostMask1<-addCovariates(TostMask1, SLCostBINARY.Tost)
names(covariates(TostMask1))[3:4] = c("binaryID","BINCODE") #Rename headers

# make NAs in BINCODE zeros:
covariates(NemegtMask1)$BINCODE[is.na(covariates(NemegtMask1)$BINCODE)] = 0
summary(covariates(NemegtMask1))

covariates(NoyonMask1)$BINCODE[is.na(covariates(NoyonMask1)$BINCODE)] = 0
summary(covariates(NoyonMask1))

covariates(TostMask1)$BINCODE[is.na(covariates(TostMask1)$BINCODE)] = 0
summary(covariates(TostMask1))

# Standarize Rgd on traps (this makes fits a bit more stable)
# -----------------------------------------------------------
TNN.cams=traps(all.data.TNN)
# look at covariates for each session:
lapply(covariates(TNN.cams),summary)
# To standardise in same way over all sessions, need to combine, standarise and then separate:
Rgds = c(covariates(TNN.cams[[1]])$Rgd,
         covariates(TNN.cams[[2]])$Rgd,
         covariates(TNN.cams[[3]])$Rgd)
stdRgd = scale(Rgds)
n1 = dim(TNN.cams[[1]])[1]
n2 = dim(TNN.cams[[2]])[1]
n3 = dim(TNN.cams[[3]])[1]
covariates(TNN.cams[[1]])$stdRgd = stdRgd[1:n1]
covariates(TNN.cams[[2]])$stdRgd = stdRgd[(n1+1):(n1+n2)]
covariates(TNN.cams[[3]])$stdRgd = stdRgd[(n1+n2+1):(n1+n2+n3)]

# Do for revised dataset
TNN.camsR=traps(all.data.TNN_R)
# look at covariates for each session:
lapply(covariates(TNN.camsR),summary)
# To standardise in same way over all sessions, need to combine, standarise and then separate:
RgdsR = c(covariates(TNN.camsR[[1]])$Rgd,
         covariates(TNN.camsR[[2]])$Rgd,
         covariates(TNN.camsR[[3]])$Rgd)
stdRgdR = scale(RgdsR)
n1 = dim(TNN.camsR[[1]])[1]
n2 = dim(TNN.camsR[[2]])[1]
n3 = dim(TNN.camsR[[3]])[1]
covariates(TNN.camsR[[1]])$stdRgd = stdRgdR[1:n1]
covariates(TNN.camsR[[2]])$stdRgd = stdRgdR[(n1+1):(n1+n2)]
covariates(TNN.camsR[[3]])$stdRgd = stdRgdR[(n1+n2+1):(n1+n2+n3)]


# look at covariates for each session again:
lapply(covariates(TNN.cams),summary)
lapply(covariates(TNN.camsR),summary)

# put traps back into capthist (if don't put back by list elements, all.data.TNN becomes a 'traps' object!)
traps(all.data.TNN[[1]]) = TNN.cams[[1]]
traps(all.data.TNN[[2]]) = TNN.cams[[2]]
traps(all.data.TNN[[3]]) = TNN.cams[[3]]

traps(all.data.TNN_R[[1]]) = TNN.camsR[[1]]
traps(all.data.TNN_R[[2]]) = TNN.camsR[[2]]
traps(all.data.TNN_R[[3]]) = TNN.camsR[[3]]
covariates(traps(all.data.TNN_R))

# Standarize GRIDCODE (in stdGC) and BINCODE (in stdBC) on mask
# ------------------------------------------------------------------------
summary(covariates(NemegtMask1))
covariates(NemegtMask1)$stdGC = scale(covariates(NemegtMask1)$GRIDCODE)
covariates(NemegtMask1)$stdBC = scale(covariates(NemegtMask1)$BINCODE)
summary(covariates(NemegtMask1))
names(covariates(NemegtMask1))

summary(covariates(NoyonMask1))
covariates(NoyonMask1)$stdGC = scale(covariates(NoyonMask1)$GRIDCODE)
covariates(NoyonMask1)$stdBC = scale(covariates(NoyonMask1)$BINCODE)
summary(covariates(NoyonMask1))
names(covariates(NoyonMask1))

summary(covariates(TostMask1))
covariates(TostMask1)$stdGC = scale(covariates(TostMask1)$GRIDCODE)
covariates(TostMask1)$stdBC = scale(covariates(TostMask1)$BINCODE)
summary(covariates(TostMask1))
names(covariates(TostMask1))

# Plot sessions together:
pdf("Allregions.pdf",h=5,w=10)
plot(bbxlim,bbylim,xlim=xlim,ylim=ylim,xlab="",ylab="",bty="n",type="n",xaxt="n",yaxt="n",asp=1) #Generates Error!
# Plot the terrain
plot(NoyonMask1, covariate="stdGC", contour = FALSE, col = terrain.colors(16), legend = FALSE, add = TRUE)
plot(NemegtMask1, covariate="stdGC", contour=FALSE, col=terrain.colors(16), legend = FALSE, add = TRUE)
plot(TostMask1, covariate="stdGC", contour=FALSE, col=terrain.colors(16), legend = FALSE, add=TRUE)
# Add the traps
plot(traps(all.data.TNN)[[1]],add=TRUE,detpar=list(col="black",pch="+"))
plot(traps(all.data.TNN)[[2]],add=TRUE,detpar=list(col="black",pch="+"))
plot(traps(all.data.TNN)[[3]],add=TRUE,detpar=list(col="black",pch="+"))
# and the borders over the top
plot(boundaryNemegt,add=TRUE)#,border=3)
plot(boundaryNoyon,add=TRUE)#,border=2)
plot(boundaryTost,add=TRUE)#,border=1)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
dev.off()

TNN.hhn<-secr.fit(all.data.TNN, 
                  model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", 
                  mask=list(TostMask1, NoyonMask1, NemegtMask1))

TNN.hhnR<-secr.fit(all.data.TNN_R, 
                  model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", 
                  mask=list(TostMask1, NoyonMask1, NemegtMask1))

coefficients(TNN.hhn)
coefficients(TNN.hhnR)


####Is the code below correct? Can't seem to add different covariates for different sessions (here areas)  
TNN.hhn.DRgd<-secr.fit(all.data.TNN, model = list(D~stdGC, lambda0~1, sigma~1), detectfn="HHN",
                       mask = list(TostMask1, NoyonMask1, NemegtMask1))

TNN.hhn.DRgd.DetTopo10<-secr.fit(all.data.TNN, model = list(D~stdGC, lambda0~Topo, sigma~1), detectfn="HHN",
                       mask = list(TostMask1, NoyonMask1, NemegtMask1))


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
                       mask = list(TostMask1, NoyonMask1, NemegtMask1),sessioncov=data.frame(sfac=sess))

TNN.hhn.DRgd.sessR<-secr.fit(all.data.TNN_R, model = list(D~stdGC+sfac, lambda0~1, sigma~1), detectfn="HHN",
                            mask = list(TostMask1, NoyonMask1, NemegtMask1),sessioncov=data.frame(sfac=sess))

TNN.hhn.DRgd.sess.DetW<-secr.fit(all.data.TNN, model = list(D~stdGC+sfac, lambda0~sfac*Water, sigma~1), detectfn="HHN",
                            mask = list(TostMask1, NoyonMask1, NemegtMask1),sessioncov=data.frame(sfac=sess))

TNN.hhn.DRgd.sess.DetWR<-secr.fit(all.data.TNN_R, model = list(D~stdGC+sfac, lambda0~sfac*Water, sigma~1), detectfn="HHN",
                                 mask = list(TostMask1, NoyonMask1, NemegtMask1),sessioncov=data.frame(sfac=sess))

TNN.hhn.DRgd.DetTopo10W<-secr.fit(all.data.TNN, model = list(D~stdGC, lambda0~Topo+Water, sigma~1), detectfn="HHN",
                                mask = list(TostMask1, NoyonMask1, NemegtMask1))
TNN.hhn.DRgd.DetTopo10WR<-secr.fit(all.data.TNN_R, model = list(D~stdGC, lambda0~Topo+Water, sigma~1), detectfn="HHN",
                                  mask = list(TostMask1, NoyonMask1, NemegtMask1))

TNN.hhn.DRgd.sess_interact<-secr.fit(all.data.TNN, model = list(D~stdGC*sfac, lambda0~1, sigma~1), detectfn="HHN",
                            mask = list(TostMask1, NoyonMask1, NemegtMask1),sessioncov=data.frame(sfac=sess))


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
TNN.hhn.DHab.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                            model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                            details = list(userdist = userdfn1),
                               start = list(noneuc = 1)) #-1 gets rid of the intercept
TNN.hhn.DHab.nonUR<-secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                            model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                            details = list(userdist = userdfn1),
                            start = list(noneuc = 1)) #-1 gets rid of the intercept

# Model with stdGC in noneuc, Topo in Detect:
# ---------------------------
TNN.hhn.DHab.DetTopo10.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                            model=list(D~stdGC, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                            details = list(userdist = userdfn1),
                            start = list(noneuc = 1)) #-1 gets rid of the intercept

TNN.hhn.DHab.DetTopo10.nonU1<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                                      model=list(D~stdGC, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                                      details = list(userdist = userdfn1),
                                      start = list(noneuc = 1)) #-1 gets rid of the intercept


TNN.hhn.DHab.DetTopo10.nonUR<-secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                                      model=list(D~stdGC, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                                      details = list(userdist = userdfn1),
                                      start = list(noneuc = 1)) #-1 gets rid of the intercept

# Model with BC in noneuc, Topo in Detect, GC in density:
# ---------------------------

TNN.hhn.DGC.DetTopo10.nonUGBR<-secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                                       model=list(D~stdGC, lambda0~Topo, sigma~1, noneuc ~ stdBC -1), 
                                       details = list(userdist = userdfn1),
                                       start = list(noneuc = 1)) #-1 gets rid of the intercept

# Model with GC in noneuc, Topo in Detect, GB in density:
# ---------------------------

TNN.hhn.DGB.DetTopo10.nonUR<-secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
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
TNN.hhn.DHab.DetToposess.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                                      model=list(D~stdGC, lambda0~Topo*sfac, sigma~1, noneuc ~ stdGC -1), 
                                      details = list(userdist = userdfn1),sessioncov=data.frame(sfac=sess),
                                      start = list(noneuc = 1)) #-1 gets rid of the intercept
TNN.hhn.DHab.DetToposess.nonUR<-secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                                        model=list(D~stdGC, lambda0~Topo*sfac, sigma~1, noneuc ~ stdGC -1), 
                                        details = list(userdist = userdfn1),sessioncov=data.frame(sfac=sess),
                                        start = list(noneuc = 1)) #-1 gets rid of the intercept


coefficients(TNN.hhn.DHab.DetToposess.nonU)


# Model with D->Rgd session interact, noneuc->rgd, Topo in lam0:
# ---------------------------
TNN.hhn.DHab.S.DetTopo10.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                                      model=list(D~stdGC*sfac, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                                      details = list(userdist = userdfn1),sessioncov=data.frame(sfac=sess),
                                      start = list(noneuc = 1)) #-1 gets rid of the intercept

TNN.hhn.DHab.S.DetTopo10.nonUR<-secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                                        model=list(D~stdGC*sfac, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                                        details = list(userdist = userdfn1),sessioncov=data.frame(sfac=sess),
                                        start = list(noneuc = 1)) #-1 gets rid of the intercept

# Model with D->Rgd session add, noneuc->rgd, Topo in lam0:
# ---------------------------
TNN.hhn.DHab_S.DetTopo10.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                                        model=list(D~stdGC+sfac, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                                        details = list(userdist = userdfn1),sessioncov=data.frame(sfac=sess),
                                        start = list(noneuc = 1)) #-1 gets rid of the intercept

TNN.hhn.DHab_S.DetTopo10.nonUR<-secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                                        model=list(D~stdGC+sfac, lambda0~Topo, sigma~1, noneuc ~ stdGC -1), 
                                        details = list(userdist = userdfn1),sessioncov=data.frame(sfac=sess),
                                        start = list(noneuc = 1)) #-1 gets rid of the intercept

coefficients(TNN.hhn.DHab.DetTopo10.nonU)
# Model with stdGC in noneuc, Topo+Water in Detect:
# ---------------------------
TNN.hhn.DHab.LamTopoWat.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                                      model=list(D~stdGC, lambda0~Topo+Water, sigma~1, noneuc ~ stdGC -1), 
                                      details = list(userdist = userdfn1),
                                      start = list(noneuc = 1)) #-1 gets rid of the intercept

TNN.hhn.DHab.LamTopoWat.nonUR<-secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                                       model=list(D~stdGC, lambda0~Topo+Water, sigma~1, noneuc ~ stdGC -1), 
                                       details = list(userdist = userdfn1),
                                       start = list(noneuc = 1)) #-1 gets rid of the intercept

coefficients(TNN.hhn.DHab.LamTopoWat.nonU)

# Model with constt D in noneuc:
# ---------------------------
TNN.hhn.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                            model=list(D~1, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                            details = list(userdist = userdfn1), start = list(noneuc = 1)) #-1 gets rid of the intercept

TNN.hhn.nonUR<-secr.fit(all.data.TNN_R, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
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
TNN.hhn.DRgd_sess.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1),
                                 model=list(D~stdGC+sfac, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                                 sessioncov=data.frame(sfac=sess),details = list(userdist = userdfn1), 
                                 start = list(noneuc = 1)) #-1 gets rid of the intercept


# Model with D dependent on Rgd & different for each session in noneuc:
# ---------------------------
TNN.hhn.DRgd.sess.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1),
                                 model=list(D~stdGC*sfac, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                                 sessioncov=data.frame(sfac=sess),details = list(userdist = userdfn1), 
                                 start = list(noneuc = 1)) #-1 gets rid of the intercept

# Model with D dependent on Rgd &  different for each session, lambda on water in noneuc:
# ---------------------------
TNN.hhn.DRgd.sess.D.W_sess.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1),
                                          model=list(D~stdGC*sfac, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                                          sessioncov=data.frame(sfac=sess), details = list(userdist = userdfn1), 
                                          start = list(noneuc = 1)) #-1 gets rid of the intercept

# Model with D dependent on Rgd, lambda on water & BOTH different for each session in noneuc:
# ---------------------------
TNN.hhn.DRgd.sess.D.W.sess.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                            model=list(D~stdGC*sfac, lambda0~Water*sfac, sigma~1, noneuc ~ stdGC -1), 
                            sessioncov=data.frame(sfac=sess),
                            details = list(userdist = userdfn1), start = list(noneuc = 1)) #-1 gets rid of the intercept

# Model with D dependent on Rgd, lambda on water+topo & BOTH different for each session in noneuc:
# ---------------------------
TNN.hhn.DRgd.sess.D_W_sess.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
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
region.N(TNN.hhn.DRgd.sess,region=TostMask1,session="1")
region.N(TNN.hhn.DRgd.sess,region=NoyonMask1,session="2")
region.N(TNN.hhn.DRgd.sess,region=NemegtMask1,session="3")

region.N(TNN.hhn.DRgd.sess,region=TostMask1,session="1")
region.N(TNN.hhn.DRgd.sess,region=NoyonMask1,session="2")
region.N(TNN.hhn.DRgd.sess,region=NemegtMask1,session="3")

region.N(TNN.hhn.DRgd.sess.D.W.sess.nonU,region=TostMask1,session="1")
region.N(TNN.hhn.DRgd.sess.D.W.sess.nonU,region=NoyonMask1,session="2")
region.N(TNN.hhn.DRgd.sess.D.W.sess.nonU, region=NemegtMask1,session="3")

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
TNN.hhn.DHabS3.lambdaSfacWater.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                              model=list(D~s(stdGC,k=3), lambda0~sfac*Water, sigma~1, noneuc ~ stdGC -1), 
                              sessioncov=data.frame(sfac=sess), details = list(userdist = userdfn1),
                              start = list(noneuc = 1)) #-1 gets rid of the intercept
# Some convergence problems with the above model; need to refit starting from this fit, but not yet done that.

# David trial Model with stdGC in noneuc and smooth (k=3) stdGC effect on D:
# ---------------------------
TNN.hhn.DHabS3.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                              model=list(D~s(stdGC,k=3), lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                              details = list(userdist = userdfn1),
                              start = list(noneuc = 1)) #-1 gets rid of the intercept
# Fit again, using estimates from above as starting values:
TNN.hhn.DHabS3a.nonU<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
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
region.N(fit,region=TostMask1,session="1")
region.N(fit,region=NoyonMask1,session="2")
region.N(fit,region=NemegtMask1,session="3")

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
masks = list(TostMask1, NoyonMask1,NemegtMask1)
quartz(h=9,w=9)
par(mfrow=c(3,2))
for(i in 1:3) {
  Dhat = covariates(fitpred[[i]])$D.0
  ord = order(Dhat)
  stdGC4plot = covariates(masks[[i]])$stdGC[ord]
  plot(stdGC4plot,sort(Dhat),type="l")
  plot(masks[[i]],covariate="stdGC",contour=FALSE,asp=1)
}

