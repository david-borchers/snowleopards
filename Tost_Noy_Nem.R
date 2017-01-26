library(secr)
library(fields)
library(maptools)
source("scrplotting.r")

# Read combined capture file and boundary for Tost, Noyon & Nemegt
all.data.TNN<-read.capthist(captfile = "./Tost_Noyon_Nemegt/TNN_Capture.csv", trapfile = "./Tost_Noyon_Nemegt/TNN_Traps.csv", detector="count", fmt = "trapID", trapcovnames = c("Effort","Rgd","Topo"))
covariates(traps(all.data.TNN))

boundaryNemegt=readShapeSpatial("./Nemegt//Habitat/Nemegt_StudyArea.shp")

plot(boundaryNemegt)

# Make 3 masks
NemegtMask=make.mask(traps(all.data.Nemegt), spacing=500, buffer = 25000, type="trapbuffer", poly=boundaryNemegt)

NoyonMask=make.mask(traps(all.data.Noyon), spacing=500, buffer = 25000, type="trapbuffer", poly=boundaryNoyon)

TostMask=make.mask(traps(all.data.Tost), spacing=500, buffer = 25000, type="trapbuffer", poly=boundaryTost)

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
summary(covariates(traps(all.data.TNN)))
names(covariates(traps(all.data.TNN)))
covariates(traps(all.data.TNN))$stdRgd = scale(covariates(traps(all.data.TNN))$Rgd)

summary(covariates(traps(all.data.TNN)))

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
plot(NoyonMask1, covariate="stdGC", contour = FALSE, col = terrain.colors(12), legend = FALSE, add = TRUE)
plot(NemegtMask1, covariate="stdGC", contour=FALSE, col=terrain.colors(16), legend = FALSE)
plot(TostMask1, covariate="stdGC", contour=FALSE, col=terrain.colors(16), legend = FALSE, add=TRUE)
###How to change xy of the screen to plot all masks together???

TNN.cams=traps(all.data.TNN)

TNN.hhn<-secr.fit(all.data.TNN, 
                  model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", 
                  mask=list(TostMask1, NoyonMask1, NemegtMask1))

####Is the code below correct? Can't seem to add different covariates for different sessions (here areas)  
TNN.hhn.DRgd<-secr.fit(all.data.TNN, model = list(D~stdGC, lambda0~1, sigma~1), detectfn="HHN",
                       mask = list(TostMask1, NoyonMask1, NemegtMask1))

coefficients(TNN.hhn)
predict(TNN.hhn)
TNNSurface.DRgd<-predictDsurface(TNN.hhn.DRgd)
windows()
plot(TNNSurface.DRgd,asp=1,contour=FALSE) #This generates an error. Something I am doing wrong here it seems...

# How to get region.N for each of the 3 areas?
 
