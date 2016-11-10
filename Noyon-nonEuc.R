library(secr)
library(fields)
library(maptools)
source("scrplotting.r")

#Running SECR for Noyon 2013

# Read capture file and boundary
all.data.Noyon<-read.capthist(captfile = "./Noyon2013/Noyon_capthist2013secr.csv", trapfile = "./Noyon2013/Noyon_trap2013secr.csv", detector="count", fmt = "trapID", trapcovnames = c("Effort",	"Topo",	"Substrate",	"Brokenness", "Rgd"))
boundaryNoyon=readShapeSpatial("./Noyon2013//Habitat/NoyonStudy_Area.shp")
# and plot it
plot(boundaryNoyon)
plot(x=all.data.Noyon, add=TRUE)

# Make mask:
NoyonMask=make.mask(traps(all.data.Noyon), spacing=500, buffer = 25000, type="trapbuffer", poly=boundaryNoyon)

# Read ruggedness covariate and put it into mask covariate GRIDCODE
SLCost.Noyon<-readShapePoly("./Noyon2013//Habitat/Noyon_Rgd500m.shp")  #ruggedness pixels averaged over 500m radius
NoyonMask1<-addCovariates(NoyonMask, SLCost.Noyon)

# Read binary habitat suitability code into mask covariate ...
SLCostBINARY.Noyon<-readShapePoly("./Noyon2013//Habitat/XXXX.shp")  #ruggedness pixels averaged over 500m radius
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