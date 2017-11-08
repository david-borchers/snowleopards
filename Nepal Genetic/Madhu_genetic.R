library(secr)
library(fields)
library(maptools)
source("scrplotting.r")
library(gdistance)

Trial.Data1<-read.capthist(captfile = "./Nepal Genetic/Madhu_Caps2.txt",
                           trapfile ="./Nepal Genetic/Madhu_T2.txt", detector="transect",fmt = "XY")

Madhu.Data1<-read.capthist(captfile = "./Nepal Genetic/Madhu_Caps_Fin.txt", binary.usage=FALSE,
                           trapfile ="./Nepal Genetic/Madhu_Traps_Fin_NoCov.txt", detector="transect",
                           fmt = "XY")

Madhu.Data2<-read.capthist(captfile = "./Nepal Genetic/Madhu_Caps_Fin.txt", binary.usage=FALSE,
                           trapfile ="./Nepal Genetic/Madhu_Traps_Fin.txt", detector="transect",fmt = "XY",
                           trapcovnames = c("Topography","Habitats","Altitude"))

boundaryMadhu=readShapeSpatial("C:/Users/Koustubh/Dropbox (Snow Leopard Trust)/CREEM/Nepal/Madhu_poly1.shp")
boundaryMadhu2=readShapeSpatial("C:/Users/Koustubh/Dropbox (Snow Leopard Trust)/CREEM/Nepal/Madhu_poly_clip.shp")

plot(boundaryMadhu)
plot(boundaryMadhu2)
plot(traps(Madhu.Data1), add=TRUE)
plot(Madhu.Data1, add=TRUE)
traps(Madhu.Data2)

# Make mask:
MadhuMask=make.mask(traps(Madhu.Data1), spacing=1000, type="traprect", buffer = 30000, 
                    poly=boundaryMadhu)
MadhuMask2=make.mask(traps(Madhu.Data2), spacing=1000, buffer = 35000, type="traprect", 
                     poly=boundaryMadhu2) #Creates a rectangular buffer around detectors, and clips to polygon
#boundaryMadhu2 excludes habitat below and above SL range

?make.mask

plot(MadhuMask2)
RPSV(Madhu.Data2,CC=TRUE)
summary(Madhu.Data2)

Madhu.hhn<-secr.fit(Madhu.Data2, model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", mask=MadhuMask)
region.N(Madhu.hhn)

summary(Madhu.Data2)
# Standarize Rgd on traps (this makes fits a bit more stable)
# -----------------------------------------------------------
summary(covariates(traps(Madhu.Data2)))
covariates(traps(Madhu.Data2))$stdAlt = scale(covariates(traps(Madhu.Data2))$Altitude)
summary(covariates(traps(Madhu.Data2)))
head(covariates(traps(Madhu.Data2)))

# Read ruggedness covariate and put it into mask covariate GRIDCODE
SLRgd.Madhu<-readShapePoly("C:/Users/Koustubh/Dropbox (Snow Leopard Trust)/CREEM/Nepal/Rugged_1000shp.shp")  #ruggedness pixels averaged over 500m radius
SLAlt.Madhu<-readShapePoly("C:/Users/Koustubh/Dropbox (Snow Leopard Trust)/CREEM/Nepal/Altitude.shp")
# using focal statistics tool

MadhuMask2a<-addCovariates(MadhuMask2, SLRgd.Madhu)
MadhuMask2a<-addCovariates(MadhuMask2a, SLAlt.Madhu)
names(covariates(MadhuMask2a))[3:4] = c("Rgd","Alt")

summary(MadhuMask2a)
plotcovariate(MadhuMask2a, covariate = "GRIDCODE", asp =1, contour = FALSE)


Madhu.hhn.lam_Topo<-secr.fit(Madhu.Data2, model=list(D~1, lambda0~"Topography", sigma~1), 
                             detectfn="HHN", mask=MadhuMask)

Madhu.hhn.lam_Alt<-secr.fit(Madhu.Data2, model=list(D~1, lambda0~"stdAlt", sigma~1), 
                            detectfn="HHN", mask=MadhuMask)

Madhu.hhn.lam_Habitat<-secr.fit(Madhu.Data2, model=list(D~1, lambda0~"Topography", sigma~1), 
                                detectfn="HHN", mask=MadhuMask)

skMadhu.hhn.D_alt<-secr.fit(Madhu.Data2, model=list(D~"StdAltitude", lambda0~"Topography", sigma~1), 
                            detectfn="HHN", mask=MadhuMask)


all.data.Zoloon<-read.capthist(captfile = "./Nepal Genetic/Zoloon_caps.csv", 
                               trapfile = "./Nepal Genetic/Zoloon_traps.csv", binary.usage=FALSE,
                               detector="count", fmt = "trapID", 
                               trapcovnames = c("Rgd",	"Topo",	"Altidute", "Water",	"Level"))
ZoloonBoundary=readShapeSpatial("./Nepal Genetic/Zoolon_UTM1.shp")
plot(ZoloonBoundary)
plot(x=all.data.Zoloon, add = TRUE, legend = TRUE)
ZoloonMask=make.mask(traps(all.data.Zoloon), spacing=500, buffer = 25000, type="trapbuffer", poly=ZoloonBoundary)
summary(all.data.Zoloon)

plot(ZoloonMask)
Zoloon.hhn<-secr.fit(all.data.Zoloon, model=list(D~1, lambda0~1, sigma~1), 
                     detectfn="HHN", mask=ZoloonMask)
derived(Zoloon.hhn)
region.N(Zoloon.hhn)
