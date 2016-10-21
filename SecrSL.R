#SEcR Noyon vs Tost & Sarychat: comparing behavior in protected & unprotected areas


#Running basic SECR models with Sl data
# Going ahead with Gobi and Sarychat. Data files and shapes.
# David added this comment
library(secr)
library(spatstat)
library(maptools)
library(fields)
library(rgdal)
library(secrlinear)
library(gdistance)
library(raster)
library(rgeos)
gpclibPermit()

#SECR models for Noyon
all.data.Noyon<-read.capthist(captfile = "./Noyon2013/Noyon_capthist2013secr.csv", trapfile = "./Noyon2013/Noyon_trap2013secr.csv", detector="count", fmt = "trapID", trapcovnames = c("Effort", "topo", "substrate", "brokenness", "Rgd"))
Noyonboundary=readShapeSpatial("./Noyon2013/Habitat/NoyonStudy_Area.shp")
plot(Noyonboundary)
plot(x=all.data.Noyon, add=TRUE)
NoyonMask=make.mask(traps(all.data.Noyon), spacing=500, buffer = 25000, type="trapbuffer", poly=Noyonboundary)
SLCost.Noyon<-readShapePoly("./Noyon2013/Habitat/Noyon_Rgd500m.shp") #ruggedness pixels averaged over 500m radius
head(SLCost.Noyon)
head(NoyonMask)

NoyonMask1<-addCovariates(NoyonMask, SLCost.Noyon)
head(covariates(NoyonMask1))
summary(covariates(NoyonMask1))
summary(covariates(traps(all.data.Noyon)))

plot(NoyonMask1, covariate="stdGRIDCODE", contour=FALSE, col=terrain.colors(15), legend = FALSE)

# Standarize Rgd (this makes fits a bit more stable)
# --------------------------------------------------
summary(covariates(traps(all.data.Noyon)))
covariates(traps(all.data.Noyon))$stdRgd = scale(covariates(traps(all.data.Noyon))$Rgd)
summary(covariates(traps(all.data.Noyon)))
head(covariates(traps(all.data.Noyon)))

# Also standarize stdGRIDCODE for good measure
# -----------------------------------------
summary(covariates(NoyonMask1))
covariates(NoyonMask1)$stdGC = scale(covariates(NoyonMask1)$stdGRIDCODE)
names(covariates(NoyonMask1))



Noyon.hhn<-secr.fit(all.data.Noyon, model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", mask=NoyonMask1)
Noyon.hhnxy<-secr.fit(all.data.Noyon, model=list(D~(x+y), lambda0~1, sigma~1), detectfn="HHN", mask=NoyonMask1)
Noyon.hhn.Det_Rgd01<-secr.fit(all.data.Noyon, model = list(D~1, lambda0~1, sigma~stdRgd), detectfn="HHN", mask = NoyonMask1)
Noyon.hhn.Det_Rgd10<-secr.fit(all.data.Noyon, model = list(D~1, lambda0~stdRgd, sigma~1), detectfn="HHN", mask = NoyonMask1)

# convergence problem with this model
Noyon.hhn.Det_Rgd11<-secr.fit(all.data.Noyon, model = list(D~1, lambda0~stdRgd, sigma~stdRgd), detectfn="HHN", mask = NoyonMask1)

# this one fits OK
counts.hhn.D.Rgd<-secr.fit(all.data.Noyon, model = list(D~stdGRIDCODE, lambda0~1, sigma~stdRgd), detectfn="HHN", mask = NoyonMask1)
# sigma.stdRgd=coefficients(counts.hhn.D.Rgd)["sigma.stdRgd",1] = 0.303697679, 
# which suggests that each additional unit of stdRgd INcreases sigma by a 
# mutiplicative factor of exp(0.303697679), i.e. about 35%
# and this wee calculation:
abs(diff(range(covariates(traps(all.data.Noyon))$stdRgd)))*exp(coefficients(counts.hhn.D.Rgd)["sigma.stdRgd",1])
# says that sigma is about 5.8 times as large at the biggest stdRgd value than the smallest!

# this one fits OK
Noyon.hhn.D.Rgd2<-secr.fit(all.data.Noyon, model = list(D~stdGRIDCODE, lambda0~1, sigma~1), detectfn="HHN", mask = NoyonMask1)

# this one fits OK
counts.hhn.D.Rgd3<-secr.fit(all.data.Noyon, model = list(D~s(stdGRIDCODE, k=3), lambda0~1, sigma~1), detectfn="HHN", mask = NoyonMask1)

AIC(counts.hhn.D.Rgd,Noyon.hhn.D.Rgd2,counts.hhn.D.Rgd3)

# David has not done anything beyond here yet
# ===========================================


AIC(Noyon.hhn, Noyon.hhnxy, Noyon.hhn.Det_Rgd01, Noyon.hhn.D.Rgd2)
coefficients(counts.hhnRgd)
coefficients(counts.hhn.D.Rgd)

NoySurface<-predictDsurface(Noyon.hhn.D.Rgd2, se.D=TRUE, cl.D=TRUE)
plot(NoySurface)

NoySurface2<-predictDsurface(counts.hhnRgd, se.D=TRUE, cl.D=TRUE)

plot(NoySurface2)

NoyonDsurface<-attr(NoySurface, "covariates")$D.0 #Create a surface of Density
NoyonSEsurface<-attr(NoySurface, "covariates")$SE.0 #Create a surface of Standard Errors
attr(NoySurface, "covariates")$CV<-NoyonSEsurface/NoyonDsurface #Create a surface of CVs
plot(NoySurface, covariate="CV", breaks = c(1:10)/10)
points(traps(all.data))


#Running SECR for Tost 2012
setwd("C:/Users/Koustubh/Dropbox (Snow Leopard Trust)/CREEM/Tost")
all.data.Tost<-read.capthist(captfile = "Tost_capthist2012.csv", trapfile = "Tost_cams_rugged2012.csv", detector="count", fmt = "trapID", trapcovnames = c("Effort",	"Topo",	"Altidute",	"Rgd"))
boundaryTost=readShapeSpatial("C:/Users/Koustubh/Dropbox (Snow Leopard Trust)/CREEM/Tost/Habitat/TostStudy_Area.shp")

plot(boundaryTost)
plot(x=all.data.Tost, add=TRUE)
TostMask=make.mask(traps(all.data.Tost), spacing=500, buffer = 25000, type="trapbuffer", poly=boundaryTost)
SLCost.Tost<-readShapePoly("C:/Users/Koustubh/Dropbox (Snow Leopard Trust)/CREEM/Tost/Habitat/Tost_Rgd500m.shp")  #ruggedness pixels averaged over 500m radius

head(SLCost.Tost)
head(TostMask)

TostMask1<-addCovariates(TostMask, SLCost.Tost)
head(covariates(TostMask1))
summary(covariates(TostMask1))
summary(covariates(traps(all.data.Tost)))

plot(TostMask1, covariate="stdGRIDCODE", contour=FALSE, col=terrain.colors(10), legend = FALSE)
head(covariates(traps(all.data.Tost)))
head(covariates(TostMask1))

Tost.hhn<-secr.fit(all.data.Tost, model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", mask=TostMask1)
Tost.hhn.detrug<-secr.fit(all.data.Tost, model=list(D~1, lambda0~Rgd, sigma~Rgd), detectfn="HHN", mask=TostMask1)
Tost.hhn.Dxy<-secr.fit(all.data.Tost, model=list(D~x+y, lambda0~1, sigma~1), detectfn="HHN", mask=TostM)
Tost.hhn.DHab<-secr.fit(all.data.Tost, model=list(D~stdGRIDCODE, lambda0~1, sigma~1), detectfn="HHN", mask=TostMask1)
Tost.hhn.Dx<-secr.fit(all.data.Tost, model=list(D~x, lambda0~1, sigma~1), detectfn="HHN", mask=TostMask1)

AIC(Tost.hhn.Dx, Tost.hhn.DHab, Tost.hhn.Dxy)
coefficients(Tost.hhn.DHab)
TostSurface<-predictDsurface(Tost.hhn.DHab, se.D=TRUE, cl.D=TRUE)
plot(TostSurface)

TostSurface<-predictDsurface(Tost.hhn.Dx, se.D=TRUE, cl.D=TRUE)

Nhat1<-region.N(Tost.hhn.DHab) #Estimates the population N of the animals within the region defined by mask


#SECR models for Nemegt
setwd("C:/Users/Koustubh/Dropbox (Snow Leopard Trust)/CREEM/Nemegt")
all.data.Nemegt<-read.capthist(captfile = "Nemegt2013_Capture.csv", trapfile = "Nemegt2013_Cams.csv", detector="count", fmt = "trapID", trapcovnames = c("Effort",	"Topo",	"Brokenness",	"Grass",	"Rgd"))
Nemegtboundary=readShapeSpatial("C:/Users/Koustubh/Dropbox (Snow Leopard Trust)/CREEM/Nemegt/Habitat/Nemegt_StudyArea.shp")
plot(Nemegtboundary)
plot(x=all.data.Nemegt, add=TRUE)
NemegtMask=make.mask(traps(all.data.Nemegt), spacing=500, buffer = 25000, type="trapbuffer", poly = Nemegtboundary)
SLHabitat.Nemegt<-readShapePoly("C:/Users/Koustubh/Dropbox (Snow Leopard Trust)/CREEM/Nemegt/Habitat/Nemegt_Rgd500m.shp")
head(SLHabitat.Nemegt)
head(NemegtMask)

NemegtMask1<-addCovariates(NemegtMask, SLHabitat.Nemegt)
head(covariates(NemegtMask1))
summary(covariates(NemegtMask1))
summary(covariates(traps(all.data.Nemegt)))

plot(NemegtMask1, covariate="stdGRIDCODE", contour=FALSE, col=terrain.colors(10), legend = FALSE)
head(covariates(traps(all.data.Nemegt)))

Nemegt.hhn<-secr.fit(all.data.Nemegt, model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", mask=NemegtMask1)
Nemegt.hhn.DRgd1<-secr.fit(all.data.Nemegt, model=list(D~stdGRIDCODE, lambda0~1, sigma~1), detectfn="HHN", mask=NemegtMask1)
Nemegt.hhn.Det_Topo<-secr.fit(all.data.Nemegt, model=list(D~stdGRIDCODE, lambda0~1, sigma~Topo), detectfn="HHN", mask=NemegtMask1)
Nemegt.hhn.Det_Topo1<-secr.fit(all.data.Nemegt, model=list(D~1, lambda0~1, sigma~Topo), detectfn="HHN", mask=NemegtMask1)
Nemegt.hhn.Det_Topo10<-secr.fit(all.data.Nemegt, model=list(D~1, lambda0~Topo, sigma~1), detectfn="HHN", mask=NemegtMask1)
Nemegt.hhn.DRgd.Det_Topo<-secr.fit(all.data.Nemegt, model=list(D~stdGRIDCODE, lambda0~Topo, sigma~Topo), detectfn="HHN", mask=NemegtMask1)
Nemegt.DRgd.Det_Topo<-Nemegt.hhn.DRgd.Det_Topo
AIC(Nemegt.hhn, Nemegt.hhn.DRgd, Nemegt.hhn.DRgd1, Nemegt.hhn.Det_Topo,Nemegt.hhn.Det_Topo1, Nemegt.DRgd.Det_Topo)

coefficients(Nemegt.hhn.DRgd.Det_Topo)
NemegtSurface<-predictDsurface(Nemegt.hhn.DRgd.Det_Topo, se.D=TRUE, cl.D=TRUE)
plot(NemegtSurface, legend=FALSE)
predict(Nemegt.hhn.DRgd.Det_Topo)

Nhat1<-region.N(Nemegt.hhn) #Estimates the population N of the animals within the region defined by mask
Nhat1
