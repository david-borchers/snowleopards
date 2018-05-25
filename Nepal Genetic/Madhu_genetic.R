library(secr)
library(fields)
library(maptools)
source("scrplotting.r")
library(gdistance)
library(raster)

Trial.Data1<-read.capthist(captfile = "./Nepal Genetic/Madhu_Caps2.txt",
                           trapfile ="./Nepal Genetic/Madhu_T2.txt", detector="transect",fmt = "XY")

Madhu.Data1<-read.capthist(captfile = "./Nepal Genetic/Madhu_Caps_Fin.txt", binary.usage=FALSE,
                           trapfile ="./Nepal Genetic/Madhu_Traps_Fin_NoCov.txt", detector="transect",
                           fmt = "XY")

Madhu.Data2<-read.capthist(captfile = "./Nepal Genetic/Madhu_Caps_Fin.txt", binary.usage=FALSE,
                           trapfile ="./Nepal Genetic/Madhu_Traps_Fin.txt", detector="transect",fmt = "XY",
                           trapcovnames = c("Topography","Habitats","Altitude"))

boundaryMadhu=readShapeSpatial("C:/Users/Koust/Dropbox (Snow Leopard Trust)/CREEM/Nepal/Madhu_poly1.shp")
boundaryMadhu2=readShapeSpatial("C:/Users/Koust/Dropbox (Snow Leopard Trust)/CREEM/Nepal/Madhu_poly_clip.shp")

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
write.csv(MadhuMask2, file = "./MadhuMask2.csv") #Save mask as CSV, to open and add covariates in ArcGIS. Much faster!
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
# -----------------------------------------------------------
#SLRgd.Madhu<-readShapePoly("C:/Users/Koust/Dropbox (Snow Leopard Trust)/CREEM/Nepal/Rugged_1000shp.shp")  #ruggedness pixels averaged over 500m radius

#SLAlt.Madhu<-readShapePoly("C:/Users/Koust/Dropbox (Snow Leopard Trust)/CREEM/Nepal/Altitude.shp")
# using focal statistics tool

#Create covariates in ArcGIS by adding columns to the exported CSV file
MadhuMask2a <- read.mask (file = "C:/Users/Koust/Dropbox (Snow Leopard Trust)/CREEM/Analyses/snowleopards/Nepal Genetic/Madhu_MaskCovs.csv", 
                        spacing = 1000, header = TRUE)

summary(covariates(MadhuMask2a))
covariates(MadhuMask2a)$StdRgd = scale(covariates(MadhuMask2a)$Rgd)
covariates(MadhuMask2a)$StdAlt = scale(covariates(MadhuMask2a)$Alt)
summary(covariates(MadhuMask2a))
names(covariates(MadhuMask2a))


plotcovariate(MadhuMask2a, covariate = "StdRgd")
plotcovariate(MadhuMask2a, covariate = "StdAlt")

summary(covariates(MadhuMask2a))


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


load("C:/Users/koust/Dropbox (Snow Leopard Trust)/CREEM/Analyses/snowleopards/Nepal Genetic/Madhu_FirstSet.RData")
AIC(Madhu.hhn, Madhu.hhn.lam_Topo, Madhu.hhn.lam_Alt, Madhu.hhn.lam_Habitat, skMadhu.hhn.D_alt, 
    skMadhu.hhn.D_alt.lam_Topo, skMadhu.hhn.D_rgd.lam_Topo, skMadhu.hhn.D_AltRgd.lam_Topo, skMadhu.hhn.D_Alt2.lam_Topo, criterion = "AIC")

coefficients(skMadhu.hhn.D_Alt2.lam_Topo)
MadhuTop<-predictDsurface(skMadhu.hhn.D_Alt2.lam_Topo)
getwd()

plot(MadhuTop,contour=FALSE)
head(covariates(MadhuTop))

write.csv(MadhuTop, file = "./MadhuPlot.csv") #Save top model surface as CSV, to open in ArcGIS!
write.csv(covariates(MadhuTop), file = "./MadhuPlot_Covs.csv") #Save top model surface as CSV, to open in ArcGIS!


model.average(skMadhu.hhn.D_Alt2.lam_Topo, skMadhu.hhn.D_alt.lam_Topo, Madhu.hhn.lam_Topo, skMadhu.hhn.D_AltRgd.lam_Topo,
              skMadhu.hhn.D_rgd.lam_Topo, Madhu.hhn.lam_Alt, realnames = "D")
predict(skMadhu.hhn.D_Alt2.lam_Topo)
coefficients(skMadhu.hhn.D_Alt2.lam_Topo)
coefficients(skMadhu.hhn.D_alt.lam_Topo)
coefficients(Madhu.hhn.lam_Topo)
coefficients(skMadhu.hhn.D_AltRgd.lam_Topo)
coefficients(skMadhu.hhn.D_rgd.lam_Topo)
coefficients(Madhu.hhn.lam_Alt)


ModAv_Pop<-(region.N(skMadhu.hhn.D_Alt2.lam_Topo)*0.33)+(region.N(skMadhu.hhn.D_alt.lam_Topo)*0.2599)+(region.N(Madhu.hhn.lam_Topo)*0.2246)+
                                                                                                         (region.N(skMadhu.hhn.D_AltRgd.lam_Topo)*.0977)+
                                                                                                         (region.N(skMadhu.hhn.D_rgd.lam_Topo)*0.0827)
ModAv_Pop
#Estimate max and min density from top model per 100 sq km
max(covariates(MadhuTop)$D.0)*10000
min(covariates(MadhuTop)$D.0)*10000
TopAbund<-region.N(skMadhu.hhn.D_Alt2.lam_Topo)
summary(MadhuMask)
D_ModAv<-ModAv_Pop*100/15233
D_ModAv



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
