library(secr)
library(fields)
library(maptools)
source("scrplotting.r")
library(gdistance)

Trial.Data1<-read.capthist(captfile = "./Nepal Genetic/Madhu_Caps2.txt",
                          trapfile ="./Nepal Genetic/Madhu_T2.txt", detector="transect",fmt = "XY")

Madhu.Data1<-read.capthist(captfile = "./Nepal Genetic/Madhu_Caps_Fin.txt", binary.usage=FALSE,
                           trapfile ="./Nepal Genetic/Madhu_Traps_Fin_NoCov.txt", detector="transect",fmt = "XY",
                          trapcovnames = c("Topography","Habitats","Altitude"))

Madhu.Data2<-read.capthist(captfile = "./Nepal Genetic/Madhu_Caps_Fin.txt", binary.usage=FALSE,
                           trapfile ="./Nepal Genetic/Madhu_Traps_Fin.txt", detector="transect",fmt = "XY",
                           trapcovnames = c("Topography","Habitats","Altitude"))

boundaryMadhu=readShapeSpatial("C:/Users/Koustubh/Dropbox (Snow Leopard Trust)/CREEM/Nepal/Madhu_poly1.shp")
plot(boundaryMadhu)

summary(Madhu.Data2)
# Standarize Rgd on traps (this makes fits a bit more stable)
# -----------------------------------------------------------
summary(covariates(traps(Madhu.Data2)))
covariates(traps(Madhu.Data2))$stdAlt = scale(covariates(traps(Madhu.Data2))$Altitude)
summary(covariates(traps(Madhu.Data2)))
head(covariates(traps(Madhu.Data2)))

summary(Trial.Data1)







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
