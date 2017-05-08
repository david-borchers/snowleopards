library(secr)
library(fields)
library(maptools)
source("scrplotting.r")
library(gdistance)

NepalCaps<-read.csv("./Nepal Genetic/Madhu_Caps1.txt",header=TRUE)
NepalTraps<-read.traps("./Nepal Genetic/Madhu_T1.txt", detector = "transect")

make.capthist(captures=NepalCaps,traps=NepalTraps,fmt="XY")
NepalTraps2<-read.traps("./Nepal Genetic/Madhu_T3.txt", detector = "transect")


Trial.Data1<-read.capthist(captfile = "./Nepal Genetic/Madhu_Caps2.txt",
                          trapfile ="./Nepal Genetic/Madhu_T2.txt", detector="transect",fmt = "XY")

summary(Trial.Data)

Trial.Data2<-read.capthist(captfile = "./Nepal Genetic/Madhu_Caps2.txt",binary.usage=FALSE,
                          trapfile ="./Nepal Genetic/Madhu_T4.txt", detector="transect",fmt = "XY",
                          trapcovnames = c("Topography","Habitats","Altitude"))

summary(Trial.Data2)

Trial.Data1<-read.capthist(captfile = "./Nepal Genetic/trialCapscsv1.csv",
                          trapfile ="./Nepal Genetic/TrialTrapscsv1.csv", detector="transect",fmt = "XY")

Trial.Data2<-read.capthist(captfile = "./Nepal Genetic/trialCapscsv1.csv",
                          trapfile ="./Nepal Genetic/TrialTrapscsv3.csv", detector="transect",fmt = "XY")

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
