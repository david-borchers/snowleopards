library(secr)
library(fields)
library(maptools)
library(gdistance)
library(raster)
library(parallel)

detectCores()

# -------------------- Start Capture Histories --------------------------
# Combined regions
getwd()
Pak.trapfiles = c("./Pakistan/Traps_Astor.csv", "./Pakistan/Traps_Basha.csv", "./Pakistan/Traps_BNPQNP.csv","./Pakistan/Traps_Chitral.csv", 
  "./Pakistan/Traps_DNP.csv", "./Pakistan/Traps_Hoper.csv", "./Pakistan/Traps_KhanBari.csv", "./Pakistan/Traps_KNP.csv",
  "./Pakistan/Traps_KNP2011.csv", "./Pakistan/Traps_MDNP.csv", "./Pakistan/Traps_Misger.csv", "./Pakistan/Traps_Misger.csv",
  "./Pakistan/Traps_Shimshal.csv", "./Pakistan/Traps_Terich.csv", "./Pakistan/Traps_Yarkhun.csv")

Pak_ch<-read.capthist(captfile = "./Pakistan/All Caps_scr.csv", 
                      binary.usage = FALSE, trapfile = Pak.trapfiles, 
                      detector="count", fmt = "trapID", 
                      trapcovnames = c("Camera","Lure", "LureType", "Habitat","Terrain"))


Pak_KNP<-read.capthist(captfile="./Pakistan/Caps_KNP2010.csv", binary.usage = FALSE,
                         trapfile = "./Pakistan/Traps_KNP.csv", detector = "count", fmt = "trapID",
                         trapcovnames = c("Camera", "Lure", "LureType", "Habitat", "Terrain"))


#Not to run
Pak_SingleSess<-read.capthist(captfile="./Pakistan/All Caps_scr_1sess.csv", binary.usage = FALSE,
                       trapfile = "./Pakistan/All Traps_scr.csv", detector = "count", fmt = "trapID",
                       trapcovnames = c("Camera", "Lure", "LureType", "Habitat", "Terrain", "Site"))

#No Chitral and Khanbari. Uses KNP data from 2 years
#Not to run
Pak_SingleSess_NoChitKhan<-read.capthist(captfile="./Pakistan/All Caps_scr_1sess.csv", binary.usage = FALSE,
                              trapfile = "./Pakistan/All Traps_scr_NoChitKhan.csv", detector = "count", fmt = "trapID",
                              trapcovnames = c("Camera", "Lure", "LureType", "Habitat", "Terrain", "Site"))

#Use this. No Chitral, Khanbari and also no data form KNP 2010
Pak_SingleSess_NoChitKhan1<-read.capthist(captfile="./Pakistan/All Caps_1sess_No2KNP.csv", binary.usage = FALSE,
                                         trapfile = "./Pakistan/All Traps_NoChitKhanKNP2010.csv", detector = "count", fmt = "trapID",
                                         trapcovnames = c("Camera", "Lure", "LureType", "Habitat", "Terrain", "Site"))


#Pak_SingleSess_NoChitKhan<-read.capthist(captfile="./Pakistan/All Caps_scr_1sess.csv", binary.usage = FALSE,
#                                         trapfile = "./Pakistan/All Traps_scr_NoChitKhan_cov1.csv", detector = "count", fmt = "trapID",
#                                         trapcovnames = c("Camera"))

summary(traps(Pak_SingleSess_NoChitKhan1))
summary(Pak_SingleSess_NoChitKhan1)
plot(Pak_SingleSess_NoChitKhan1, add = TRUE)
plot(traps(Pak_SingleSess_NoChitKhan1), add = TRUE)

# Make mask:
PakMask= make.mask(traps(Pak_SingleSess_NoChitKhan1), buffer = 35000, spacing = 500, type = "trapbuffer")
KNPMask = make.mask(traps(Pak_KNP), buffer = 30000, spacing = 500, type = "trapbuffer")

plot(PakMask)
#Write mask into CSV file to add covariates
write.csv(PakMask, file = "./Pakistan/PakMask_file.csv") #Save mask as CSV, to open and add covariates in ArcGIS. Much faster!

#Bring mask back with covariates to the exported CSV file
PakMask1 <- read.mask (file = "./Pakistan/PakMask_file_zCov.csv", spacing = 500, header = TRUE)
head(covariates(PakMask1))

plot(PakMask1)
plot(PakMask1, covariate="altitude", contour=FALSE, legend = FALSE)
plot(PakMask1, covariate="settlement", contour=FALSE, legend = FALSE)
plot(PakMask1, covariate="road", contour=FALSE, legend = FALSE)

plot(traps(Pak_SingleSess_NoChitKhan1), add = TRUE)

summary(covariates(traps(Pak_SingleSess_NoChitKhan1)))
head(covariates(traps(Pak_SingleSess_NoChitKhan1)))
head(covariates(PakMask1))

KNP_hhn<-secr.fit(Pak_KNP, model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", mask=KNPMask)
Pak_hhn<-secr.fit(Pak_SingleSess_NoChitKhan, model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", mask=PakMask)

Pak_hhn1<-secr.fit(Pak_SingleSess_NoChitKhan1, model=list(D~1, lambda0~1, sigma~1), detectfn="HHN", mask=PakMask)

Pak_hhn.D_alt<-secr.fit(Pak_SingleSess_NoChitKhan1, model=list(D~altitude, lambda0~1, sigma~1), detectfn="HHN", mask=PakMask1)
Pak_hhn.sigSite<-secr.fit(Pak_SingleSess_NoChitKhan1, model=list(D~1, lambda0~1, sigma~Site), detectfn="HHN", mask=PakMask1)

coefficients(Pak_hhn.D_alt)
AltD_Abundance<-region.N(Pak_hhn.D_alt)
AltSurface<-predictDsurface(Pak_hhn.D_alt)
plot(AltSurface, col=terrain.colors(40), contour = FALSE)
min(covariates(AltSurface))
max(covariates(AltSurface))
head(covariates(AltSurface))

?Parallel

region.N(Pak_hhn)
coefficients(Pak_hhn)
predict(Pak_hhn)

DensityConstt<-region.N(Pak_hhn1)


summary(PakMask1)

min(PakMask$x)
max(PakMask$x)
min(PakMask$y)
max(PakMask$y)

min(traps(Pak_SingleSess_NoChitKhan))
max(traps(Pak_SingleSess_NoChitKhan))
