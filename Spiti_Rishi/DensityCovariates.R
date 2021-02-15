#Spatially Explicit Capture Recapture (SECR) 13 March 2017
library(secr)
library(maptools)
library(raster)
gpclibPermit()
library(rgdal)
library(MuMIn)
library(fields)
library(rgdal)


### Working directory and scrplotting
source("scrplotting.r") #Needed for plotcovariate

# Don't know what these are, just kept here in case need to find out at some point:
#load("Rishi-fits.RData")
#load("Rishi-New (Koustubh Sharma's conflicted copy 2019-09-27).RData")

#setwd("C:/Users/koust/Dropbox (Snow Leopard Trust)/CREEM/Analyses/snowleopards/Spiti_Rishi")

### read capture history and trap coordinates for David
Rishi.data<-read.capthist('rishi_caphist.txt', 
                          'rishi_trapid_count_covt.txt',
                          trapcovnames = c("Z"), detector = "count", binary.usage=FALSE) #read capture history with a covariate "Z" which is terrain where trap is placed

### read capture history and trap coordinates for Koustubh
Rishi.data<-read.capthist('rishi_caphist.txt', 'rishi_trapid_count_covt.txt',
                          trapcovnames = c("Z"), detector = "count", binary.usage=FALSE) #read capture history with a covariate "Z" which is terrain where trap is placed
summary(Rishi.data)


plot(Rishi.data)

### Create a primary mask with buffer 25k
RishiMask1<-make.mask(traps(Rishi.data), buffer = 25000, spacing = 500, type="trapbuffer")
plot(RishiMask1)
plot(traps(Rishi.data), add = TRUE)
write.csv(RishiMask1, file = "RishiMask1csv.csv") 
TestModel1<-secr.fit(capthist=Rishi.data, model = list(D~1,lambda0~1, sigma~1),mask=RishiMask1, detectfn="HHN")
predict(TestModel1)
Suggest.buffer.Rishi<- suggest.buffer(Rishi.data, detectfn="HHN",
                                    detectpar=list(sigma=6.3625e+03, lambda0=1.7634),
                                    RBtarget = 0.001)
Suggest.buffer.Rishi

### Reading the Mask File with covariates into Density
RishiMask2<-read.mask(file = "Rishimask_WPRgdAlt1.csv",header = TRUE) #Mask with Wild Prey, Ruggedness and Altitude
RishiMask3<-read.mask(file = "Rishimask_WPRgdAltLSVill.csv",header = TRUE) #Mask with above and Livestock & village data. Prepared later

head(covariates(RishiMask2))
summary(covariates(RishiMask2))
summary(covariates(traps(Rishi.data)))

head(covariates(RishiMask3))
summary(covariates(RishiMask3))


plot(RishiMask2, covariate = "Alt", contour = FALSE, legend = FALSE)

plot(RishiMask3, covariate = "LS_Biomass", contour = FALSE, legend = FALSE)
plot(RishiMask3, covariate = "Small_LS", contour = FALSE, legend = FALSE)
plot(RishiMask3, covariate = "LargeLS", contour = FALSE, legend = FALSE)

#### Run models with Rishi data
fRishi.null1<-secr.fit(capthist=Rishi.data1, model = list(D~1,lambda0~1, sigma~1),mask=RishiMask2, detectfn="HHN")

predict(Rishi.null)
predict(Rishi.null1)

Rishi.l0_z<-secr.fit(capthist=Rishi.data, model = list(D~1,lambda0~Z, sigma~1),mask=RishiMask2, detectfn="HHN")
AIC(Rishi.l0_z)
summary(covariates(RishiMask2))
covariates(RishiMask2)$StdAlt = scale(covariates(RishiMask2)$Alt)
covariates(RishiMask2)$StdRgd = scale(covariates(RishiMask2)$Rgd)
covariates(RishiMask2)$StdWPreyD = scale(covariates(RishiMask2)$WPreyD)

names(covariates(RishiMask3))
covariates(RishiMask3)$StdAlt = scale(covariates(RishiMask3)$Alt)
covariates(RishiMask3)$StdRgd = scale(covariates(RishiMask3)$Rgd)
covariates(RishiMask3)$StdWPreyD = scale(covariates(RishiMask3)$WPreyD)
covariates(RishiMask3)$StdSmall_LS = scale(covariates(RishiMask3)$Small_LS)
covariates(RishiMask3)$StdLarge_LS = scale(covariates(RishiMask3)$LargeLS)
covariates(RishiMask3)$StdLS_Biomass = scale(covariates(RishiMask3)$LS_Biomass)
covariates(RishiMask3)$StdVillLCD = scale(covariates(RishiMask3)$Village_LC)

summary(covariates(RishiMask3))

summary(covariates(RishiMask2))
names(covariates(RishiMask2))
cor(covariates(RishiMask2)$Rgd, covariates(RishiMask2)$Alt)

Rishi.D_Alt.l0_z<-secr.fit(capthist=Rishi.data, model = list(D~StdAlt,lambda0~Z, sigma~1),mask=RishiMask2, detectfn="HHN")
Rishi.D_Rgd.l0_z<-secr.fit(capthist=Rishi.data, model = list(D~StdRgd,lambda0~Z, sigma~1),mask=RishiMask2, detectfn="HHN")
Rishi.D_WPrey.l0_z<-secr.fit(capthist=Rishi.data, model = list(D~StdWPreyD,lambda0~Z, sigma~1),mask=RishiMask2, detectfn="HHN")
Rishi.D_Alt2.l0_z<-secr.fit(capthist=Rishi.data, model = list(D~StdAlt+I(StdAlt^2),lambda0~Z, sigma~1),mask=RishiMask2, detectfn="HHN")
Rishi.D_WP_Alt.l0_z<-secr.fit(capthist=Rishi.data, model = list(D~StdAlt+StdWPreyD,lambda0~Z, sigma~1),mask=RishiMask2, detectfn="HHN")


# Non-Euclidian fits
# define functions for non-Euclidian distance calculation:
# -------------------------------------------------------

# Create a user-defined tranistion function:
myConductanceFun = function(x) exp(mean(log(x)))

# Create a user-defined function to calculate least-cost distances, from the
# mask covariate "noneuc" with the user-defined function "myConductanceFun"
mydistFun = function (xy1, xy2, mask) {
  if (missing(xy1)) return("noneuc")
  if (!require(gdistance))
    stop ("install package gdistance to use this function")
  Sraster <- raster(mask, "noneuc")
  tr <- transition(Sraster, transitionFunction=myConductanceFun, directions=16)
  tr <- geoCorrection(tr, type = "c", multpl = FALSE)
  costDistance(tr, as.matrix(xy1), as.matrix(xy2))
}

Rishi.D_WP_Alt_Rgd.l0_z<-secr.fit(capthist=Rishi.data, model = list(D~StdAlt+StdWPreyD+StdRgd,lambda0~Z, sigma~1),
                                  mask=RishiMask2, detectfn="HHN")

Rishi.AltNonU<-secr.fit(Rishi.data, detectfn="HHN", mask=RishiMask2,
                        model=list(D~1, lambda0~Z, sigma~1, noneuc ~ StdAlt -1), details = list(userdist = mydistFun),
                        start = list(noneuc = 1)) #-1 gets rid of the intercept

Rishi.WPrey.AltNonU<-secr.fit(Rishi.data, detectfn="HHN", mask=RishiMask2,
                        model=list(D~StdWPreyD, lambda0~Z, sigma~1, noneuc ~ StdAlt -1), 
                        details = list(userdist = mydistFun),
                        start = list(noneuc = 1)) #-1 gets rid of the intercept

####New Set of models with Livestock data####
names(covariates(RishiMask3))
Rishi.LargeLS.AltNonU<-secr.fit(Rishi.data, detectfn="HHN", mask=RishiMask3,
                              model=list(D~StdLarge_LS, lambda0~Z, sigma~1, noneuc ~ StdAlt -1), 
                              details = list(userdist = mydistFun),
                              start = list(noneuc = 1))

Rishi.SmallLS.AltNonU<-secr.fit(Rishi.data, detectfn="HHN", mask=RishiMask3,
                                model=list(D~StdSmall_LS, lambda0~Z, sigma~1, noneuc ~ StdAlt -1), 
                                details = list(userdist = mydistFun),
                                start = list(noneuc = 1))

Rishi.LSBio.AltNonU<-secr.fit(Rishi.data, detectfn="HHN", mask=RishiMask3,
                                model=list(D~StdLS_Biomass, lambda0~Z, sigma~1, noneuc ~ StdAlt -1), 
                                details = list(userdist = mydistFun),
                                start = list(noneuc = 1))

Rishi.WP_LSBio.AltNonU<-secr.fit(Rishi.data, detectfn="HHN", mask=RishiMask3,
                              model=list(D~StdWPreyD + StdLS_Biomass, lambda0~Z, sigma~1, noneuc ~ StdAlt -1), 
                              details = list(userdist = mydistFun),
                              start = list(noneuc = 1))

Rishi.WPxLSBio.AltNonU<-secr.fit(Rishi.data, detectfn="HHN", mask=RishiMask3,
                                 model=list(D~StdWPreyD * StdLS_Biomass, lambda0~Z, sigma~1, noneuc ~ StdAlt -1), 
                                 details = list(userdist = mydistFun),
                                 start = list(noneuc = 1))

Rishi.Vill_LC.AltNonU<-secr.fit(Rishi.data, detectfn="HHN", mask=RishiMask3,
                                 model=list(D~StdVillLCD, lambda0~Z, sigma~1, noneuc ~ StdAlt -1), 
                                 details = list(userdist = mydistFun),
                                 start = list(noneuc = 1))

Rishi.WP_Vill_LC.AltNonU<-secr.fit(Rishi.data, detectfn="HHN", mask=RishiMask3,
                                model=list(D~StdWPreyD+StdVillLCD, lambda0~Z, sigma~1, noneuc ~ StdAlt -1), 
                                details = list(userdist = mydistFun),
                                start = list(noneuc = 1))

Rishi.WPxVill_LC.AltNonU<-secr.fit(Rishi.data, detectfn="HHN", mask=RishiMask3,
                                   model=list(D~StdWPreyD*StdVillLCD, lambda0~Z, sigma~1, noneuc ~ StdAlt -1), 
                                   details = list(userdist = mydistFun),
                                   start = list(noneuc = 1))


RishiAIC<-AIC(Rishi.l0_z, Rishi.null, Rishi.D_Alt.l0_z, Rishi.D_Rgd.l0_z, Rishi.D_WPrey.l0_z, Rishi.D_Alt2.l0_z, 
    Rishi.D_WP_Alt.l0_z,Rishi.AltNonU, Rishi.D_WP_Alt_Rgd.l0_z, Rishi.WPrey.AltNonU,
    Rishi.LargeLS.AltNonU, Rishi.SmallLS.AltNonU, Rishi.LSBio.AltNonU, Rishi.WP_LSBio.AltNonU,
    Rishi.WPxLSBio.AltNonU, criterion = "AIC")
write.csv(RishiAIC, file="./RishiAIC.csv")

#Not using Village least cost models in the AIC table.
#Agreed with Rishi to not use any village data as it is not really ecologically explainable.
# Moreover least cost distance from village is somehow negatively correlated (r= -0.27) with wild prey density!
# Instead, using livestock as a surrogate for human 'disturbance' or 'impact'.

save(Rishi.l0_z, Rishi.null, Rishi.D_Alt.l0_z, Rishi.D_Rgd.l0_z, Rishi.D_WPrey.l0_z, Rishi.D_Alt2.l0_z, 
     Rishi.D_WP_Alt.l0_z,Rishi.AltNonU, Rishi.D_WP_Alt_Rgd.l0_z, Rishi.WPrey.AltNonU,
     Rishi.LargeLS.AltNonU, Rishi.SmallLS.AltNonU, Rishi.LSBio.AltNonU, Rishi.WP_LSBio.AltNonU,
     Rishi.WPxLSBio.AltNonU, Rishi.Vill_LC.AltNonU, Rishi.WP_Vill_LC.AltNonU, Rishi.WPxVill_LC.AltNonU, 
     file="Rishi-New.RData") #Saves model outputs
load("Rishi-New.RData")

coefficients(Rishi.AltNonU)
coefficients(Rishi.WPrey.AltNonU)


######################## Debugging NAs ###########################

library(tcltk2) # for progress bar
Rishi.WPrey.AltNonU.refit<-secr.fit(Rishi.data, detectfn="HHN", mask=RishiMask2,
                              model=list(D~StdWPreyD, lambda0~Z, sigma~1, noneuc ~ StdAlt -1), 
                              details = list(userdist = mydistFun),
                              start = Rishi.WPrey.AltNonU) #-1 gets rid of the intercept
coefficients(Rishi.WPrey.AltNonU.refit)
pmask = predictDsurface(Rishi.WPrey.AltNonU.refit,parameter="noneuc")
plotcovariate(pmask,"noneuc.0")
range(covariates(pmask)$noneuc.0)

ll = secr.fit(Rishi.data, detectfn="HHN", mask=RishiMask2,
              model=list(D~StdWPreyD, lambda0~Z, sigma~1, noneuc ~ StdAlt -1), 
              details = list(userdist = mydistFun,LLonly=TRUE),
              start = Rishi.WPrey.AltNonU.refit)
# Hmmm, that gives warning:
# In valid.userdist(userdist, detector(traps), xy1 = traps, xy2 = mask,  :
#                     replacing infinite, negative and NA userdist values with 1e10
covariates(RishiMask2)$noneuc = covariates(pmask)$noneuc.0 # add "noneuc" because mydistFun needs it
test = secr:::valid.userdist(mydistFun, detector(traps(Rishi.data)), xy1=traps(Rishi.data), xy2=RishiMask2, 
                             mask=RishiMask2, sessnum=1)
bads = test*0 # binary to mark bad distances
inds = which(test>=1e+10,arr.ind=TRUE)
bads[inds[,1],inds[,2]] = 1 # mark bad distances

cams = traps(Rishi.data)
starti = 1
fn = paste("traps_",starti,"_to_",(starti+3),".pdf",sep="")
pdf(file=fn)
par(mfrow=c(2,2))
for(i in starti:(starti+3)) { # choose which trap to consider
  bmask = RishiMask2 # Don't want to be adding a bunch of junk to the actual mask - add it to rm intead
  covariates(bmask)$bad = bads[i,]
  plotcovariate(bmask,"bad") # ... and plot them
  points(cams[i,],pch=19,col="red")
} 
dev.off()
# Check where inaccessible points are:
plot(bmask)
badpts = covariates(bmask)$bad==1
points(bmask$x[badpts],bmask$y[badpts],pch=19,col="red")
# Remove inaccessible points from mask
newmask = subset(RishiMask2,subset=!badpts)
plot(newmask)
# Try refitting with new mask
Rishi.WPrey.AltNonU.newmask<-secr.fit(Rishi.data, detectfn="HHN", mask=newmask,
                                      model=list(D~StdWPreyD, lambda0~Z, sigma~1, noneuc ~ StdAlt -1), 
                                      details = list(userdist = mydistFun),
                                      start = list(noneuc = 1))
# refit with last values of parameters
Rishi.WPrey.AltNonU.newmask.refit<-secr.fit(Rishi.data, detectfn="HHN", mask=newmask,
                                      model=list(D~StdWPreyD, lambda0~Z, sigma~1, noneuc ~ StdAlt -1), 
                                      details = list(userdist = mydistFun),
                                      start = coef(Rishi.WPrey.AltNonU.newmask)$beta)
coefficients(Rishi.WPrey.AltNonU.newmask.refit)
lam0 = coefficients(Rishi.WPrey.AltNonU.newmask.refit)["lambda0",1]
nu = coefficients(Rishi.WPrey.AltNonU.newmask.refit)["noneuc.StdAlt",1]
nval = 20
lam0s = seq(lam0-0.5,lam0+5,length=nval)
nus = seq(nu-5,nu+0.2,length=nval)

# Examine teh likelihood surface in lambda0 and noneuc.StdAlt dimensions to see what is going on:
system.time(tmp <- LLsurface(Rishi.WPrey.AltNonU.newmask.refit, betapar=c("lambda0", "noneuc.StdAlt"), 
                             realscale=FALSE,xval=nus, yval=lam0s, plot=TRUE))
pdf("llik_NU_lambda.pdf")
image(nus,lam0s,tmp)
contour(nus,lam0s,tmp,add=TRUE)
dev.off()

# Now refit, starting at what looks like maximum of likeliohood surface.
newstart = coef(Rishi.WPrey.AltNonU.newmask.refit)$beta # vector to hold new starting values
newstart[3] = 0.2 # new starting value for lambda0, from above plot
newstart[7] = -2.5 # new starting value for noneuc.StdAlt, from above plot
Rishi.WPrey.AltNonU.newmask.restart<-secr.fit(Rishi.data, detectfn="HHN", mask=newmask,
                                            model=list(D~StdWPreyD, lambda0~Z, sigma~1, noneuc ~ StdAlt -1), 
                                            details = list(userdist = mydistFun),
                                            start = newstart)
coefficients(Rishi.WPrey.AltNonU.newmask.restart)
# That did not work; try restarting form where left off:
Rishi.WPrey.AltNonU.newmask.restart2<-secr.fit(Rishi.data, detectfn="HHN", mask=newmask,
                                              model=list(D~StdWPreyD, lambda0~Z, sigma~1, noneuc ~ StdAlt -1), 
                                              details = list(userdist = mydistFun),
                                              start = coef(Rishi.WPrey.AltNonU.newmask.restart)$beta)
# That looks like it worked!
coefficients(Rishi.WPrey.AltNonU.newmask.restart2)

######################## Debugging NAs ###########################




coefficients(Rishi.WP_LSBio.AltNonU)
coefficients(Rishi.WPxVill_LC.AltNonU)

derived(Rishi.WPrey.AltNonU)
MapSurface<-predictDsurface(Rishi.WPrey.AltNonU)
max(covariates(MapSurface)$D.0)*10000
min(covariates(MapSurface)$D.0)*10000

AbundTop2<-region.N(Rishi.AltNonU)
AbundTop<-region.N(Rishi.WPrey.AltNonU)

write.csv(covariates(MapSurface), file = "./RishiPlot_Covs.csv") #Save top model surface as CSV, to open in ArcGIS!

hist(covariates(MapSurface)$D.0, breaks = c(0.000025, 0.00005,0.0001))

summary(RishiMask2)

region.N(Rishi.D_WPrey.l0_z)

region.N(Rishi.l0_z)

AbundNull<-region.N(Rishi.null)

AbundTop*100/5143.587

names(covariates(RishiMask3))
cor(covariates(RishiMask3)$StdAlt, covariates(RishiMask3)$StdWPreyD)
cor(covariates(RishiMask3)$StdVillLCD, covariates(RishiMask3)$StdWPreyD)
cor(covariates(RishiMask3)$StdLS_Biomass, covariates(RishiMask3)$StdWPreyD)
cor(covariates(RishiMask3)$StdLS_Biomass, covariates(RishiMask3)$StdWPreyD)

VillageLCD<-(covariates(RishiMask3)$StdVillLCD)
max(covariates(RishiMask3)$LS_Biomass)
min(covariates(RishiMask3)$LS_Biomass)
mean(covariates(RishiMask3)$LS_Biomass)
sd(covariates(RishiMask3)$LS_Biomass)

WildPrey<-covariates(RishiMask3)$StdWPreyD
Livestock<-covariates(RishiMask3)$StdLS_Biomass

cor(WildPrey, Livestock)


# ####
# #covariates(RishiMask)$Ruggedness[is.na(covariates(RishiMask)$Ruggedness)] = 0 #Make any null values on mask 0
# #DENSITY MODELS ( potential g0 combinations = g0 ~ 1, g0 ~ b, g0 ~ h2, g0~Z)
# 
# #Models with variation in g0 and sigma, to pick what to use with Density models later
# 
# ModelA<-secr.fit(Rishi, mask=RishiMask, model = g0~1, sigma~1)
# ModelB<-secr.fit(Rishi, mask=RishiMask, model = g0~b, sigma~1)
# ModelC<-secr.fit(Rishi, mask=RishiMask, model = g0~h2, sigma~1)
# ModelD<-secr.fit(Rishi, mask=RishiMask, model = g0~Z, sigma~1)
# ModelE<-secr.fit(Rishi, mask=RishiMask, model = g0~t, sigma~1)
# 
# 
# AIC(ModelA, ModelB, ModelC, ModelD, ModelE)
# 
# 
# ###Models used once 
# 
# Model1 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness, g0~t))
# Model2 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ WildPrey, g0~t))
# Model3 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ LargeLS, g0~t))
# Model4 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ SmallLS, g0~t))
# Model5 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ DisFromVill, g0~t))
# 
# Model6 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness + WildPrey, g0~t))
# Model7 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness + LargeLS, g0~t))
# Model8 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness + SmallLS, g0~t))
# Model9 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness + DisFromVill, g0~t))
# 
# Model10 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ WildPrey + LargeLS, g0~t))
# Model11 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ WildPrey + SmallLS, g0~t))
# Model12 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ WildPrey + DisFromVill, g0~t))
# 
# Model13 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ LargeLS + SmallLS, g0~t))
# Model14 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ LargeLS + DisFromVill, g0~t))
# Model15 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ SmallLS + DisFromVill, g0~t))
# 
# Model16 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness + WildPrey + LargeLS, g0~t))
# Model17 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness + WildPrey + SmallLS, g0~t))
# Model18 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness + WildPrey + DisFromVill, g0~t))
# Model19 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ WildPrey + LargeLS + SmallLS, g0~t))
# Model20 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ WildPrey + LargeLS + DisFromVill, g0~t))
# Model21 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ LargeLS + SmallLS + DisFromVill, g0~t))
# 
# 
# 
# Model22 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness + WildPrey + LargeLS + SmallLS, g0~t))
# Model23 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness + WildPrey + LargeLS + DisFromVill, g0~t))
# Model24 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ WildPrey + LargeLS + SmallLS + DisFromVill, g0~t))
# 
# Model25 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness + WildPrey + LargeLS + SmallLS + DisFromVill, g0~t))
# 
# 
# 
# AIC (Model1, Model2, Model3, Model4, Model5, Model6,Model7, Model8, Model9, 
#      Model10, Model11, Model12,Model13, Model14, Model15,  Model16, Model17, 
#      Model18, Model19, Model20, Model21, Model22, Model23, Model24, Model25)
# 
# save(Model1, Model2, Model3, Model4, Model5, Model6,Model7, Model8, Model9, 
#      Model10, Model11, Model12,Model13, Model14, Model15,  Model16, Model17, 
#      Model18, Model19, Model20, Model21, Model22, Model23, Model24, Model25, file = DensityModels.R)
# 
# # Plotting Density Surface Straightforward method
# 
# surfaceDxy <- predictDsurface(Model2, mask = RishiMask, se.D = TRUE, 
#                               cl.D = TRUE, alpha = 0.05)
# plot(surfaceDxy,covariate = 'D.0', plottype = 'shaded', scale = 2)
# plot(surfaceDxy,covariate = 'D', plottype = 'contour')
# plot(surfaceDxy,covariate = 'D', plottype = 'persp')
# 
# plot(surfaceDxy, covariate = 'D', plottype = 'dots')
# plot(RishiTraps, add = TRUE)
# 
# DRange<-predictDsurface(Model2, cl.D=TRUE, alpha=0.05) 
# Z <- covariates(DRange)
# 
# X <- covariates(DRange)[which(covariates(DRange)$D.0==(max(covariates(DRange)$D.0))),]
# 
# Y <- covariates(DRange)[which(covariates(DRange)$D.0==(min(covariates(DRange)$D.0))),]
# write.csv(Z, file = "D:/MEGA/PhD/r/spitisecr/densityZCIs.csv", row.names = FALSE)
# 
# write.csv(Z, file = "C:/Users/Rishi/Google Drive/PhD_/r/spitisecr/densityZCIs.csv", row.names = FALSE)
# # Making raster of a density surface
# RasterRishi<-raster(surfaceDxy, covariate= "D.0") #Create raster from Mask)
# 
# plot(RasterRishi)
# plot(RishiTraps, add = TRUE)
# 
# writeRaster(RasterRishi,"Rishi_D.tif", format="GTiff", covariate="D.0", overwrite=TRUE)
# 
# 
# 
# 
# 
# AIC (Mode1, Mode2, Mode3, Mode4)
# 
# AIC (Model1, Model2, Model3, Model4, Model5, Model6)
# 
# 
# ####RishiSECR<-secr.fit(Rishi, mask=RishiMask)
# mask.check(Model1)   #uses log-likelihood to arrive at a better estimate of the buffer
# suggest.buffer(Model1)
# 
# Model2 #Use this to see results, including density estimates
# 
# SurfaceRishi<-fx.total(Model1) #fx.total to be used to plot the density surface, i.e.
# # summed distribution of estimated range centres from the homogenous Poisson model (THIS IS NOT DENSITY)
# ###New Code Inserted###
# summary(SurfaceRishi)
# 
# covariates(SurfaceRishi)$D.sum10000<-covariates(SurfaceRishi)$D.sum*10000
# 
# RasterRishi<-raster(SurfaceRishi, covariate="D.sum10000") #Create raster from Mask)
# 
# plot(RasterRishi)
# plot(RishiTraps, add = TRUE)
# 
# writeRaster(RasterRishi,"Rishi_Dfinal.tif", format="GTiff", covariate="D.sum10000", overwrite=TRUE)
# 
# ##########################Getting N for the area of integration###
# 
# region.N (object = Model2, region = NULL, spacing = NULL, session = NULL, group = NULL, se.N = TRUE, alpha = 0.05, loginterval = TRUE,keep.region = FALSE, nlowerbound = TRUE, RN.method = "poisson",pooled.RN = FALSE)
# 
# ##Closed population estimator
# closedN(Rishi, estimator = NULL, level = 0.95, maxN = 50, dmax = 10 )
# 
# ??rm
# 
# esa(Model2)
# derived(Model2, sessnum = NULL, groups = NULL, alpha = 0.05, 
#         se.esa = FALSE, se.D = TRUE, loginterval = TRUE, 
#         distribution = NULL, ncores = 1)
# esa(Model2, sessnum = 1, beta = NULL, real = NULL, noccasions = NULL)
# ?esa
# coefficients(Model2) #Lists coefficients of covariates
# region.N(Model2) #Provides abundance estimates for the entire study area!
# 
# save(Model1, Model2, Model3, Model4, Model5, Model6,Model7, Model8, Model9, 
#       Model10, Model11, Model12,Model13, Model14, Model15,  Model16, Model17, 
#       Model18, Model19, Model20, Model21, Model22, Model23, Model24, Model25, 
#      file= "/Users/rishi/Google Drive/PhD /r/spitisecr/Rishi-fits.RData") #Saves model outputs
# save(Model1,Model2, file="C:/Users/Rishi/Rishi-fits.RData") #Saves model outputs