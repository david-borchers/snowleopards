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
source("D:/MEGA/PhD/r/spitisecr/scrplotting.r") #Needed for plotcovariate
setwd("C:/Users/koust/Dropbox (Snow Leopard Trust)/CREEM/Analyses/snowleopards/Spiti_Rishi")

### read capture history and trap coordinates
Rishi.data<-read.capthist('rishi_caphist.txt', 'rishi_trapid_count_covt.txt',
                     trapcovnames = c("Z"), detector = "count", binary.usage=FALSE) #read capture history with a covariate "Z" which is terrain where trap is placed

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
RishiMask2<-read.mask(file = "Rishimask_WPRgdAlt1.csv",header = TRUE)
head(covariates(RishiMask2))
summary(covariates(RishiMask2))
summary(covariates(traps(Rishi.data)))

plot(RishiMask2, covariate = "Alt", contour = FALSE, legend = FALSE)
plot(RishiMask2)

#### Run models with Rishi data
Rishi.null<-secr.fit(capthist=Rishi.data, model = list(D~1,lambda0~1, sigma~1),mask=RishiMask2, detectfn="HHN")
Rishi.null1<-secr.fit(capthist=Rishi.data1, model = list(D~1,lambda0~1, sigma~1),mask=RishiMask2, detectfn="HHN")

predict(Rishi.null)
predict(Rishi.null1)

Rishi.l0_z<-secr.fit(capthist=Rishi.data, model = list(D~1,lambda0~Z, sigma~1),mask=RishiMask2, detectfn="HHN")
AIC(Rishi.l0_z)
summary(covariates(RishiMask2))
covariates(RishiMask2)$StdAlt = scale(covariates(RishiMask2)$Alt)
covariates(RishiMask2)$StdRgd = scale(covariates(RishiMask2)$Rgd)
covariates(RishiMask2)$StdWPreyD = scale(covariates(RishiMask2)$WPreyD)
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

Rishi.WPrey.AltNonU<-secr.fit(Rishi.data, detectfn="HHN", mask=RishiMask2,
                              model=list(D~StdWPreyD, lambda0~Z, sigma~1, noneuc ~ StdAlt -1), 
                              details = list(userdist = mydistFun),
                              start = list(noneuc = 1)) #-1 gets rid of the intercept

Rishi.WPrey.RgdNonU<-secr.fit(Rishi.data, detectfn="HHN", mask=RishiMask2,
                              model=list(D~StdWPreyD, lambda0~Z, sigma~1, noneuc ~ StdRgd -1), 
                              details = list(userdist = mydistFun),
                              start = list(noneuc = 1)) #-1 gets rid of the intercept

AIC(Rishi.l0_z, Rishi.null, Rishi.D_Alt.l0_z, Rishi.D_Rgd.l0_z, Rishi.D_WPrey.l0_z, Rishi.D_Alt2.l0_z, 
    Rishi.D_WP_Alt.l0_z,Rishi.AltNonU, Rishi.D_WP_Alt_Rgd.l0_z, Rishi.WPrey.AltNonU, Rishi.WPrey.RgdNonU, 
    criterion = "AIC")
coefficients(Rishi.WPrey.AltNonU)

save(Rishi.l0_z, Rishi.null, Rishi.D_Alt.l0_z, Rishi.D_Rgd.l0_z, Rishi.D_WPrey.l0_z, Rishi.D_Alt2.l0_z, 
     Rishi.D_WP_Alt.l0_z,Rishi.AltNonU, Rishi.D_WP_Alt_Rgd.l0_z, Rishi.WPrey.AltNonU, file="Rishi-New.RData") #Saves model outputs

SL_WPD=predictDsurface(Rishi.D_WPrey.l0_z)
plot(SL_WPD)
covariates(SL_WPD)$D.0

write.csv(SL_WPD, file = "Rishi_WP.csv")
write.csv(covariates(SL_WPD), file = "Rishi_WP1.csv")

# SIGNED geometric mean (identity link)
# -------------------------------------
# need to use link = list(noneuc="identity")
# positive parameter means moving from higher to lower covariate is easier (more conductive)
# and effect is stronger at higher covariate values

signedgeommean = function(x) exp(sign(diff(x))*mean(x)) 

signedgeomLCdist <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') # required by secr
  require(gdistance) # to load transition and geoCorrection functions
  if(is.element("noneuc",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc') # Make raster from mesh with covariate 'noneuc'
  } else if(is.element("noneuc.0",names(covariates(mask)))) {
    Sraster <- raster(mask, 'noneuc.0') # Make raster from mesh with covariate 'noneuc'
  } else stop("Got to have covariate named `noneuc` or `noneuc.0` on mask.")  
  # Calculate all the conductances, using mytransfun
  trans <- transition(Sraster, transitionFunction = signedgeommean, directions = 16,symm=FALSE)
  # Adjust for distance difference between square and diagonal neighbours
  trans <- geoCorrection(trans)
  # calculate the least-cost distance (trans must be conductances)
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}

Rishi.signedGeomid<-secr.fit(capthist=Rishi.data,mask=RishiMask2, detectfn="HHN",
                             model = list(D~1,lambda0~1, sigma~1, noneuc~StdAlt-1),
                             details = list(userdist = signedgeomLCdist),
                             start=list(D=exp(-9.8),lambda0=exp(-2.1),sigma=exp(8.8),noneuc=1),
                             link = list(noneuc="identity"))

AIC(Rishi.l0_z, Rishi.null, Rishi.D_Alt.l0_z, Rishi.D_Rgd.l0_z, Rishi.D_WPrey.l0_z, Rishi.D_Alt2.l0_z, 
    Rishi.D_WP_Alt.l0_z,Rishi.D_WP_Alt_Rgd.l0_z, Rishi.WPrey.RgdNonU, Rishi.signedGeomid, 
    Rishi.AltNonU, Rishi.WPrey.AltNonU,criterion = "AIC")

save(Rishi.l0_z, Rishi.null, Rishi.D_Alt.l0_z, Rishi.D_Rgd.l0_z, Rishi.D_WPrey.l0_z, Rishi.D_Alt2.l0_z, 
     Rishi.D_WP_Alt.l0_z,Rishi.AltNonU, Rishi.D_WP_Alt_Rgd.l0_z, Rishi.WPrey.AltNonU,Rishi.WPrey.RgdNonU,
     Rishi.signedGeomid, file="Rishi-New.RData") #Saves model outputs


load("Rishi-New.RData")

coefficients(Rishi.AltNonU)
coefficients(Rishi.WPrey.AltNonU)

region.N(Rishi.AltNonU)

region.N(Rishi.D_WPrey.l0_z)

region.N(Rishi.l0_z)

region.N(Rishi.null)
plot(traps(Rishi.data), add = TRUE)


################## End of current Code August 29, 2019 ##########################
#################################################################################


#covariates(RishiMask)$Ruggedness[is.na(covariates(RishiMask)$Ruggedness)] = 0 #Make any null values on mask 0
#DENSITY MODELS ( potential g0 combinations = g0 ~ 1, g0 ~ b, g0 ~ h2, g0~Z)

#Models with variation in g0 and sigma, to pick what to use with Density models later

ModelA<-secr.fit(Rishi, mask=RishiMask, model = g0~1, sigma~1)
ModelB<-secr.fit(Rishi, mask=RishiMask, model = g0~b, sigma~1)
ModelC<-secr.fit(Rishi, mask=RishiMask, model = g0~h2, sigma~1)
ModelD<-secr.fit(Rishi, mask=RishiMask, model = g0~Z, sigma~1)
ModelE<-secr.fit(Rishi, mask=RishiMask, model = g0~t, sigma~1)


AIC(ModelA, ModelB, ModelC, ModelD, ModelE)


###Models used once 

Model1 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness, g0~t))
Model2 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ WildPrey, g0~t))
Model3 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ LargeLS, g0~t))
Model4 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ SmallLS, g0~t))
Model5 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ DisFromVill, g0~t))

Model6 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness + WildPrey, g0~t))
Model7 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness + LargeLS, g0~t))
Model8 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness + SmallLS, g0~t))
Model9 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness + DisFromVill, g0~t))

Model10 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ WildPrey + LargeLS, g0~t))
Model11 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ WildPrey + SmallLS, g0~t))
Model12 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ WildPrey + DisFromVill, g0~t))

Model13 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ LargeLS + SmallLS, g0~t))
Model14 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ LargeLS + DisFromVill, g0~t))
Model15 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ SmallLS + DisFromVill, g0~t))

Model16 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness + WildPrey + LargeLS, g0~t))
Model17 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness + WildPrey + SmallLS, g0~t))
Model18 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness + WildPrey + DisFromVill, g0~t))
Model19 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ WildPrey + LargeLS + SmallLS, g0~t))
Model20 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ WildPrey + LargeLS + DisFromVill, g0~t))
Model21 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ LargeLS + SmallLS + DisFromVill, g0~t))



Model22 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness + WildPrey + LargeLS + SmallLS, g0~t))
Model23 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness + WildPrey + LargeLS + DisFromVill, g0~t))
Model24 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ WildPrey + LargeLS + SmallLS + DisFromVill, g0~t))

Model25 <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ Ruggedness + WildPrey + LargeLS + SmallLS + DisFromVill, g0~t))



AIC (Model1, Model2, Model3, Model4, Model5, Model6,Model7, Model8, Model9, 
     Model10, Model11, Model12,Model13, Model14, Model15,  Model16, Model17, 
     Model18, Model19, Model20, Model21, Model22, Model23, Model24, Model25)

save(Model1, Model2, Model3, Model4, Model5, Model6,Model7, Model8, Model9, 
     Model10, Model11, Model12,Model13, Model14, Model15,  Model16, Model17, 
     Model18, Model19, Model20, Model21, Model22, Model23, Model24, Model25, file = DensityModels.R)

# Plotting Density Surface Straightforward method

surfaceDxy <- predictDsurface(Model2, mask = RishiMask, se.D = TRUE, 
                              cl.D = TRUE, alpha = 0.05)
plot(surfaceDxy,covariate = 'D.0', plottype = 'shaded', scale = 2)
plot(surfaceDxy,covariate = 'D', plottype = 'contour')
plot(surfaceDxy,covariate = 'D', plottype = 'persp')

plot(surfaceDxy, covariate = 'D', plottype = 'dots')
plot(RishiTraps, add = TRUE)

DRange<-predictDsurface(Model2, cl.D=TRUE, alpha=0.05) 
Z <- covariates(DRange)

X <- covariates(DRange)[which(covariates(DRange)$D.0==(max(covariates(DRange)$D.0))),]

Y <- covariates(DRange)[which(covariates(DRange)$D.0==(min(covariates(DRange)$D.0))),]
write.csv(Z, file = "D:/MEGA/PhD/r/spitisecr/densityZCIs.csv", row.names = FALSE)

write.csv(Z, file = "C:/Users/Rishi/Google Drive/PhD_/r/spitisecr/densityZCIs.csv", row.names = FALSE)
# Making raster of a density surface
RasterRishi<-raster(surfaceDxy, covariate= "D.0") #Create raster from Mask)

plot(RasterRishi)
plot(RishiTraps, add = TRUE)

writeRaster(RasterRishi,"Rishi_D.tif", format="GTiff", covariate="D.0", overwrite=TRUE)





AIC (Mode1, Mode2, Mode3, Mode4)

AIC (Model1, Model2, Model3, Model4, Model5, Model6)


####RishiSECR<-secr.fit(Rishi, mask=RishiMask)
mask.check(Model1)   #uses log-likelihood to arrive at a better estimate of the buffer
suggest.buffer(Model1)

Model2 #Use this to see results, including density estimates

SurfaceRishi<-fx.total(Model1) #fx.total to be used to plot the density surface, i.e.
# summed distribution of estimated range centres from the homogenous Poisson model (THIS IS NOT DENSITY)
###New Code Inserted###
summary(SurfaceRishi)

covariates(SurfaceRishi)$D.sum10000<-covariates(SurfaceRishi)$D.sum*10000

RasterRishi<-raster(SurfaceRishi, covariate="D.sum10000") #Create raster from Mask)

plot(RasterRishi)
plot(RishiTraps, add = TRUE)

writeRaster(RasterRishi,"Rishi_Dfinal.tif", format="GTiff", covariate="D.sum10000", overwrite=TRUE)

##########################Getting N for the area of integration###

region.N (object = Model2, region = NULL, spacing = NULL, session = NULL, group = NULL, se.N = TRUE, alpha = 0.05, loginterval = TRUE,keep.region = FALSE, nlowerbound = TRUE, RN.method = "poisson",pooled.RN = FALSE)

##Closed population estimator
closedN(Rishi, estimator = NULL, level = 0.95, maxN = 50, dmax = 10 )

??rm

esa(Model2)
derived(Model2, sessnum = NULL, groups = NULL, alpha = 0.05, 
        se.esa = FALSE, se.D = TRUE, loginterval = TRUE, 
        distribution = NULL, ncores = 1)
esa(Model2, sessnum = 1, beta = NULL, real = NULL, noccasions = NULL)
?esa
coefficients(Model2) #Lists coefficients of covariates
region.N(Model2) #Provides abundance estimates for the entire study area!

save(Model1, Model2, Model3, Model4, Model5, Model6,Model7, Model8, Model9, 
      Model10, Model11, Model12,Model13, Model14, Model15,  Model16, Model17, 
      Model18, Model19, Model20, Model21, Model22, Model23, Model24, Model25, 
     file= "/Users/rishi/Google Drive/PhD /r/spitisecr/Rishi-fits.RData") #Saves model outputs
save(Model1,Model2, file="C:/Users/Rishi/Rishi-fits.RData") #Saves model outputs