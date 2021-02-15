#Spatially Explicit Capture Recapture (SECR) Rishi
library(secr)
library(fields)
library(maptools)
library(raster)
gpclibPermit()
library(rgdal)


Rishi.Data<-read.capthist(captfile="C:/Users/Koustubh/Dropbox (Snow Leopard Trust)/CREEM/Spiti/rishi_caphist.txt",
                          trapfile = "C:/Users/Koustubh/Dropbox (Snow Leopard Trust)/CREEM/Spiti/rishi_trapid_count.txt",
                          detector = "count", binary.usage = FALSE, fmt = "trapID")

RishiTraps<-traps(Rishi.Data)

#Reduced.Rishi<-reduce(Rishi, by = 'ALL', outputdetector = 'count')# Reduce function to speeden SECR fi



# Readmask with covariates included
RishiMask <- read.mask (file = "C:/Users/Koustubh/Dropbox (Snow Leopard Trust)/CREEM/Spiti/mask_clean.txt",
                        spacing = 1000, header = TRUE)
plot(RishiMask)
summary(RishiMask)
plotcovariate(RishiMask,covariate="rugged",contour=FALSE,asp=1) #Plots covariates by reading the script scrplotting.r
plotcovariate(RishiMask,covariate="Altitude", contour=FALSE, asp=1)

# Standarize Ruggedness on mask
# ------------------------------------------------------------------------
summary(covariates(RishiMask))
covariates(RishiMask)$rugged[is.na(covariates(RishiMask)$rugged)] = 0 #Make any null values on mask 0
covariates(RishiMask)$StdRgd = scale(covariates(RishiMask)$rugged) #Scales covariates to make fits a bit more stable,

# Standarize Altitude on mask
# ------------------------------------------------------------------------
summary(covariates(RishiMask))
covariates(RishiMask)$Altitude[is.na(covariates(RishiMask)$Altitude)] = 0 #Make any null values on mask 0
covariates(RishiMask)$StdAlt = scale(covariates(RishiMask)$Altitude) #Scales covariates to make fits a bit more stable,
#also makes coefficients comparable


# DENSITY MODELS with scaling of covariates
#NULL Model

windows() #Opens a new window to plot the map

# Plotting Density Surface
plot(surfaceDxy,covariate = 'D', plottype = 'shaded', scale = 2)
plot(RishiTraps, add = TRUE)
covariates(RishiTraps)

# Making raster of a density surface
RasterRishi<-raster(surfaceDxy, covariate= "D.0") #Create raster from Mask)

#Null model
Spiti.hhn<-secr.fit(Rishi.Data, detectfn="HHN", mask=RishiMask,model=list(D~1, lambda0~1, sigma~1))

#Density as function of Ruggedness
Spiti.hhn.DRgd<-secr.fit(Rishi.Data, detectfn="HHN", mask=RishiMask,model=list(D~StdRgd, lambda0~1, sigma~1))

#Density as function of Altitude
Spiti.hhn.DAlt<-secr.fit(Rishi.Data, detectfn="HHN", mask=RishiMask,model=list(D~StdAlt, lambda0~1, sigma~1))

AIC(Spiti.hhn.DRgd, Spiti.hhn.DAlt, Spiti.hhn)

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

#Slightly clunky, exponential function according to David, but still works...
# ==================

userdfn1 <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') #When function is called, more of a jargon. Tells that it is a non-euclidean function
  require(gdistance) #to load transition and geoCorrection functions
  Sraster <- raster(mask, 'noneuc') #Creates a raster from a set of coordinates and attributes and turn that into a raster. noneuc needs to be made in advance in the mask that is being used in the analysis
  ## conductance is inverse of friction
  trans <- transition(Sraster, transitionFunction = function(x) 1/mean(x),? directions = 16)
  trans <- geoCorrection(trans) #takes care of earth's curvature and also the distance differences between square and diagonally neighbouring cells
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}

summary(covariates(RishiMask))

#Constant Density and detection, non-Euclidean ranging patterns
#Using conductance from literature
Spiti.AltNonU<-secr.fit(Rishi.Data, detectfn="HHN", mask=RishiMask,
                        model=list(D~1, lambda0~1, sigma~1, noneuc ~ StdAlt -1),
                        details = list(userdist = mydistFun),
                        start = list(noneuc = 1)) #-1 gets rid of the intercept

#Using exponential cost function
Spiti.AltNonUx<-secr.fit(Rishi, detectfn="HHN", mask=RishiMask,
                         model=list(D~1, lambda0~1, sigma~1, noneuc ~ StdAlt -1),
                         details = list(userdist = userdfn1),
                         start = list(noneuc = 1)) #-1 gets rid of the intercept

#Density function of Ruggedness, Non Euclidean ranging as a function of altitude
#Using conductance from literature
Spiti.DRgd.AltNonU<-secr.fit(Rishi.Data, detectfn="HHN", mask=RishiMask,
                             model=list(D~StdRgd, lambda0~1, sigma~1, noneuc ~ StdAlt -1),
                             details = list(userdist = mydistFun),
                             start = list(noneuc = 1)) #-1 gets rid of the intercept

#Density function of Ruggedness, Non-Euc function of altitude
#Using cost function as exponential
Spiti.DRgd.AltNonUx<-secr.fit(Rishi.Data, detectfn="HHN", mask=RishiMask,
                              model=list(D~StdRgd, lambda0~1, sigma~1, noneuc ~ StdAlt -1),
                              details = list(userdist = userdfn1),
                              start = list(noneuc = 1)) #-1 gets rid of the intercept

AICSpiti<-AIC(Spiti.hhn, Spiti.AltNonU, Spiti.DRgd.AltNonU, Spiti.DRgd.AltNonUx, Spiti.AltNonUx,
              Spiti.hhn.DRgd, Spiti.hhn.DAlt)

AICSpiti

save(Spiti.hhn, Spiti.AltNonU, Spiti.DRgd.AltNonU, Spiti.DRgd.AltNonUx, Spiti.AltNonUx,
     Spiti.hhn.DRgd, Spiti.hhn.DAlt,
     file="C:/Users/Koustubh/Dropbox (Snow Leopard Trust)/CREEM/Spiti/Rishi-fits.RData") #Saves model outputs

Model2.AltNonU<-secr.fit(Rishi.Data, detectfn="HHN", mask=RishiMask,
                        model=list(D~WildPrey, lambda0~1, sigma~1, noneuc ~ StdAlt -1),
                        details = list(userdist = mydistFun),
                        start = list(noneuc = 1)) #-1 gets rid of the intercept

Model2.AltNonU <- secr.fit(Rishi, mask = RishiMask,model = list(D ~ WildPrey, g0~Z, noneuc ~ StdAlt -1),
                           details = list(userdist = mydistFun),
                           start = list(noneuc = 1))

AIC (Model2,Model2.AltNonU)
surfaceDxy <- predictDsurface(Model2.AltNonU, mask = RishiMask, se.D = FALSE, 
                              cl.D = FALSE, alpha = 0.05)
RasterRishi<-raster(surfaceDxy, covariate= "D.0") #Create raster from Mask)

plot(RasterRishi)
plot(RishiTraps, add = TRUE)