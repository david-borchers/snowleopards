# save(all.data.Tost,all.data.Tost,TostMask1,userdfn1,Tost.hhn.DHab.nonU.GBx,file="TostBestFit")
load("TostBestFit") # load all objects needed

library(secr)
library(fields)
library(maptools)
source("scrplotting.r")
library(gdistance)

# Here's userdfn1 (should not need - should have loaded above)
#userdfn1 <- function (xy1, xy2, mask) {
#  if (missing(xy1)) return('noneuc') #When function is called, more of a jargon. Tells that it is a non-euclidean function
#  require(gdistance) #to load transition and geoCorrection functions
#  Sraster <- raster(mask, 'noneuc') #Creates a raster from a set of coordinates and attributes and turn that into a raster. noneuc needs to be made in advance in the mask that is being used in the analysis
#  ## conductance is inverse of friction
#  trans <- transition(Sraster, transitionFunction = function(x) 1/mean(x),  directions = 16)
#  trans <- geoCorrection(trans) #takes care of earth's curvature and also the distance differences between square and diagonally neighbouring cells
#  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
#}

# Here's the code that created Tost.hhn.DHab.nonU.GBx
#Tost.hhn.DHab.nonU.GBx<-secr.fit(all.data.Tost, detectfn="HHN", mask=TostMask1,
#                                 model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdBC-1), 
#                                 details = list(userdist = userdfn1),
#                                 start = list(noneuc = 1),trace=0)

FXTost<-fx.total(Tost.hhn.DHab.nonU.GBx)
# Plot using plot.Dsurface from secr 
# (Note scrplotting.r also has a function called plot.Dsurface, but it has diffent arguments.)
secr:::plot.Dsurface(FXTost, covariate = 'D.sum', breaks = seq(0,10e-5,1e-5), poly = FALSE)
plot(Tost.hhn.DHab.nonU.GBx$capthist,add=TRUE,tracks=TRUE)
plot(traps(Tost.hhn.DHab.nonU.GBx$capthist),add=TRUE)

# Plot using plotcovariate from secrplotting.r
plotcovariate(FXTost, covariate = 'D.sum',xaxt="n",yaxt="n",xlab="",ylab="",
              bty="n",contour=FALSE,asp=1)
plotMaskEdge(Tost.hhn.DHab.nonU.GBx$mask,add=TRUE)
plot(Tost.hhn.DHab.nonU.GBx$capthist,add=TRUE,tracks=TRUE)
plot(traps(Tost.hhn.DHab.nonU.GBx$capthist),add=TRUE)
# Also consider log(D.sum)
covariates(FXTost)$logD.sum = log(covariates(FXTost)$D.sum)
plotcovariate(FXTost, covariate = 'logD.sum',xaxt="n",yaxt="n",xlab="",ylab="",
              bty="n",contour=FALSE,asp=1)
plotMaskEdge(Tost.hhn.DHab.nonU.GBx$mask,add=TRUE)
plot(Tost.hhn.DHab.nonU.GBx$capthist,add=TRUE,tracks=TRUE)
plot(traps(Tost.hhn.DHab.nonU.GBx$capthist),add=TRUE)

# Compare to density or log(density) plot:
TostDhat = predictDsurface(Tost.hhn.DHab.nonU.GBx)
covariates(TostDhat)$logD.0 = log(covariates(TostDhat)$D.0)
plotcovariate(TostDhat, covariate = 'D.0',xaxt="n",yaxt="n",xlab="",ylab="",
              bty="n",contour=FALSE,asp=1)
plotcovariate(TostDhat, covariate = 'logD.0',xaxt="n",yaxt="n",xlab="",ylab="",
              bty="n",contour=FALSE,asp=1)
plotMaskEdge(Tost.hhn.DHab.nonU.GBx$mask,add=TRUE)
plot(Tost.hhn.DHab.nonU.GBx$capthist,add=TRUE,tracks=TRUE)
plot(traps(Tost.hhn.DHab.nonU.GBx$capthist),add=TRUE)


