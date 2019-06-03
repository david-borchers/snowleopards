library(secr)
library(fields)
library(maptools)
library(pals)
source("scrplotting.r")

#----------------------- Load the RData objects made by Make_TNN_RData.r --------------------------
load("./Analysis4paper/TNN_boundaries.RData") # Tostboundary,Noyonboundary,Nemegtboundary
load("./Analysis4paper/TNN_masks.RData") # TostMask,NoyonMask,NemegtMask
load("./Analysis4paper/TNN_caphists.RData") # Tost_ch,Noyon_ch,Nemegt_chTNN_ch
#----------------------- ----------------------------------------------- --------------------------

#----------------------- Do some plots to check data seems OK --------------------------

# Find extent of bounding boxes of boundarys:
bbox.Nemegt = bbox(Nemegtboundary)
bbox.Noyon = bbox(Noyonboundary)
bbox.Tost = bbox(Tostboundary)
bbxlim = range(bbox.Nemegt["x",],bbox.Noyon["x",],bbox.Tost["x",]) #Set plot limit for all 3 areas together
bbylim = range(bbox.Nemegt["y",],bbox.Noyon["y",],bbox.Tost["y",])

dxlim = diff(bbxlim)
xlim = c(bbxlim[1],bbxlim[2]+0.15*dxlim)
ylim = bbylim

# ----------------------- Individual region models ---------------------------
# This taken straight from secr vignette:
userdfn1 <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') #When function is called, more of a jargon. Tells that it is a non-euclidean function
  require(gdistance) #to load transition and geoCorrection functions
  Sraster <- raster(mask, 'noneuc') #Creates a raster from a set of coordinates and attributes and turn that into a raster. noneuc needs to be made in advance in the mask that is being used in the analysis
  ## conductance is inverse of friction
  trans <- transition(Sraster, transitionFunction = function(x) 1/mean(x),  directions = 16)
  trans <- geoCorrection(trans) #takes care of earth's curvature and also the distance differences between square and diagonally neighbouring cells
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}

# D~1, euc
TNNfit.DGrid.a0Waterxsess_topo.sig0<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                              model=list(D~1, a0~Topo+Water*session, sigma~1), ncores=2)

# D~stdGC, euc
TNNfit.DGrid.a0Waterxsess_topo.sig1<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                              model=list(D~stdGC, a0~Topo+Water*session, sigma~1), ncores=2)

# D~1, noneuc~stdGC
TNNfit.DGrid.a0Waterxsess_topo.sig2 <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                              model=list(D~1, a0~Topo+Water*session, sigma~1, noneuc ~ stdGC -1), 
                              details = list(userdist = userdfn1),
                              start = TNNfit.DGrid.a0Waterxsess_topo.sig0)

# D~stdGC, noneuc~stdGC
TNNfit.DGrid.a0Waterxsess_topo.sig3 <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                                model=list(D~stdGC, a0~Topo+Water*session, sigma~1, noneuc ~ stdGC -1), 
                                                details = list(userdist = userdfn1),
                                                start = TNNfit.DGrid.a0Waterxsess_topo.sig0)

#save(TNNfit.DGrid.a0Waterxsess_topo.sig0, TNNfit.DGrid.a0Waterxsess_topo.sig1,
#     TNNfit.DGrid.a0Waterxsess_topo.sig2,TNNfit.DGrid.a0Waterxsess_topo.sig3, 
#     file = "Analysis4paper/nonEuc-mdls.RData")

load("Analysis4paper/nonEuc-mdls.RData")

AIC(TNNfit.DGrid.a0Waterxsess_topo.sig0,TNNfit.DGrid.a0Waterxsess_topo.sig1,
    TNNfit.DGrid.a0Waterxsess_topo.sig2,TNNfit.DGrid.a0Waterxsess_topo.sig3)

summary(TNNfit.DGrid.a0Waterxsess_topo.sig3)

zlim = range(covariates(predictDsurface(TNNfit.DGrid.a0Waterxsess_topo.sig3)[[1]])$`D.0`,
             covariates(predictDsurface(TNNfit.DGrid.a0Waterxsess_topo.sig3)[[2]])$`D.0`,
             covariates(predictDsurface(TNNfit.DGrid.a0Waterxsess_topo.sig3)[[3]])$`D.0`)

plot(bbxlim,bbylim,xlim=xlim,ylim=ylim,xlab="",ylab="",bty="n",type="n",xaxt="n",yaxt="n",asp=1) 
plotcovariate(predictDsurface(TNNfit.DGrid.a0Waterxsess_topo.sig3)[[1]],"D.0",col=parula(40), asp = 1, zlim = zlim, contour = FALSE, add = TRUE)
plotcovariate(predictDsurface(TNNfit.DGrid.a0Waterxsess_topo.sig3)[[2]],"D.0",col=parula(40), asp = 1, zlim = zlim, contour = FALSE, add = TRUE)
plotcovariate(predictDsurface(TNNfit.DGrid.a0Waterxsess_topo.sig3)[[3]],"D.0",col=parula(40), asp = 1, zlim = zlim, contour = FALSE, add = TRUE)

# not sure how to show noneuc results when D~1
plotcovariate(predictDsurface(TNNfit.DGrid.a0Waterxsess_topo.sig2)[[1]],"D.0",col=parula(40))
