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

# Weird test, just to see if there are convergence problems when start from where left off (as other models seem to fail to converge in this case):
testfit <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                                model=list(D~stdGC, a0~Topo+Water*session, sigma~1, noneuc ~ stdGC -1), 
                                                details = list(userdist = userdfn1),
                                                start = TNNfit.DGrid.a0Waterxsess_topo.sig3)

# Check if this is any better than model with lambda0 instead of a0:
TNNfit.DGrid.lam0Waterxsess_topo.sig3 <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                                model=list(D~stdGC, lambda0~Topo+Water*session, sigma~1, noneuc ~ stdGC -1), 
                                                details = list(userdist = userdfn1), ncores=ncores,
                                                start = list(noneuc = 1))
# Try again, from where ended up:
TNNfit.DGrid.lam0Waterxsess_topo.sig3 <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                                  model=list(D~stdGC, lambda0~Topo+Water*session, sigma~1, noneuc ~ stdGC -1), 
                                                  details = list(userdist = userdfn1), ncores=ncores,
                                                  start = TNNfit.DGrid.lam0Waterxsess_topo.sig3)
# Worked fine

# Compare to equivalent a0 model:
AIC(TNNfit.DGrid.lam0Waterxsess_topo.sig3,TNNfit.DGrid.a0Waterxsess_topo.sig3)
# They are identical!


# Check if the lambda0 versions of the other models are identical:
# D~1, euc
TNNfit.DGrid.lam0Waterxsess_topo.sig0<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                              model=list(D~1, lambda0~Topo+Water*session, sigma~1), ncores=ncores)
# Conveged OK

# D~stdGC, euc
TNNfit.DGrid.lam0Waterxsess_topo.sig1<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                              model=list(D~stdGC, lambda0~Topo+Water*session, sigma~1), ncores=ncores)
#Warning messages:
#1: In autoini(ch, msk, binomN = details$binomN, ignoreusage = details$ignoreusage) :
#  'autoini' failed to find g0; setting initial g0 = 0.1
#2: In secr.fit(TNN_ch, detectfn = "HHN", list(TostMask, NoyonMask,  :
#                                                at least one variance calculation failed 


# D~1, noneuc~stdGC
TNNfit.DGrid.lam0Waterxsess_topo.sig2 <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                                model=list(D~1, lambda0~Topo+Water*session, sigma~1, noneuc ~ stdGC -1), ncores=ncores,
                                                details = list(userdist = userdfn1),
                                                start = TNNfit.DGrid.lam0Waterxsess_topo.sig0)
# Warning message:
# In secr.fit(TNN_ch, detectfn = "HHN", mask = list(TostMask, NoyonMask,  :
#                                                     possible maximization error: nlm returned code 3. See ?nlm


save(TNNfit.DGrid.a0Waterxsess_topo.sig0, TNNfit.DGrid.a0Waterxsess_topo.sig1,
     TNNfit.DGrid.a0Waterxsess_topo.sig2,TNNfit.DGrid.a0Waterxsess_topo.sig3, 
     TNNfit.DGrid.lam0Waterxsess_topo.sig0, TNNfit.DGrid.lam0Waterxsess_topo.sig1,
     TNNfit.DGrid.lam0Waterxsess_topo.sig2,TNNfit.DGrid.lam0Waterxsess_topo.sig3, 
     file = "Analysis4paper/nonEuc-mdls.RData")

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


# ====== try something similar to previous noneuc model that gave problems ==========
# First, this is the previous problematic model
TNN.GCmean_dev.WW <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                              model=list(D~rmeanGC+rmeanGCdev, lambda0~Water:Winter, sigma~1, noneuc ~ stdGC -1), 
                              details = list(userdist = userdfn1),
                              start = list(noneuc = 1))
# So now try this instead:
TNN.GCmean_dev.NE <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                              model=list(D~rmeanGC+rmeanGCdev, a0~Topo+Water*session, sigma~1, noneuc ~ stdGC -1), 
                              details = list(userdist = userdfn1),
                              start = list(noneuc = 1))

AIC(TNNfit.DGrid.a0Waterxsess_topo.sig0,TNNfit.DGrid.a0Waterxsess_topo.sig1,
    TNNfit.DGrid.a0Waterxsess_topo.sig2,TNNfit.DGrid.a0Waterxsess_topo.sig3,
    TNN.GCmean_dev.NE,TNN.GCmean_dev.WW)

# Try best model, with noneuc, separately in each region:
Tostfit.DGrid.a0Waterxsess_topo.sig <- secr.fit(Tost_ch, detectfn="HHN", mask=TostMask,
                                                model=list(D~stdGC, a0~Topo+Water, sigma~1, noneuc ~ stdGC -1), 
                                                details = list(userdist = userdfn1),
                                                start = list(noneuc = 1))
# refit, starting where left off:
Tostfit.DGrid.a0Waterxsess_topo.sig <- secr.fit(Tost_ch, detectfn="HHN", mask=TostMask,
                                                model=list(D~stdGC, a0~Topo+Water, sigma~1, noneuc ~ stdGC -1), 
                                                details = list(userdist = userdfn1),
                                                start = Tostfit.DGrid.a0Waterxsess_topo.sig)

Noyonfit.DGrid.a0Waterxsess_topo.sig <- secr.fit(Noyon_ch, detectfn="HHN", mask=NoyonMask,
                                                 model=list(D~stdGC, a0~Topo+Water, sigma~1, noneuc ~ stdGC -1), 
                                                 details = list(userdist = userdfn1),
                                                 start = list(noneuc = 1))
# refit, starting where left off:
Noyonfit.DGrid.a0Waterxsess_topo.sig <- secr.fit(Noyon_ch, detectfn="HHN", mask=NoyonMask,
                                                 model=list(D~stdGC, a0~Topo+Water, sigma~1, noneuc ~ stdGC -1), 
                                                 details = list(userdist = userdfn1),
                                                 start = Noyonfit.DGrid.a0Waterxsess_topo.sig)
# re-refit, starting where left off:
Noyonfit.DGrid.a0Waterxsess_topo.sig <- secr.fit(Noyon_ch, detectfn="HHN", mask=NoyonMask,
                                                 model=list(D~stdGC, a0~Topo+Water, sigma~1, noneuc ~ stdGC -1), 
                                                 details = list(userdist = userdfn1),
                                                 start = Noyonfit.DGrid.a0Waterxsess_topo.sig)
# Problem estimating variance; seems to arise from Steppe level of Topo factor

Nemegtfit.DGrid.a0Waterxsess_topo.sig <- secr.fit(Nemegt_ch, detectfn="HHN", mask=NemegtMask,
                                                  model=list(D~stdGC, a0~Topo+Water, sigma~1, noneuc ~ stdGC -1), 
                                                  details = list(userdist = userdfn1),
                                                  start = list(noneuc = 1))
#===================

### water:winter on a0

# D~1, euc
TNNfit.DGrid.a0WaterxWinter.sig0 <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                             model=list(D~1, a0~Water:Winter, sigma~1), 
                                             ncores = 3)

# D~stdGC, euc
TNNfit.DGrid.a0WaterxWinter.sig1 <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                             model=list(D~stdGC, a0~Water:Winter, sigma~1), 
                                             start = TNNfit.DGrid.a0WaterxWinter.sig0,
                                             ncores = 3)

# D~1, noneuc~stdGC
TNNfit.DGrid.a0WaterxWinter.sig2 <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                             model=list(D~1, a0~Water:Winter, sigma~1, noneuc ~ stdGC -1), 
                                             details = list(userdist = userdfn1),
                                             start = TNNfit.DGrid.a0WaterxWinter.sig0,
                                             ncores = 3)

# D~stdGC, noneuc~stdGC
TNNfit.DGrid.a0WaterxWinter.sig3 <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                             model=list(D~stdGC, a0~Water:Winter, sigma~1, noneuc ~ stdGC -1), 
                                             details = list(userdist = userdfn1),
                                             start = TNNfit.DGrid.a0WaterxWinter.sig0,
                                             ncores = 3)

AIC(TNNfit.DGrid.a0WaterxWinter.sig0, TNNfit.DGrid.a0WaterxWinter.sig1, 
    TNNfit.DGrid.a0WaterxWinter.sig2, TNNfit.DGrid.a0WaterxWinter.sig3)


### no model on a0 (a0~1)

# D~1, noneuc
TNNfit.check1 <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                          model=list(D~1, a0 ~ 1, sigma~1, noneuc ~ stdGC -1), 
                          details = list(userdist = userdfn1),
                          ncores = 3)

# D~1, neuc
TNNfit.check2 <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                          model=list(D~1, a0 ~ 1, sigma~1), 
                          ncores = 3)

# D~stdGC, euc
TNNfit.check3 <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                          model=list(D~stdGC, a0 ~ 1, sigma~1),
                          ncores = 3)

# D~stdGC, noneuc
TNNfit.check4 <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                          model=list(D~stdGC, a0 ~ 1, sigma~1, noneuc ~ stdGC -1), 
                          details = list(userdist = userdfn1),
                          ncores = 3)

AIC(TNNfit.check1, TNNfit.check2, TNNfit.check3, TNNfit.check4)

### separate regions

Tost.a0Waterxsess_topo.sig3 <- secr.fit(Tost_ch, detectfn="HHN", mask=TostMask,
                                        model=list(D~stdGC, a0~Topo+Water, sigma~1, noneuc ~ stdGC -1), 
                                        details = list(userdist = userdfn1))

Noyon.a0Waterxsess_topo.sig3 <- secr.fit(Noyon_ch, detectfn="HHN", mask=NoyonMask,
                                         model=list(D~stdGC, a0~Topo+Water, sigma~1, noneuc ~ stdGC -1), 
                                         details = list(userdist = userdfn1),
                                         start = Tost.a0Waterxsess_topo.sig3)

Nemegt.a0Waterxsess_topo.sig3 <- secr.fit(Nemegt_ch, detectfn="HHN", mask=NemegtMask,
                                          model=list(D~stdGC, a0~Topo+Water, sigma~1, noneuc ~ stdGC -1), 
                                          details = list(userdist = userdfn1),
                                          start = Tost.a0Waterxsess_topo.sig3)

All.a0Waterxsess_topo.sig3 <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                       model=list(D~stdGC, a0~Topo+Water*session, sigma~1, noneuc ~ stdGC*session -1), 
                                       details = list(userdist = userdfn1),
                                       start = TNNfit.DGrid.a0Waterxsess_topo.sig0)

AIC(TNNfit.DGrid.a0WaterxWinter.sig0, TNNfit.DGrid.a0WaterxWinter.sig1, 
    TNNfit.DGrid.a0WaterxWinter.sig2, TNNfit.DGrid.a0WaterxWinter.sig3,
    All.a0Waterxsess_topo.sig3)

### save all
# save(TNNfit.DGrid.a0Waterxsess_topo.sig0, TNNfit.DGrid.a0Waterxsess_topo.sig1,
#      TNNfit.DGrid.a0Waterxsess_topo.sig2,TNNfit.DGrid.a0Waterxsess_topo.sig3, 
#      TNNfit.DGrid.a0WaterxWinter.sig0, TNNfit.DGrid.a0WaterxWinter.sig1,
#      TNNfit.DGrid.a0WaterxWinter.sig2,TNNfit.DGrid.a0WaterxWinter.sig3, 
#      TNNfit.check1, TNNfit.check2, TNNfit.check3, TNNfit.check4,
#      Tost.a0Waterxsess_topo.sig3, Noyon.a0Waterxsess_topo.sig3, Nemegt.a0Waterxsess_topo.sig3,
#      All.a0Waterxsess_topo.sig3,
#      file = "Analysis4paper/nonEuc-mdls.RData")

# reduced (no a0)
zlim = range(covariates(predictDsurface(TNNfit.check4)[[1]])$`D.0`,
             covariates(predictDsurface(TNNfit.check4)[[2]])$`D.0`,
             covariates(predictDsurface(TNNfit.check4)[[3]])$`D.0`)

plot(bbxlim,bbylim,xlim=xlim,ylim=ylim,xlab="",ylab="",bty="n",type="n",xaxt="n",yaxt="n",asp=1) 
plotcovariate(predictDsurface(TNNfit.check4)[[1]],"D.0",col=parula(40), asp = 1, zlim = zlim, contour = FALSE, add = TRUE)
plotcovariate(predictDsurface(TNNfit.check4)[[2]],"D.0",col=parula(40), asp = 1, zlim = zlim, contour = FALSE, add = TRUE)
plotcovariate(predictDsurface(TNNfit.check4)[[3]],"D.0",col=parula(40), asp = 1, zlim = zlim, contour = FALSE, add = TRUE)



