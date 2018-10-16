library(secr)
library(pals)
library(raster)
library(spatstat)
library(maptools)
library(gdistance) 


if(.Platform$OS.type=="windows") { # this to make the command quartz( ) work on windows machines
  quartz<-function() windows()
}

# Get some user functions for using noneuc in secr:
# ================================================
source("./Analysis4paper/noneuc-utils.R")

# Get the Tost, Noyon and Nemegt data
# ===================================
##----------------------- Load the RData objects made by Make_TNN_RData.r --------------------------
#load("./Analysis4paper/TNN_boundaries.RData") # Tostboundary,Noyonboundary,Nemegtboundary
#load("./Analysis4paper/TNN_masks.RData") # TostMask,NoyonMask,NemegtMask
#load("./Analysis4paper/TNN_caphists.RData") # Tost_ch,Noyon_ch,Nemegt_sh,TNN_ch)
##----------------------- ----------------------------------------------- --------------------------
## Bundle them together in a single Rds file for convenience:
#Tost = list(capthist=Tost_ch,mask=TostMask,boundary=Tostboundary)
#Noyon = list(capthist=Noyon_ch,mask=NoyonMask,boundary=Noyonboundary)
#Nemegt = list(capthist=Nemegt_ch,mask=NemegtMask,boundary=Nemegtboundary)
#saveRDS(list(Tost=Tost,Noyon=Noyon,Nemegt=Nemegt),file="./Analysis4paper/TNN.Rds")
dat = readRDS("./Analysis4paper/TNN.Rds")
Tost = dat$Tost
Nemegt = dat$Nemegt
Noyon = dat$Noyon

# Some plotting of data:
# Find extent of bounding boxes of boundarys:
bbox.Nemegt = bbox(Nemegt$boundary)
bbox.Noyon = bbox(Noyon$boundary)
bbox.Tost = bbox(Tost$boundary)
bbxlim = range(bbox.Noyon["x",],bbox.Tost["x",],bbox.Nemegt["x",]) #Set plot limit for all 3 areas together
bbylim = range(bbox.Noyon["y",],bbox.Tost["y",],bbox.Nemegt["y",])
dxlim = diff(bbxlim)
xlim = c(bbxlim[1],bbxlim[2]+0.15*dxlim)
ylim = bbylim
# Plot the 3 regions 
quartz(w=10,h=6)
plot(xlim,ylim,xlim=xlim,ylim=ylim,type="n",asp=1,bty="n",xlab="Easting",ylab="Northing") 
plotcovariate(Tost$mask,covariate="stdGC",what="image",add=TRUE,col=terrain.colors(30))
plotcovariate(Noyon$mask,covariate="stdGC",what="image",add=TRUE,col=terrain.colors(30))
plotcovariate(Nemegt$mask,covariate="stdGC",what="image",add=TRUE,col=terrain.colors(30))
plot(traps(Tost$capthist),add=TRUE)
plot(traps(Noyon$capthist),add=TRUE)
plot(traps(Nemegt$capthist),add=TRUE)
plot(Tost$capthist, tracks=TRUE, add=TRUE,cappar=list(cex=0.5),varycol=FALSE)
plot(Noyon$capthist, tracks=TRUE, add=TRUE,cappar=list(cex=0.5),varycol=FALSE)
plot(Nemegt$capthist, tracks=TRUE, add=TRUE,cappar=list(cex=0.5),varycol=FALSE)



# Consider Tost data
# ==================
ch = dat$Tost$capthist
mask = dat$Tost$mask
boundary = dat$Tost$boundary
cams = traps(ch)

## Fit a LC distance model with flat D:
## ===================================
## (startvals from previous fit)
#startvals = list(D=exp(-9.4515904),lambda0=exp(-4.2505931),sigma=exp(8.6914951),noneuc=0.3314881)
#sl.ne <-secr.fit(ch, detectfn="HHN", mask=mask,
#                 model=list(D~1, lambda0~1, sigma~1, noneuc~stdGC-1), 
#                 details = list(userdist = geomLCdist),
#                 start=startvals,
#                 link = list(noneuc="identity"))
## Also fit a Euclidian distance model with varying D
## ==================================================
## (startvals from previous fit)
#startvals = list(D=exp(-9.5227750),lambda0=exp(-4.4012107),sigma=exp(8.8538228))
#sl.nuD <-secr.fit(ch, detectfn="HHN", mask=mask, start=startvals,
#                  model=list(D~stdGC, lambda0~1, sigma~1))
#nuD = predictDsurface(sl.nuD)
#plotcovariate(nuD,"D.0",col=parula(40))
#region.N(sl.nuD)
##estimate SE.estimate       lcl      ucl  n
##E.N 15.48028    4.158663  9.227452 25.97025 14
##R.N 15.48035    1.346921 14.323074 20.78309 14
## Save these two density surfaces
#saveRDS(sl.nuD,file="sl.nuD.Rds")
#saveRDS(sl.ne,file="sl.ne.Rds")

# Load these two density surfaces
sl.nuD = readRDS("sl.nuD.Rds")
sl.ne = readRDS("sl.ne.Rds")

# get the predictd density surface
nuD = predictDsurface(sl.nuD)

# Simulate a population using above density surface from sl.nuD
# =============================================================
E.N=45 # Set the expected number of snow leopards in the region
scale = E.N/region.N(sl.nuD)["R.N","estimate"]
# Plot the density
simD = covariates(nuD)$"D.0" * scale
# Simulate a population:
seed=31 # for repeatability of this one population
pop = sim.popn(simD,core=mask,poly=boundary,model2D="IHP",Ndbn="fixed",seed=seed)
plot(boundary)
plot(pop,add=TRUE,pch="+")
dim(pop)[1]

# Set noneuc and detection prarameters to give stronger/weaker dependince on stdGC.
# (Parameter sl.fudge$fit$par[4] is the \alpha_2 parameter of Sutherland et al (2015).)
# (Fiddle with sl.fudge$fit$par[2] (lambda_0) and sl.fudge$fit$par[3] (sigma) 
# to give right sort of sample size:)
# ============================================================================================
sl.ne$fit$par
sl.fudge = sl.ne
sl.fudge$fit$par[4] = 2
sl.fudge$fit$par[3] = sl.ne$fit$par[3] - log(3)
sl.fudge$fit$par[2] = sl.ne$fit$par[2] - log(15)
# Plot to see if you have the kind of habitat use that you're happy with:
ne.fudge = predictDsurface(sl.fudge,parameter="noneuc")
par(mfrow=c(1,1))
# Note: lcusageplot requires you to click on a point to create an activity centre; click again to reset plot
# Do this repeatedly; hit esc when you have had enough (you get an error message - have not figured out how to exit more elegantly yet)
lcusageplot(sl.ne,mask=ne.fudge,col=parula(15))
# Note: lcusageplot requires you to click on two points to create start and end points
# Do this repeatedly; hit esc when you have had enough (you get an error message - have not figured out how to exit more elegantly yet)
lcpathplot(mask=ne.fudge,"geommean",col=parula(20),main="noneuc",what="image")


# Simulate capture histories
# ==========================
udist = geomLCdist(cams,pop,ne.fudge) # first set least-cost distances
covariates(ne.fudge)$noneuc = covariates(ne.fudge)$noneuc.0 # create the noneuc variable needed for LC dist
seed = 1 # for reproducibility
ch = sim.capthist(cams,pop,detectfn="HHN",detectpar=detectpar(sl.fudge),noccasions=1,userdist=udist,seed=seed)
plotcovariate(ne.fudge,"noneuc",col=parula(40))
plot(boundary,add=TRUE)
plot(ch,tracks=TRUE,add=TRUE,rad=10000000,type="petal")
plot(cams,add=TRUE)
summary(ch)

# Now fit some models and look at AICs and N estimates:

## Fit constant density model with LC distance:
#sl.ne.fit <-secr.fit(ch, detectfn="HHN", mask=mask,
#                     model=list(D~1, lambda0~1, sigma~1, noneuc~stdGC-1), 
#                     details = list(userdist = geomLCdist),
#                     start=list(lambda0=exp(sl.fudge$fit$par[2]),
#                                sigma=exp(sl.fudge$fit$par[3]),
#                                noneuc=sl.fudge$fit$par[4]),
#                     link = list(noneuc="identity"))
#sl.noneuc = predictDsurface(sl.ne.fit,parameter="noneuc")
## Look at some usage and least-cost paths:
#par(mfrow=c(1,1))
#lcusageplot(sl.ne.fit,mask=sl.noneuc,col=parula(15))
#lcpathplot(mask=sl.noneuc,"geommean",col=parula(20),main="noneuc",what="image")

## Fit varying density model with LC distance:
#sl.ne.fit.D <-secr.fit(ch, detectfn="HHN", mask=mask,
#                       model=list(D~stdGC, lambda0~1, sigma~1, noneuc~stdGC-1), 
#                       details = list(userdist = geomLCdist),
#                       start=list(lambda0=exp(sl.fudge$fit$par[2]),
#                                  sigma=exp(sl.fudge$fit$par[3]),
#                                  noneuc=sl.fudge$fit$par[4]),
#                       link = list(noneuc="identity"))
#sl.noneuc.D = predictDsurface(sl.ne.fit.D)
#plotcovariate(sl.noneuc.D,"D.0",col=parula(40))
#sl.noneuc.D = predictDsurface(sl.ne.fit,parameter="noneuc")
## Look at some usage and least-cost paths:
#par(mfrow=c(1,1))
#lcpathplot(mask=sl.noneuc.D,"geommean",col=parula(20),main="noneuc",what="image")
#lcusageplot(sl.ne.fit.D,mask=sl.noneuc.D,col=parula(15))

## Fit a null model:
#sl.ne.0 <-secr.fit(ch, detectfn="HHN", mask=mask,
#                   model=list(D~1, lambda0~1, sigma~1))

## Fit varying density model with Euclidian distance:
#sl.fit.D <-secr.fit(dat$Tost$capthist, detectfn="HHN", mask=dat$Tost$mask,
#                    model=list(D~stdGC, lambda0~1, sigma~1))
#sl.D = predictDsurface(sl.fit.D)
#plotcovariate(sl.D,"D.0",col=parula(40))

## Save these fits
#saveRDS(sl.ne.fit,file="sl.ne.fit.Rds")
#saveRDS(sl.ne.fit.D,file="sl.ne.fit.D.Rds")
#saveRDS(sl.ne.0,file="sl.ne.0.Rds")
#saveRDS(sl.fit.D,file="sl.fit.D.Rds")
# Read these fits
sl.ne.fit = readRDS("sl.ne.fit.Rds")
sl.ne.fit.D = readRDS("sl.ne.fit.D.Rds")
sl.ne.0 = readRDS("sl.ne.0.Rds")
sl.fit.D = readRDS("sl.fit.D.Rds")

# Compare AICs:
AIC(sl.ne.0,sl.ne.fit,sl.ne.fit.D,sl.fit.D)
#> AIC(sl.ne.0,sl.ne.fit,sl.ne.fit.D,sl.fit.D)
#model          detectfn npar    logLik     AIC    AICc   dAICc AICcwt
#sl.ne.fit       D~1 lambda0~1 sigma~1 noneuc~stdGC - 1 hazard halfnormal    4 -145.1693 298.339 300.692   0.000 0.8424
#sl.ne.fit.D D~stdGC lambda0~1 sigma~1 noneuc~stdGC - 1 hazard halfnormal    5 -145.1472 300.294 304.044   3.352 0.1576
#sl.ne.0                          D~1 lambda0~1 sigma~1 hazard halfnormal    3 -152.7691 311.538 312.871  12.179 0.0000
#sl.fit.D                     D~stdGC lambda0~1 sigma~1 hazard halfnormal    4 -210.1210 428.242 432.686 131.994 0.0000

region.N(sl.ne.fit)
region.N(sl.ne.0)
region.N(sl.ne.fit.D)
region.N(sl.fit.D)


### ===========================================================================================================
### Ian: The stuff below here I did quite a while ago - can't really remember what it was or how finished it is
###      Might be worth looking at - I am not sure
### ===========================================================================================================


# Do a bunch of simulations with flat density surface and Euclidian distance:
require(tcltk)

# TOST Non-uniform density, alpha=4, N=45
E.N=45 # Set the expected number of snow leopards in the region
scale = E.N/region.N(sl.nuD)["R.N","estimate"]
# Plot the density
simD = covariates(nuD)$"D.0" * scale
nsim=100
Nest = cover = rep(NA,nsim)
# create progress bar
pb <- tkProgressBar(title="progress bar", min=0, max=nsim, width=300)
seed=31 # for repeatability of this one population
for(i in 1:nsim) {
  if(i==1) pop = sim.popn(simD,core=mask,poly=boundary,model2D="IHP",Ndbn="fixed",seed=seed)
  else   pop = sim.popn(simD,core=mask,poly=boundary,model2D="IHP",Ndbn="fixed")
  udist = geomLCdist(cams,pop,ne.fudge) # first set least-cost distances
  covariates(ne.fudge)$noneuc = covariates(ne.fudge)$noneuc.0 # create the noneuc variable needed for LC dist
  ch = sim.capthist(cams,pop,detectfn="HHN",detectpar=detectpar(sl.fudge),noccasions=1,userdist=udist)
  simfit.0 <-secr.fit(ch, detectfn="HHN", mask=mask,  model=list(D~1, lambda0~1, sigma~1),trace=0)
  Nest.0 = region.N(simfit.0)
  Nest[i] = Nest.0["R.N","estimate"]
  cover[i] = (Nest.0["R.N","lcl"]<=E.N & E.N<=Nest.0["R.N","ucl"])
  setTkProgressBar(pb, i, title=paste( round(i/nsim*100, 0),"% done"))
}
close(pb)
hist(Nest,nclass=30)
mean(Nest)
sum(cover)/nsim
sim100_N45_alpha4 = list(Nest=Nest,cover=cover,trueD=ne.fudge,seed=seed)
saveRDS(sim100_N45_alpha4,file="./Analysis4paper/sims/Tost_sim100_N45_alpha4.Rds")


# TOST Non-uniform density, alpha=0, N=16
# Fit varying density model with Euclidian distance:
sl.fit.D <-secr.fit(dat$Tost$capthist, detectfn="HHN", mask=dat$Tost$mask,
                    model=list(D~stdGC, lambda0~1, sigma~1))
detpars = detectpar(sl.fit.D)
sl.D = predictDsurface(sl.fit.D)
meanD = mean(covariates(sl.D)$D.0)
plotcovariate(sl.D,"D.0",col=parula(40))
E.N=ceiling(region.N(sl.nuD)["R.N","estimate"]) 
# Set the expected number of snow leopards in the region
scale = E.N/region.N(sl.nuD)["R.N","estimate"]
# Plot the density
simD = covariates(sl.D)$"D.0" * scale
meanD = mean(covariates(sl.D)$"D.0")
nsim=100
Nest = cover = rep(NA,nsim)
# create progress bar
pb <- tkProgressBar(title="progress bar", min=0, max=nsim, width=300)
seed=31 # for repeatability of this one population
for(i in 1:nsim) {
  if(i==1) pop = sim.popn(simD,core=mask,poly=boundary,model2D="IHP",Ndbn="fixed",seed=seed)
  else   pop = sim.popn(simD,core=mask,poly=boundary,model2D="IHP",Ndbn="fixed")
#  udist = geomLCdist(cams,pop,ne.fudge) # first set least-cost distances
#  covariates(ne.fudge)$noneuc = covariates(ne.fudge)$noneuc.0 # create the noneuc variable needed for LC dist
  ch = sim.capthist(cams,pop,detectfn="HHN",detectpar=detpars,noccasions=1)
  simfit.0 <-secr.fit(ch, detectfn="HHN", mask=mask,  model=list(D~1, lambda0~1, sigma~1), 
                      start=list(D=meanD, lambda=detpars$lambda0, sigma~detpars$sigma), trace=0)
  Nest.0 = region.N(simfit.0)
  Nest[i] = Nest.0["R.N","estimate"]
  cover[i] = (Nest.0["R.N","lcl"]<=E.N & E.N<=Nest.0["R.N","ucl"])
  setTkProgressBar(pb, i, title=paste( round(i/nsim*100, 0),"% done"))
}
close(pb)
hist(Nest,nclass=30)
100*(mean(Nest)-E.N)/E.N
sum(cover)/nsim
sim100_N16_alpha0 = list(Nest=Nest,cover=cover,trueD=simD,seed=seed)
saveRDS(sim100_N16_alpha0,file="./Analysis4paper/sims/Tost_sim100_N16_alpha0.Rds")



# ALL 3 REGIONS  -- still only partially coded !!!
TNN.ch<-read.capthist(captfile = "./Tost_Noyon_Nemegt/TNN_Capture.csv", 
                            binary.usage = FALSE, trapfile = TNN.trapfiles, 
                            detector="count", fmt = "trapID", 
                            trapcovnames = c("Rgd","Topo", "Water", "Winter"))
mask=list(dat$Tost$mask, dat$Noyon$mask, dat$Nemegt$mask)
TNN.GCmean <- secr.fit(TNN.ch, detectfn="HHN", mask=mask, model=list(D~rmeanGC, lambda0~Water, sigma~1))

detpars = detectpar(TNN.GCmean)
TNN.GC.D = predictDsurface(TNN.GCmean)
estceiling = function(x) x["R.N","estimate"]
E.N=lapply(region.N(TNN.GCmean),estceiling) 
# Set the expected number of snow leopards in the region
scale = 1 # amoun by which to multiply sample size
scalefn = function(x,scale) x*scale
simD = lapply(TNN.GC.D,scalefn,scale=scale)
Dmean = function(x) mean(covariates(x)$D.0)
meanD = lapply(simD,Dmean)
nsim=100
Nest = cover = list(Tost=rep(NA,nsim), Noyon=rep(NA,nsim),Nemegt=rep(NA,nsim))
# create progress bar
pb <- tkProgressBar(title="progress bar", min=0, max=nsim, width=300)
seed=31 # for repeatability of this one population
for(i in 1:nsim) {
  if(i==1) {
    pop1 = sim.popn(simD,core=mask,poly=boundary,model2D="IHP",Ndbn="fixed",seed=seed)
  else   pop = sim.popn(simD,core=mask,poly=boundary,model2D="IHP",Ndbn="fixed")
  #  udist = geomLCdist(cams,pop,ne.fudge) # first set least-cost distances
  #  covariates(ne.fudge)$noneuc = covariates(ne.fudge)$noneuc.0 # create the noneuc variable needed for LC dist
  ch = sim.capthist(cams,pop,detectfn="HHN",detectpar=detpars,noccasions=1)
  simfit.0 <-secr.fit(ch, detectfn="HHN", mask=mask,  model=list(D~1, lambda0~1, sigma~1), 
                      start=list(D=meanD, lambda=detpars$lambda0, sigma~detpars$sigma), trace=0)
  Nest.0 = region.N(simfit.0)
  Nest[i] = Nest.0["R.N","estimate"]
  cover[i] = (Nest.0["R.N","lcl"]<=E.N & E.N<=Nest.0["R.N","ucl"])
  setTkProgressBar(pb, i, title=paste( round(i/nsim*100, 0),"% done"))
}
close(pb)
hist(Nest,nclass=30)
100*(mean(Nest)-E.N)/E.N
sum(cover)/nsim
sim100_N16_alpha0 = list(Nest=Nest,cover=cover,trueD=simD,seed=seed)
saveRDS(sim100_N16_alpha0,file="./Analysis4paper/sims/Tost_sim100_N16_alpha0.Rds")

