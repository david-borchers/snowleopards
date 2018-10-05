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
#----------------------- Load the RData objects made by Make_TNN_RData.r --------------------------
load("./Analysis4paper/TNN_boundaries.RData") # Tostboundary,Noyonboundary,Nemegtboundary
load("./Analysis4paper/TNN_masks.RData") # TostMask,NoyonMask,NemegtMask
load("./Analysis4paper/TNN_caphists.RData") # Tost_ch,Noyon_ch,Nemegt_sh,TNN_ch)
#----------------------- ----------------------------------------------- --------------------------
## Bundle them together in a single Rds file for convenience:
#Tost = list(capthist=Tost_ch,mask=TostMask,boundary=Tostboundary)
#Noyon = list(capthist=Noyon_ch,mask=NoyonMask,boundary=Noyonboundary)
#Nemegt = list(capthist=Nemegt_ch,mask=NemegtMask,boundary=Nemegtboundary)
#saveRDS(list(Tost=Tost,Noyon=Noyon,Nemegt=Nemegt),file="./Analysis4paper/TNN.Rds")
dat = readRDS(TNNdat, file="./Analysis4paper/TNN.Rds")
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

# Fit a LC distance model with flat D:
# ===================================
# (startvals from previous fit)
startvals = list(D=exp(-9.4515904),lambda0=exp(-4.2505931),sigma=exp(8.6914951),noneuc=0.3314881)
sl.ne <-secr.fit(ch, detectfn="HHN", mask=mask,
                 model=list(D~1, lambda0~1, sigma~1, noneuc~stdGC-1), 
                 details = list(userdist = geomLCdist),
                 start=startvals,
                 link = list(noneuc="identity"))
# Also fit a Euclidian distance model with varying D
# ==================================================
# (startvals from previous fit)
startvals = list(D=exp(-9.5227750),lambda0=exp(-4.4012107),sigma=exp(8.8538228))
sl.nuD <-secr.fit(ch, detectfn="HHN", mask=mask, start=startvals,
                  model=list(D~stdGC, lambda0~1, sigma~1))
nuD = predictDsurface(sl.nuD)
plotcovariate(nuD,"D.0",col=parula(40))
region.N(sl.nuD)
#estimate SE.estimate       lcl      ucl  n
#E.N 15.48028    4.158663  9.227452 25.97025 14
#R.N 15.48035    1.346921 14.323074 20.78309 14

# Simulate a population using above density surface from sl.nuD
# =============================================================
E.N=45 # Set the expected number of snow leopards in the region
scale = N/region.N(sl.nuD)["R.N","estimate"]
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
# (Fiddle with sl.fudge$fit$par[2] and sl.fudge$fit$par[3] to give right sort of sample size:)
# ============================================================================================
sl.ne$fit$par
sl.fudge = sl.ne
sl.fudge$fit$par[4] = 2
sl.fudge$fit$par[3] = sl.ne$fit$par[3] - log(3)
sl.fudge$fit$par[2] = sl.ne$fit$par[2] - log(15)
# Plot to see if you have the kind of habitat use that you're happy with:
ne.fudge = predictDsurface(sl.fudge,parameter="noneuc")
par(mfrow=c(1,1))
lcusageplot(sl.ne,mask=ne,col=parula(15))
lcpathplot(mask=ne.fudge,"geommean",col=parula(20),main="noneuc",what="image")


# Simulate capture historise
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

# Fit constant density model:
sl.ne.fit <-secr.fit(ch, detectfn="HHN", mask=mask,
                     model=list(D~1, lambda0~1, sigma~1, noneuc~stdGC-1), 
                     details = list(userdist = geomLCdist),
                     start=list(lambda0=exp(sl.fudge$fit$par[2]),
                                sigma=exp(sl.fudge$fit$par[3]),
                                noneuc=sl.fudge$fit$par[4]),
                     link = list(noneuc="identity"))
sl.noneuc = predictDsurface(sl.ne.fit,parameter="noneuc")
# Look at some usage and least-cost paths:
par(mfrow=c(1,1))
lcusageplot(sl.ne.fit,mask=sl.noneuc,col=parula(15))
lcpathplot(mask=sl.noneuc,"geommean",col=parula(20),main="noneuc",what="image")

# Fit varying density model with LC distance:
sl.ne.fit.D <-secr.fit(ch, detectfn="HHN", mask=mask,
                       model=list(D~stdGC, lambda0~1, sigma~1, noneuc~stdGC-1), 
                       details = list(userdist = geomLCdist),
                       start=list(lambda0=exp(sl.fudge$fit$par[2]),
                                  sigma=exp(sl.fudge$fit$par[3]),
                                  noneuc=sl.fudge$fit$par[4]),
                       link = list(noneuc="identity"))
sl.noneuc.D = predictDsurface(sl.ne.fit.D)
plotcovariate(sl.noneuc.D,"D.0",col=parula(40))
sl.noneuc.D = predictDsurface(sl.ne.fit,parameter="noneuc")
# Look at some usage and least-cost paths:
par(mfrow=c(1,1))
lcpathplot(mask=sl.noneuc.D,"geommean",col=parula(20),main="noneuc",what="image")
lcusageplot(sl.ne.fi.Dt,mask=sl.noneuc.D,col=parula(15))

# Fit a null model:
sl.ne.0 <-secr.fit(ch, detectfn="HHN", mask=mask,
                   model=list(D~1, lambda0~1, sigma~1))


# Fit varying density model with Euclidian distance:
sl.fit.D <-secr.fit(sl.ch, detectfn="HHN", mask=sl.mask,
                    model=list(D~stdGC, lambda0~1, sigma~1))
sl.D = predictDsurface(sl.fit.D)
plotcovariate(sl.noneuc.D,"D.0",col=parula(40))

# Compare AICs:
AIC(sl.ne.0,sl.ne.fit,sl.ne.fit.D,sl.fit.D)
# model          detectfn npar    logLik     AIC    AICc  dAICc AICcwt
# sl.ne.fit       D~1 lambda0~1 sigma~1 noneuc~stdGC - 1 hazard halfnormal    4 -136.8904 281.781 284.003  0.000 0.8323
# sl.ne.fit.D D~stdGC lambda0~1 sigma~1 noneuc~stdGC - 1 hazard halfnormal    5 -136.8387 283.677 287.207  3.204 0.1677
# sl.fit.D                     D~stdGC lambda0~1 sigma~1 hazard halfnormal    4 -143.4723 294.945 297.167 13.164 0.0000
# sl.ne.0                          D~1 lambda0~1 sigma~1 hazard halfnormal    3 -145.4336 296.867 298.130 14.127 0.0000


region.N(sl.ne.fit)
# estimate SE.estimate      lcl      ucl  n
# E.N 52.80341    13.46921 32.28109 86.37253 23
# R.N 45.78204    11.34091 32.05916 80.29244 23
region.N(sl.ne.0)
# estimate SE.estimate      lcl      ucl  n
# E.N 28.55434    6.512775 18.36467 44.39775 23
# R.N 28.55443    3.723157 24.68312 41.33008 23
region.N(sl.ne.fit.D)
# estimate SE.estimate      lcl      ucl  n
# E.N 44.07274    23.20644 16.71454 116.2105 23
# R.N 39.29320    22.23659 25.18305 144.6044 23
region.N(sl.fit.D)
# estimate SE.estimate      lcl      ucl  n
# E.N 27.25502    5.919049 17.89428 41.51247 23
# R.N 27.25503    2.789288 24.31818 36.73504 23

