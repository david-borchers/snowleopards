# ================== MODELS WITHOUT NONEUC ===================================

library(secr)
library(fields)
library(maptools)
source("scrplotting.r")

#----------------------- Load the RData objects made by Make_TNN_RData.r --------------------------
load("./Analysis4paper/TNN_boundaries.RData") # Tostboundary,Noyonboundary,Nemegtboundary
load("./Analysis4paper/TNN_masks.RData") # TostMask,NoyonMask,NemegtMask
#load("./Analysis4paper/TN_caphists.RData") # Tost_ch,Noyon_ch,TN_ch
load("./Analysis4paper/TNN_caphists.RData") # Tost_ch,Noyon_ch,Nemegt_sh,TNN_ch
#----------------------- ----------------------------------------------- --------------------------

#----------------------- Do some plots to check data seems OK --------------------------

# Find extent of bounding boxes of boundarys:
bbox.Nemegt = bbox(Nemegtboundary)
bbox.Noyon = bbox(Noyonboundary)
bbox.Tost = bbox(Tostboundary)
bbxlim = range(bbox.Noyon["x",],bbox.Tost["x",],bbox.Nemegt["x",]) #Set plot limit for all 3 areas together
bbylim = range(bbox.Noyon["y",],bbox.Tost["y",],bbox.Nemegt["y",])

dxlim = diff(bbxlim)
xlim = c(bbxlim[1],bbxlim[2]+0.15*dxlim)
ylim = bbylim


if(.Platform$OS.type=="windows") { # this to make the command quartz( ) work on windows machines
  quartz<-function() windows()
}
quartz(w=8,h=4)

# Plot boundaries
plot(xlim,ylim,xlim=xlim,ylim=ylim,type="n",asp=1,bty="n",xlab="Easting",ylab="Northing") 
plot(Tostboundary,add=TRUE)
plot(Noyonboundary,add=TRUE)
plot(Nemegtboundary,add=TRUE)

# Plot GC means
pdf("./ANalysis4paper/Plots/TNNrmeanGC.pdf",h=5,w=10)
zlim = range(covariates(NemegtMask)$rmeanGC,
             covariates(NoyonMask)$rmeanGC,
             covariates(TostMask)$rmeanGC)
plot(xlim,ylim,xlim=xlim,ylim=ylim,type="n",asp=1,bty="n",xlab="Easting",ylab="Northing") 
plotcovariate(NemegtMask,covariate="rmeanGC",add=TRUE,zlim=zlim,contour=FALSE)
plotcovariate(NoyonMask,covariate="rmeanGC",add=TRUE,zlim=zlim,contour=FALSE)
plotcovariate(TostMask,covariate="rmeanGC",add=TRUE,zlim=zlim,contour=FALSE)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
dev.off()

# Plot GC sums
pdf("./ANalysis4paper/Plots/TNNrsumGC.pdf",h=5,w=10)
zlim = range(covariates(NemegtMask)$rsumGC,
             covariates(NoyonMask)$rsumGC,
             covariates(TostMask)$rsumGC)
plot(xlim,ylim,xlim=xlim,ylim=ylim,type="n",asp=1,bty="n",xlab="Easting",ylab="Northing") 
plotcovariate(NemegtMask,covariate="rsumGC",add=TRUE,zlim=zlim,contour=FALSE)
plotcovariate(NoyonMask,covariate="rsumGC",add=TRUE,zlim=zlim,contour=FALSE)
plotcovariate(TostMask,covariate="rsumGC",add=TRUE,zlim=zlim,contour=FALSE)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
dev.off()

# Plot deviations from GC mean
pdf("./ANalysis4paper/Plots/TNNstdGCdev.pdf",h=5,w=10)
zlim = range(covariates(NemegtMask)$rmeanGCdev,
             covariates(NoyonMask)$rmeanGCdev,
             covariates(TostMask)$rmeanGCdev)
plot(xlim,ylim,xlim=xlim,ylim=ylim,type="n",asp=1,bty="n",xlab="Easting",ylab="Northing") 
plotcovariate(NemegtMask,covariate="rmeanGCdev",add=TRUE,zlim=zlim,contour=FALSE,col=terrain.colors(40))
plotcovariate(NoyonMask,covariate="rmeanGCdev",add=TRUE,zlim=zlim,contour=FALSE,col=terrain.colors(40))
plotcovariate(TostMask,covariate="rmeanGCdev",add=TRUE,zlim=zlim,contour=FALSE,col=terrain.colors(40))
# add cameras
plot(TNN.cams,detpar=list(col="black"),add=TRUE)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
dev.off()

quartz(h=6,w=3)
par(mfrow=c(3,1))
hist(covariates(NemegtMask)$stdGC)
hist(covariates(NoyonMask)$stdGC)
hist(covariates(TostMask)$stdGC)
par(mfrow=c(1,1))


zlim = range(covariates(NoyonMask)$GC,
             covariates(TostMask)$GC)
pdf("./ANalysis4paper/Plots/TNgridcode.pdf",h=5,w=10)
plot(bbxlim,bbylim,xlim=xlim,ylim=ylim,xlab="",ylab="",bty="n",type="n",xaxt="n",yaxt="n",asp=1) 
# Plot the terrain
plotcovariate(NoyonMask, covariate="GC", contour = FALSE, col = terrain.colors(16), zlim=zlim, add = TRUE)
plotcovariate(TostMask, covariate="GC", contour=FALSE, col=terrain.colors(16), zlim=zlim, add=TRUE)
# Add the traps
plot(traps(TN_ch)[[1]],add=TRUE,detpar=list(col="black",pch="+"))
plot(traps(TN_ch)[[2]],add=TRUE,detpar=list(col="black",pch="+"))

# and the borders over the top
plot(Noyonboundary,add=TRUE)#,border=2)
plot(Tostboundary,add=TRUE)#,border=1)
# add region labels
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
dev.off()

# ----------------------- End Plotting ----------------------------


# ===================== Tost and Noyon only ====================
TN_ch = TNN_ch[1:2] # extract Tost and Noyon capture histories only
# --------------- Fits to indiviual regions --------------------
# Models like this did not converge!
TN.Noyon.stdGC<-secr.fit(Noyon_ch, detectfn="HHN", mask=NoyonMask,
                         model=list(D~stdGC, lambda0~1, sigma~stdGC))
# ... so I removed stdGC from sigma:
TN.Tost.stdGC<-secr.fit(Tost_ch, detectfn="HHN", mask=TostMask,
                        model=list(D~stdGC, lambda0~1, sigma~1))
TN.Noyon.stdGC<-secr.fit(Noyon_ch, detectfn="HHN", mask=NoyonMask,
                         model=list(D~stdGC, lambda0~1, sigma~1))
coefficients(TN.Tost.stdGC)
coefficients(TN.Noyon.stdGC)

# Try same models, but with centred stdGC instead of intercept:
# first just check that stdGC = rmeanGC + rmeanGCdev
stdGCs = covariates(TostMask)$stdGC
rmeanGCs = covariates(TostMask)$rmeanGC
rmeanGCdevs = covariates(TostMask)$rmeanGCdev
sum(stdGCs-rmeanGCs-rmeanGCdevs)
# Now fit same model, but with rmeanGCdev
TN.Tost.stdGCmean<-secr.fit(Tost_ch, detectfn="HHN", mask=TostMask, model=list(D~rmeanGCdev, lambda0~1, sigma ~ 1))
# check intersept term is as expected
coefficients(TN.Tost.stdGCmean)["D","beta"]
as.numeric(coefficients(TN.Tost.stdGC)["D","beta"]) + as.numeric(D.rmeanGCdev*rmeanGCs[1])
# and slope is as expected:
as.numeric(coefficients(TN.Tost.stdGC)["D.stdGC","beta"]) 
as.numeric(coefficients(TN.Tost.stdGCmean)["D.rmeanGCdev","beta"]) 
# Same for Noyon:
TN.Noyon.stdGCmean<-secr.fit(Noyon_ch, detectfn="HHN", mask=NoyonMask, model=list(D~rmeanGCdev, lambda0~1, sigma ~ 1))


# ----------------------- Combined region models ---------------------------
# Try same models, fitting all at once, using rmeanGCdev:
TNfit <- secr.fit(TN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask),
                  model=list(D~rmeanGCdev*session, lambda0~session, sigma~session))
# Compare coefficients for Tost
coefficients(TNfit)[c(1,2,5,7),"beta"]
as.numeric(coefficients(TN.Tost.stdGCmean)[,"beta"])
# Compare coefficients for Noyon
coefficients(TNfit)[c(3,4,6,8),"beta"] + coefficients(TNfit)[c(1,2,5,7),"beta"] 
as.numeric(coefficients(TN.Noyon.stdGCmean)[,"beta"])
# Noyon coefficients different from TN.Noyon.stdGCmean, but TN.Noyon.stdGCmean looks rubbish anywhay.
# Check if AICs add up: 
AIC(TNfit)$AIC
AIC(TN.Noyon.stdGCmean)$AIC + AIC(TN.Tost.stdGCmean)$AIC
# They do
# CONCLUSION so far: get identical fits to Tost and Noyon by fitting separately, or combined, with interactions
# (This is no surprise, but worth checking before moving on.)

# Now try with
# (1) rmeanGCV to explain difference in D due to different average stdGC available in two regions, 
# (2) rmeanGCdev to explain effect of stdGC on D within region
TNfit.rmeanGC <- secr.fit(TN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask),
                  model=list(D~rmeanGC + rmeanGCdev:session, lambda0~session, sigma~session))
# Check if AICs same for both combined models: 
AIC(TNfit)$AIC; AIC(TNfit.rmeanGC)$AIC
# Check that coefficients make sense"
# Lambda and sigma parameters:
coefficients(TNfit)[5:8,"beta"] 
coefficients(TNfit.rmeanGC)[5:8,"beta"]
# Check predicted density
Nhat = region.N(TNfit)
Nhat.rmeanGC = region.N(TNfit.rmeanGC)
Nhat; Nhat.rmeanGC
# OK, so far so good

# Now add 
# (3) session factor to explain difference in D over and above (1) and (2) - which might be conservation status effect
TNfit.rmeanGC.sess <- secr.fit(TN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask),
                          model=list(D~rmeanGC + rmeanGCdev:session + session, lambda0~session, sigma~session))

# Compare models with and without the region effect, over and above the GC effects:
AIC(TNfit.rmeanGC.sess,TNfit.rmeanGC)
# CONCLUSION: Model with no session effect is preferable (zero improvement in likelihood), 
# i.e. differenc in density between Tost and Noyon is explained by different mean(GC) in each region



# ===================== Tost, Noyon and Nemegt ====================
# --------------- Fits to indiviual regions --------------------

# No need to refit to Tost and Noyon, unless did not fit with code above (in which case, remove # from next 2 lines:)
#TN.Tost.stdGCmean<-secr.fit(Tost_ch, detectfn="HHN", mask=TostMask, model=list(D~rmeanGCdev, lambda0~1, sigma ~ 1))
#TN.Noyon.stdGCmean<-secr.fit(Noyon_ch, detectfn="HHN", mask=NoyonMask, model=list(D~rmeanGCdev, lambda0~1, sigma ~ 1))
TN.Nemegt.stdGCmean<-secr.fit(Nemegt_ch, detectfn="HHN", mask=NemegtMask, model=list(D~rmeanGCdev, lambda0~1, sigma ~ 1))


# ----------------------- Combined region models ---------------------------
# Try same models, fitting all at once, using rmeanGCdev:
TNNfit <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                  model=list(D~rmeanGCdev*session, lambda0~session, sigma~session))
coefficients(TNNfit)
# Check if AICs add up: 
AIC(TNfit)$AIC
AIC(TN.Noyon.stdGCmean)$AIC + AIC(TN.Tost.stdGCmean)$AIC + AIC(TN.Nemegt.stdGCmean)$AIC
# They do
# CONCLUSION so far: get identical fits to Tost and Noyon by fitting separately, or combined, with interactions
# (This is no surprise, but worth checking before moving on.)

# Now try with
# (1) rmeanGCV to explain difference in D due to different average stdGC available in three regions, 
# (2) rmeanGCdev to explain effect of stdGC on D within region
TNNfit.rmeanGC <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                          model=list(D~rmeanGC + rmeanGCdev:session, lambda0~session, sigma~session))
coefficients(TNNfit.rmeanGC)
coefficients(TNNfit)
# Check AICs: TNNfit has one more parameter (one for each session, compared to intercept & slope for rmeanGC)
AIC(TNNfit,TNNfit.rmeanGC)
# Check predicted density
TNN.Nhat = region.N(TNNfit)
TNN.Nhat.rmeanGC = region.N(TNNfit.rmeanGC)
TNN.Nhat; TNN.Nhat.rmeanGC

# Now add 
# (3) session factor to explain difference in D over and above (1) and (2) - which might be conservation status effect
TNNfit.rmeanGC.sess <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                               model=list(D~rmeanGC + rmeanGCdev:session + session, lambda0~session, sigma~session))

# Compare models with and without the region effect, over and above the GC effects:
AIC(TNNfit.rmeanGC.sess,TNNfit.rmeanGC)
# CONCLUSION: No conservation status effect

# Save the fits:
save(TN.Tost.stdGCmean,TN.Noyon.stdGCmean,TN.Nemegt.stdGCmean,TNfit,TNNfit,TNNfit.rmeanGC.sess,TNNfit.rmeanGC,file="./Analysis4paper/TNNnoNEfits.RData")

# Predict density surface from best combined model:
predD = predictDsurface(TNNfit.rmeanGC)
covariates(TostMask)$Dhat = covariates(predD[[1]])$D.0
covariates(NoyonMask)$Dhat = covariates(predD[[2]])$D.0
covariates(NemegtMask)$Dhat = covariates(predD[[3]])$D.0

# Plot Density estimates
pdf("./ANalysis4paper/Plots/TNNDhats.pdf",h=5,w=10)
zlim = range(covariates(NemegtMask)$Dhat,
             covariates(NoyonMask)$Dhat,
             covariates(TostMask)$Dhat)
plot(xlim,ylim,xlim=xlim,ylim=ylim,type="n",asp=1,bty="n",xlab="Easting",ylab="Northing") 
plotcovariate(NemegtMask,covariate="Dhat",add=TRUE,zlim=zlim,contour=FALSE)
plotcovariate(NoyonMask,covariate="Dhat",add=TRUE,zlim=zlim,contour=FALSE)
plotcovariate(TostMask,covariate="Dhat",add=TRUE,zlim=zlim,contour=FALSE)
# add cameras
#plot(TNN.cams,detpar=list(col="red",cex=0.5),add=TRUE)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
dev.off()



# Predict density surface from separate models:
covariates(TostMask)$Dhat.ind = covariates(predictDsurface(TN.Tost.stdGCmean))$D.0
covariates(NoyonMask)$Dhat.ind = covariates(predictDsurface(TN.Noyon.stdGCmean))$D.0
covariates(NemegtMask)$Dhat.ind = covariates(predictDsurface(TN.Nemegt.stdGCmean))$D.0

# Plot Density estimates
pdf("./ANalysis4paper/Plots/TNNDhatsIndividModels.pdf",h=5,w=10)
zlim = range(covariates(NemegtMask)$Dhat.ind,
             covariates(NoyonMask)$Dhat.ind,
             covariates(TostMask)$Dhat.ind)
plot(xlim,ylim,xlim=xlim,ylim=ylim,type="n",asp=1,bty="n",xlab="Easting",ylab="Northing") 
plotcovariate(NemegtMask,covariate="Dhat.ind",add=TRUE,zlim=zlim,contour=FALSE)
plotcovariate(NoyonMask,covariate="Dhat.ind",add=TRUE,zlim=zlim,contour=FALSE)
plotcovariate(TostMask,covariate="Dhat.ind",add=TRUE,zlim=zlim,contour=FALSE)
# add cameras
#plot(TNN.cams,detpar=list(col="red",cex=0.5),add=TRUE)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
dev.off()

