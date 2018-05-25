# ================== MODELS WITHOUT NONEUC ===================================

library(secr)
library(fields)
library(maptools)
source("scrplotting.r")
library(parallel)

detectCores()

#----------------------- Load the RData objects made by Make_TNN_RData.r --------------------------
load("./Analysis4paper/TNN_boundaries.RData") # Tostboundary,Noyonboundary,Nemegtboundary
load("./Analysis4paper/TNN_masks.RData") # TostMask,NoyonMask,NemegtMask
#load("./Analysis4paper/TN_caphists.RData") # Tost_ch,Noyon_ch,TN_ch
load("./Analysis4paper/TNN_caphists.RData") # Tost_ch,Noyon_ch,Nemegt_sh,TNN_ch)
#----------------------- ----------------------------------------------- --------------------------

TNN.trapfiles = c("./Tost_Noyon_Nemegt/Tost_Traps.txt",
                  "./Tost_Noyon_Nemegt/Noyon_Traps.txt",
                  "./Tost_Noyon_Nemegt/Nemegt_Traps.txt")
all.data.TNN<-read.capthist(captfile = "./Tost_Noyon_Nemegt/TNN_Capture.csv", 
                            binary.usage = FALSE, trapfile = TNN.trapfiles, 
                            detector="count", fmt = "trapID", 
                            trapcovnames = c("Rgd","Topo", "Water", "Winter"))
summary(TNN_ch)


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
pdf("./ANalysis4paper/Plots/TNNCameras.pdf",h=5,w=10)
plot(xlim,ylim,xlim=xlim,ylim=ylim,type="n",asp=1,bty="n",xlab="Easting",ylab="Northing") 
plot(Tostboundary,add=TRUE)
plot(Noyonboundary,add=TRUE)
plot(Nemegtboundary,add=TRUE)
plot(x=Tost.cams, add=TRUE)
plot(x=Noyon.cams, add=TRUE)
plot(x=Nemegt.cams, add=TRUE)


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
summary(covariates(traps(Nemegt_ch)))
summary(covariates(traps(Tost_ch)))
summary(covariates(traps(Noyon_ch)))

summary(covariates(traps(TNN_ch)))
names(covariates(TostMask))
# ----------------------- Combined region models ---------------------------
# Try same models, fitting all at once, using rmeanGCdev:
#Linear regression on log scale. Intercept depends on rmeanGC. Different intercept for each region because each of them have a different
#value of rmeanGC (average GC across whole of study area). At any point within the region, value of stdGC is the mean+deviation from the mean
#rmeanGC is the mean for the region. Rmeangcdev is the difference between the mean and value at that point.
#Each region will have a different slope!
#  X*session is the same as X+x:session
#TNNfit <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),binary.usage = FALSE,
#                  model=list(D~rmeanGC + rmeanGCdev:session, lambda0~session, sigma~session)
#coefficients(TNNfit)

#1  Density a function of ruggedness and varies differently (has a different slope) for different areas. Lambda0 and sigma different, but constant for each area
TNNfit.DGridSess.lamSess.SigSess <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                           model=list(D~rmeanGC + rmeanGCdev:session + session, lambda0~session, sigma~session), ncores=16)

#2  Density a function of ruggedness and varies differently for different areas. Lambda0 a function of topography but constant across areas
#   Sigma constant across all areas
TNN.DGrid.LamTopo.Sig1<-secr.fit(TNN_ch, detectfn="HHN", mask = list(TostMask, NoyonMask, NemegtMask),
                            model = list(D~rmeanGC + rmeanGCdev: session + session, lambda0~Topo, sigma~1), ncores=16)

#3  Density function of ruggedness and varies differently for 3 areas. Lambda0 a function of topography and different for 3 areas. Sigma constant
TNN.DGridXsess.LamTopoXsess<-secr.fit(TNN_ch, detectfn="HHN", mask = list(TostMask, NoyonMask, NemegtMask),
                            model = list(D~rmeanGC + rmeanGCdev: session + session, lambda0~Topo*session, sigma~1), ncores=16)

#4  Density function of ruggedness and varies differently for 3 areas. Lambda0 function of topography & different for 3 areas.
#   Sigma constant, but different for the three regions
TNN.DGridXsess.DetTopoXsess<-secr.fit(TNN_ch, detectfn="HHN", mask = list(TostMask, NoyonMask, NemegtMask),
                                      model = list(D~rmeanGC + rmeanGCdev: session + session, lambda0~Topo*session, sigma~session), ncores=16)

#5  Density constant and same across 3 areas. Lambda0 different for 3 areas, sigma constant and same for 3 areas
TNNfit.D1.LamSess.sig1 <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                   model=list(D~1, lambda0~session, sigma~1), ncores=16)

#6  Density, lambda0 and sigma constant and same across 3 areas. 
TNNfit.D1.Lam1s.sig1 <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                   model=list(D~1, lambda0~1, sigma~1), ncores=16)

### Temporary test of par.secr.fit
##  MULTIPLE CORES for several models listed together using par.secr.fit
#   IT WORKS!!!
#TNNfit.D1.LamSess.sig1000 <- list(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
#                                   model=list(D~1, lambda0~session, sigma~1))
#TNNfit.D1.Lam1s.sig1000 <- list(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
#                                 model=list(D~1, lambda0~1, sigma~1))
#fits000<-par.secr.fit (c('TNNfit.D1.LamSess.sig1000', 'TNNfit.D1.Lam1s.sig1000'), ncores = 4)


#7  Density, lambda0 and sigma constant but vary between three areas.
TNNfit.Dsess.Lamsess.sigsess <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                                 model=list(D~session, lambda0~session, sigma~session), ncores=16)

#8  Density varies as a function of ruggedness but is the same across three areas. a0 a function of topography and sigma is constant 
###   can this not simply use D~stdGC?
TNNfit.DGrid.a0Topo.sig1<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                    model=list(D~stdGC, a0~Topo, sigma~1), ncores=16)

#9  Density varies as a function of ruggedness and is the same for the three areas. a0 a function of topography and is different across 3 areas
#   Sigma different for the three study areas
###   can this not simply use D~stdGC?
TNNfit.DGrid.a0TopoSess.sig1<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                   model=list(D~stdGC, a0~Topo*session, sigma~1), ncores=16)

#10 Density a function of ruggedness, but not varying across 3 areas. a0 a function of topography and varies across 3 regions.
#   Sigma is constant but different for the three areas
### can  this not simply use D~stdGC?
TNNfit.DGrid.a0TopoSess.sigsess<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                       model=list(D~stdGC, a0~Topo*session, sigma~session), ncores=16)

#11 Density constant but varies between the three areas. a0 a function of topography and varies between 3 areas. Sigma different for each area
TNNfit.DSess.a0TopoSess.sigSess<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                          model=list(D~session, a0~Topo*session, sigma~session), ncores=16)

#12 Density constant across 3 areas. a0 varies between three areas and also a function of topography. Sigma constant
#Not on Thinkpad
TNNfit.D1.a0TopoSess.sig1<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                          model=list(D~1, a0~Topo*session, sigma~1), ncores=16)

#13 Density a function of ruggedness and varies differently in 3 areas. a0 varies with topography differently for 3 areas. Sigma is constant
# Not on Thinkpad
TNNfit.DGridSess.a0TopoSess.sig1<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                       model=list(D~rmeanGC + rmeanGCdev: session + session, a0~Topo*session, sigma~1), ncores=16)

#14 Density a function of ruggedness and varies differently between 3 areas. a0 varies with topography differently for 3 areas.
#   Sigma different for each area. NOT ON Thinkpad
TNNfit.DGridSess.a0TopoSess.sigsess<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                           model=list(D~rmeanGC + rmeanGCdev: session + session, a0~Topo*session, sigma~session), ncores=16)


#15 Density a function of ruggedness and varies differently between 3 areas. a0 different for three areas, sigma different for 3 areas.
# Not on Thinkpad
TNNfit.DGridSess.a0Sess.sigsess<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                              model=list(D~rmeanGC + rmeanGCdev: session + session, a0~session, sigma~session), ncores=16)

save(TNNfit.DGridSess.lamSess.SigSess, TNN.DGrid.LamTopo.Sig1, TNN.DGridXsess.LamTopoXsess,
     TNN.DGridXsess.DetTopoXsess, TNNfit.D1.LamSess.sig1, TNNfit.D1.Lam1s.sig1, TNNfit.Dsess.Lamsess.sigsess, 
     TNNfit.DSess.a0TopoSess.sigSess, TNNfit.DGrid_Sess.a0Topo.Sess.sigsess,
     file="./Analysis4paper/TNNnoNEfits_Thinkpad.RData")


#16 Density a function of ruggedness and varies with a different slope between 3 areas. a0 different for three areas, sigma different for 3 areas.
TNNfit.DGridXSess.a0Sess.sigsess<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                          model=list(D~stdGC*session, a0~session, sigma~session), ncores=16)

#17 Density a function of ruggedness and varies with a different slope between 3 areas. a0 different for three areas, sigma different for 3 areas.
TNNfit.DGrid_Sess.a0Sess.sigsess<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                           model=list(D~stdGC+session, a0~session, sigma~session), ncores=16)

#18 Density a function of ruggedness and varies with a different slope between 3 areas. a0 different for three areas, sigma different for 3 areas.
TNNfit.DGrid_Sess.a0TopoXSess.sigsess<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                           model=list(D~stdGC+session, a0~Topo*session, sigma~session), ncores=16)

#19 Density a function of ruggedness and varies with a different slope between 3 areas. a0 different for three areas, sigma different for 3 areas.
TNNfit.DGrid_Sess.a0Topo.Sess.sigsess<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                                model=list(D~stdGC+session, a0~session, sigma~session), ncores=16)



AIC(TNNfit.DGrid.a0TopoSess.sigsess,TNNfit.DGrid.a0TopoSess.sig1,TNNfit.DGrid.a0Topo.sig1, TNNfit.Dsess.Lamsess.sigsess,
    TNNfit.D1.Lam1s.sig1, TNNfit.D1.LamSess.sig1, TNN.DGridXsess.DetTopoXsess, TNN.DGridXsess.LamTopoXsess, TNN.DGrid.LamTopo,
    TNNfit, TNNfit.rmeanGC, TNNfit.DSess.a0TopoSess.sigSess, TNNfit.DGridXSess.a0Sess.sigsess, TNNfit.DGrid_Sess.a0Sess.sigsess,
    TNNfit.DGrid_Sess.a0Topo.Sess.sigsess)
# Save the fits:
save(TNNfit.DGrid.a0TopoSess.sigsess,TNNfit.DGrid.a0TopoSess.sig1,TNNfit.DGrid.a0Topo.sig1, TNNfit.Dsess.Lamsess.sigsess,
     TNNfit.D1.Lam1s.sig1, TNNfit.D1.LamSess.sig1, TNN.DGridXsess.DetTopoXsess, TNN.DGridXsess.LamTopoXsess, TNN.DGrid.LamTopo,
     TNNfit, TNNfit.rmeanGC, TNNfit.DSess.a0TopoSess.sigSess, TNNfit.DGridXSess.a0Sess.sigsess, TNNfit.DGrid_Sess.a0Sess.sigsess,
     TNNfit.DGrid_Sess.a0Topo.Sess.sigsess,file="./Analysis4paper/TNNnoNEfits_all.RData")

    
region.N(TNNfit.DGrid.a0Topo.sig1)
region.N(TNN.DGrid.LamTopo)

coefficients(TNNfit.DGrid.a0Topo.sig1)
coefficients(TNN.DGrid.LamTopo)

# Check if AICs add up: 
AIC(TNfit)$AIC
AIC(TN.Noyon.stdGCmean)$AIC + AIC(TN.Tost.stdGCmean)$AIC + AIC(TN.Nemegt.stdGCmean)$AIC
# They do
# CONCLUSION so far: get identical fits to Tost and Noyon by fitting separately, or combined, with interactions
# (This is no surprise, but worth checking before moving on.)

# Now try with
# (1) rmeanGCV to explain difference in D due to different average stdGC available in three regions, 
# (2) rmeanGCdev to explain effect of stdGC on D within region

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
# To use this one to estimate differences!
TNNfit.rmeanGC.sess <- secr.fit(TNN_ch, detectfn="HHN", mask=list(TostMask, NoyonMask, NemegtMask),
                               model=list(D~rmeanGC + rmeanGCdev:session + session, lambda0~session, sigma~session))
#This explains what's going on due to the difference in means. There could be a difference between the two is based on session!
#*session to define effect of covariate is not uniform

summary(covariates(traps(TNN_ch)))

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

load("./Analysis4paper/TNNnoNEfits_all.Rdata")
TNNSurface<-predictDsurface(TNNfit.DGrid.a0Waterxsess_topo.sig1, se.D=TRUE, cl.D=TRUE)
TostSurface<-predictDsurface(Tost.hhn.DHab, se.D=TRUE, cl.D=TRUE)
plot(TostSurface, contour = FALSE)

coefficients(TNNfit.DGrid.a0Waterxsess_topo.sig1)
coefficients(TNNfit.DGrid.a0Waterxsess.sig1)
  