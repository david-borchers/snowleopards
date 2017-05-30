library(secr)
library(fields)
source("scrplotting.r")

# Best models from Koustubh (load these using code in "david's summary stuff.r"):
Tost.hhn.DHab.nonU.GBx
Noyon.hhn.DHab.nonU.GBx
Nemegt.hhn.DHab.nonU.LamW2xR

# Get cameras, capture histories and masks out of these
# Capture histories
ch.Tost = Tost.hhn.DHab.nonU.GBx$capthist
ch.Noyon = Noyon.hhn.DHab.nonU.GBx$capthist
ch.Nemegt = Nemegt.hhn.DHab.nonU.LamW2xR$capthist
# Masks
mask.Tost = Tost.hhn.DHab.nonU.GBx$mask
mask.Noyon = Noyon.hhn.DHab.nonU.GBx$mask
mask.Nemegt = Nemegt.hhn.DHab.nonU.LamW2xR$mask
# Cameras
cams.Tost = traps(Tost.hhn.DHab.nonU.GBx$capthist)
cams.Noyon = traps(Noyon.hhn.DHab.nonU.GBx$capthist)
cams.Nemegt = traps(Nemegt.hhn.DHab.nonU.LamW2xR$capthist)

# But Nemegt is better if use stdBC than is use stdGC:
Nemegt.hhn.DHab.nonUB.LamW2xR<-secr.fit(ch.Nemegt, detectfn="HHN", mask=mask.Nemegt,
                                       model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdBC -1), 
                                       details = list(userdist = userdfn1),
                                       start = list(noneuc = 1))
coefficients(Nemegt.hhn.DHab.nonUB.LamW2xR)
AIC(Nemegt.hhn.DHab.nonU.LamW2xR,Nemegt.hhn.DHab.nonUB.LamW2xR)

# Best models have Nemegt with stdBC instead of stdGC:
Tost.hhn.DHab.nonU.GBx
Noyon.hhn.DHab.nonU.GBx
Nemegt.hhn.DHab.nonUB.LamW2xR

# Predict density, log(density) and noneuc, and put into covariates of Dhat.<region>
Dhat.Tost = predictDsurface(Tost.hhn.DHab.nonU.GBx)
covariates(Dhat.Tost)$logD = log(covariates(Dhat.Tost)$D.0)
covariates(Dhat.Tost)$noneuc =  covariates(predictDsurface(Tost.hhn.DHab.nonU.GBx,parameter="noneuc"))$noneuc.0
Dhat.Noyon = predictDsurface(Noyon.hhn.DHab.nonU.GBx)
covariates(Dhat.Noyon)$logD = log(covariates(Dhat.Noyon)$D.0)
covariates(Dhat.Noyon)$noneuc =  covariates(predictDsurface(Noyon.hhn.DHab.nonU.GBx,parameter="noneuc"))$noneuc.0
Dhat.Nemegt = predictDsurface(Nemegt.hhn.DHab.nonUB.LamW2xR)
covariates(Dhat.Nemegt)$logD = log(covariates(Dhat.Nemegt)$D.0)
covariates(Dhat.Nemegt)$noneuc =  covariates(predictDsurface(Nemegt.hhn.DHab.nonUB.LamW2xR,parameter="noneuc"))$noneuc.0



# Plotting all regions together:
# =============================
# get x and y ranges by finding extent of bounding boxes:
# Read boundary files
boundaryNemegt=readShapeSpatial("./Nemegt//Habitat/Nemegt_StudyArea.shp")
boundaryNoyon=readShapeSpatial("./Noyon2013/Habitat/NoyonStudy_Area.shp")
boundaryTost=readShapeSpatial("./Tost//Habitat/TostStudy_Area.shp")
bbox.Nemegt = bbox(boundaryNemegt)
bbox.Noyon = bbox(boundaryNoyon)
bbox.Tost = bbox(boundaryTost)
bbxlim = range(bbox.Nemegt["x",],bbox.Noyon["x",],bbox.Tost["x",]) 
bbylim = range(bbox.Nemegt["y",],bbox.Noyon["y",],bbox.Tost["y",])


# stdGC Plots
# -----------
# Color scale for same stdGC color across regions:
GCrange = range(covariates(Dhat.Tost)$stdGC,covariates(Dhat.Noyon)$stdGC,covariates(Dhat.Nemegt)$stdGC)
nGC = 40
GCcol = terrain.colors(nGC)
GCbreaks = seq(GCrange[1],GCrange[2],length=nGC+1)

# Do the plots of stdGC, with capture locations and traps
pdf(file="capthists.pdf",h=3.75,w=9)
par(mar=c(1, 1, 0, 4) + 0.1)
plot(bbxlim,bbylim,xlim=bbxlim,ylim=bbylim,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",asp=1) 
plotcovariate(Dhat.Tost, covariate='stdGC',xaxt="n",yaxt="n",xlab="",ylab="",bty="n",contour=FALSE,add=TRUE,asp=1,key=FALSE,col=GCcol,breaks=GCbreaks)
plotcovariate(Dhat.Noyon, covariate='stdGC',xaxt="n",yaxt="n",xlab="",ylab="",bty="n",contour=FALSE,add=TRUE,asp=1,key=FALSE,col=GCcol,breaks=GCbreaks)
plotcovariate(Dhat.Nemegt, covariate='stdGC',xaxt="n",yaxt="n",xlab="",ylab="",bty="n",contour=FALSE,add=TRUE,asp=1,key=TRUE,col=GCcol,breaks=GCbreaks)
#plot(boundaryNemegt,add=TRUE,border=1)
#plot(boundaryNoyon,add=TRUE,border=1)
#plot(boundaryTost,add=TRUE,border=1)
plotMaskEdge(Tost.hhn.DHab.nonU.GBx$mask,add=TRUE)
plotMaskEdge(Noyon.hhn.DHab.nonU.GBx$mask,add=TRUE)
plotMaskEdge(Nemegt.hhn.DHab.nonUB.LamW2xR$mask,add=TRUE)
plot(ch.Tost,tracks=TRUE,add=TRUE,varycol=FALSE,trkpar=list(col="black"),cappar=list(col="black",cex=1))
plot(ch.Noyon,tracks=TRUE,add=TRUE,varycol=FALSE,trkpar=list(col="black"),cappar=list(col="black",cex=1))
plot(ch.Nemegt,tracks=TRUE,add=TRUE,varycol=FALSE,trkpar=list(col="black"),cappar=list(col="black",cex=1))
plot(cams.Tost,add=TRUE,detpar=list(col="red"))
plot(cams.Noyon,add=TRUE,detpar=list(col="red"))
plot(cams.Nemegt,add=TRUE,detpar=list(col="red"))
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
dev.off()




# stdBC Plots
# -----------
# Color scale for sane stdBC color across regions:
BCrange = range(covariates(Dhat.Tost)$stdBC,covariates(Dhat.Noyon)$stdBC,covariates(Dhat.Nemegt)$stdBC)
nBC = 8
BCcol = terrain.colors(nBC)
BCbreaks = seq(BCrange[1],BCrange[2],length=nBC+1)

# Do the plots of stdGC, with capture locations and traps
pdf(file="BCs.pdf",h=3.75,w=9)
par(mar=c(1, 1, 0, 4) + 0.1)
plot(bbxlim,bbylim,xlim=bbxlim,ylim=bbylim,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",asp=1) 
plotcovariate(Dhat.Tost, covariate='stdBC',xaxt="n",yaxt="n",xlab="",ylab="",bty="n",contour=FALSE,add=TRUE,asp=1,key=FALSE,col=BCcol,breaks=BCbreaks)
plotcovariate(Dhat.Noyon, covariate='stdBC',xaxt="n",yaxt="n",xlab="",ylab="",bty="n",contour=FALSE,add=TRUE,asp=1,key=FALSE,col=BCcol,breaks=BCbreaks)
plotcovariate(Dhat.Nemegt, covariate='stdBC',xaxt="n",yaxt="n",xlab="",ylab="",bty="n",contour=FALSE,add=TRUE,asp=1,key=TRUE,col=BCcol,breaks=BCbreaks)
plotMaskEdge(Tost.hhn.DHab.nonU.GBx$mask,add=TRUE)
plotMaskEdge(Noyon.hhn.DHab.nonU.GBx$mask,add=TRUE)
plotMaskEdge(Nemegt.hhn.DHab.nonUB.LamW2xR$mask,add=TRUE)
plot(ch.Tost,tracks=TRUE,add=TRUE,varycol=FALSE,trkpar=list(col="black"),cappar=list(col="black",cex=1))
plot(ch.Noyon,tracks=TRUE,add=TRUE,varycol=FALSE,trkpar=list(col="black"),cappar=list(col="black",cex=1))
plot(ch.Nemegt,tracks=TRUE,add=TRUE,varycol=FALSE,trkpar=list(col="black"),cappar=list(col="black",cex=1))
plot(cams.Tost,add=TRUE,detpar=list(col="red"))
plot(cams.Noyon,add=TRUE,detpar=list(col="red"))
plot(cams.Nemegt,add=TRUE,detpar=list(col="red"))
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
dev.off()



# noneuc Plots
# -----------
# Color scale for uniform noneuc color across regions:
NUrange = range(covariates(Dhat.Tost)$stdBC,covariates(Dhat.Noyon)$stdBC,covariates(Dhat.Nemegt)$stdBC)
nNU = 40
NUcol = terrain.colors(nNU)
NUbreaks = seq(NUrange[1],NUrange[2],length=nNU+1)

# Do the plots of stdGC, with capture locations and traps
pdf(file="noneucs.pdf",h=3.75,w=9)
par(mar=c(1, 1, 0, 4) + 0.1)
plot(bbxlim,bbylim,xlim=bbxlim,ylim=bbylim,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n",asp=1) 
plotcovariate(Dhat.Tost, covariate='noneuc',xaxt="n",yaxt="n",xlab="",ylab="",bty="n",contour=FALSE,add=TRUE,asp=1,key=FALSE,col=NUcol,breaks=NUbreaks)
plotcovariate(Dhat.Noyon, covariate='noneuc',xaxt="n",yaxt="n",xlab="",ylab="",bty="n",contour=FALSE,add=TRUE,asp=1,key=FALSE,col=NUcol,breaks=NUbreaks)
plotcovariate(Dhat.Nemegt, covariate='noneuc',xaxt="n",yaxt="n",xlab="",ylab="",bty="n",contour=FALSE,add=TRUE,asp=1,key=TRUE,col=NUcol,breaks=NUbreaks)
plotMaskEdge(Tost.hhn.DHab.nonU.GBx$mask,add=TRUE)
plotMaskEdge(Noyon.hhn.DHab.nonU.GBx$mask,add=TRUE)
plotMaskEdge(Nemegt.hhn.DHab.nonUB.LamW2xR$mask,add=TRUE)
plot(ch.Tost,tracks=TRUE,add=TRUE,varycol=FALSE,trkpar=list(col="black"),cappar=list(col="black",cex=1))
plot(ch.Noyon,tracks=TRUE,add=TRUE,varycol=FALSE,trkpar=list(col="black"),cappar=list(col="black",cex=1))
plot(ch.Nemegt,tracks=TRUE,add=TRUE,varycol=FALSE,trkpar=list(col="black"),cappar=list(col="black",cex=1))
plot(cams.Tost,add=TRUE,detpar=list(col="red"))
plot(cams.Noyon,add=TRUE,detpar=list(col="red"))
plot(cams.Nemegt,add=TRUE,detpar=list(col="red"))
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
dev.off()



# Some investigation of the LC paths

lcpathplot(Dhat.Nemegt,"stdBC",transitionFunction=function(x) 1/x,n=1,contour=FALSE,col=c("orange","brown"),key=FALSE, linecol="white",
           xaxt="n",yaxt="n",xlab="",ylab="",bty="n",asp=1)

covariate = "stdBC"
noneucdists(Dhat.Tost,userdfn1,covariate)
noneucdists(Dhat.Noyon,userdfn1,covariate)
noneucdists(Dhat.Nemegt,userdfn1,covariate)






unique(covariates(Dhat.Tost)$BINCODE)
unique(covariates(Dhat.Noyon)$BINCODE)
unique(covariates(Dhat.Nemegt)$BINCODE)
unique(covariates(Dhat.Tost)$stdBC)
unique(covariates(Dhat.Noyon)$stdBC)
unique(covariates(Dhat.Nemegt)$stdBC)
diff(range(unique(covariates(Dhat.Tost)$stdBC)))
diff(range(unique(covariates(Dhat.Noyon)$stdBC)))
diff(range(unique(covariates(Dhat.Nemegt)$stdBC)))

# Try fitting best models with BINCODE instead of stdBC:
Nemegt.hhn.DHab.nonUBIN.LamW2xR<-secr.fit(ch.Nemegt, detectfn="HHN", mask=mask.Nemegt,
                                        model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ BINCODE -1), 
                                        details = list(userdist = userdfn1),
                                        start = list(noneuc = 1))
coefficients(Nemegt.hhn.DHab.nonUBIN.LamW2xR)
AIC(Nemegt.hhn.DHab.nonU.LamW2xR,Nemegt.hhn.DHab.nonUB.LamW2xR,Nemegt.hhn.DHab.nonUBIN.LamW2xR)



