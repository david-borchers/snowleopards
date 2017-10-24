NemegtADDITIONAL_ch<-read.capthist(captfile = "./Analysis4paper/Data/Nemegt_capthist.csv", 
                                   trapfile = "./Analysis4paper/Data/Nemegt_ADDITIONALcams.csv", 
                                   detector="count", binary.usage=FALSE, fmt = "trapID", 
                                   trapcovnames = c("Rgd", "Topo", "Water", "Winter"))
# get cameras
acams = traps(NemegtADDITIONAL_ch)
# add stdGC to cameras:
acams = addCovariates(acams,NemegtMask,columns="stdGC")
animalN = apply(NemegtADDITIONAL_ch[,1,],1,sum)
trapN = apply(NemegtADDITIONAL_ch[,1,],2,sum)
trap.zeros.rgd = covariates(acams)$stdGC[trapN==0]
trap.nonzeros.rgd = covariates(acams)$stdGC[trapN>0]
xrange = range(trap.zeros.rgd,trap.nonzeros.rgd)
dx = diff(xrange)/10
yrange = range(trapN)
par(mfrow=c(2,1))
plot(trap.zeros.rgd,trapN[trapN==0],xlim=xrange,ylim=yrange,xlab="stdGC",ylab="Number captures",pch="+",main="")
points(trap.nonzeros.rgd,trapN[trapN>0],col="red",pch="+")
hist(covariates(acams)$stdGC,breaks=seq(xrange[1],xrange[2],length=21),xlab="stdGC",main="")


nz.tinystdGC = which(trap.nonzeros.rgd<0.002 & covariates(acams)$stdGC[trapN>0])
tinystdGCcam = acams[nz.tinystdGC,]

# Plot Nemegt data with ADDITOINAL traps added, with camera that is at almost zero stdGC with a capture:
quartz(h=5,w=10)
plotcovariate(NemegtMask, covariate="GC", contour = FALSE, col = terrain.colors(16), zlim=zlim,asp=1)
plot(traps(NemegtADDITIONAL_ch),add=TRUE)
plot(NemegtADDITIONAL_ch,tracks=TRUE,add=TRUE)
points(tinystdGCcam$x,tinystdGCcam$y,pch="*",col="white")

rownames(tinystdGCcam)
tinytrap = which(rownames(acams)==rownames(tinystdGCcam))

keep = 1:length(rownames(acams))
keep = keep[-tinytrap]
testCH = subset(NemegtADDITIONAL_ch,traps=keep)

# Try Nemegt with best model but with ADDITIONAL traps in Steppe, but without trap NemegtSLT017:
NemegtADD.testCH<-secr.fit(testCH, detectfn="HHN", mask=NemegtMask,
                            model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                            details = list(userdist = userdfn1),
                            start = list(noneuc = 1))
# That did not change anything!

# Look at the trap with 14 out of the 50 detections
quartz(h=5,w=10)
plotcovariate(NemegtMask, covariate="GC", contour = FALSE, col = terrain.colors(16), zlim=zlim,asp=1)
plot(traps(NemegtADDITIONAL_ch),add=TRUE)
plot(NemegtADDITIONAL_ch,tracks=TRUE,add=TRUE)
Ncap = as.character(trapN)
text(acams$x,acams$y,labels=Ncap)

bigcapstrap = which(trapN==max(trapN))
NemegtADDITIONAL_ch[,1,bigcapstrap]

# Plot of Scrawny's capture history
cats = c("IdereGerle","Pundit","Scrawny","Sumbe")
quartz(h=8,w=8)
par(mfrow=c(2,2))
for(catnum in 1:length(cats)) {
  catCH = subset(NemegtADDITIONAL_ch,subset=which(row.names(NemegtADDITIONAL_ch)==cats[catnum]))
  plotcovariate(NemegtMask, covariate="GC", contour = FALSE, col = terrain.colors(16), xlim=c(730000,740000),zlim=zlim,asp=1,
                main=cats[catnum])
  plot(traps(catCH),add=TRUE)
  plot(catCH,tracks=TRUE,add=TRUE)
  points(traps(catCH)$x,traps(catCH)$y,pch=19,cex=1.5,col="white")
  ncap = as.character(catCH[1,1,])
  text(traps(catCH)$x,traps(catCH)$y,labels=ncap,cex=0.75)
}


# Remove traps that 4 cats above were caught on, and all further east
which(is.element(row.names(NemegtADDITIONAL_ch),cats))
catstraps = which(apply(NemegtADDITIONAL_ch[cats,1,],2,sum)>0)
mineast = min(traps(NemegtADDITIONAL_ch)[catstraps,]$x)
westraps = which(traps(NemegtADDITIONAL_ch)$x<mineast)
eastraps = which(traps(NemegtADDITIONAL_ch)$x>=mineast)
NemegtADDwest_ch = subset(NemegtADDITIONAL_ch,traps=westraps)
NemegtADDeast_ch = subset(NemegtADDITIONAL_ch,traps=eastraps)
xlim.west = range(c(traps(NemegtADDwest_ch)$x,min(NemegtMask$x)))
xlim.east = range(traps(NemegtADDeast_ch)$x)
# Look at the two subsets
quartz(h=8,w=8)
par(mfrow=c(2,1))
plotcovariate(NemegtMask, covariate="GC", contour=FALSE, col=terrain.colors(16), 
              zlim=zlim,asp=1,xlim=xlim.west)
plot(traps(NemegtADDwest_ch),add=TRUE)
plot(NemegtADDwest_ch,tracks=TRUE,add=TRUE)
plotcovariate(NemegtMask, covariate="GC", contour=FALSE, col=terrain.colors(16), 
              zlim=zlim,asp=1,xlim=xlim.east)
plot(traps(NemegtADDeast_ch),add=TRUE)
plot(NemegtADDeast_ch,tracks=TRUE,add=TRUE)

# Fit same model to east and west bits:
NemegtADDwest<-secr.fit(NemegtADDwest_ch, detectfn="HHN", mask=NemegtMask,
                           model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                           details = list(userdist = userdfn1),
                           start = list(noneuc = 1))
#Fuck - that STILL gives negative signs!!

# Can't determine staring values for data below - leave for now
#NemegtADDeast<-secr.fit(NemegtADDeast_ch, detectfn="HHN", mask=NemegtMask,
                        model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                        details = list(userdist = userdfn1),
                        start = list(noneuc = 1))

# Try cutting at x=680000 because recatpures around x=700000 have same features as those farther east:
mineast2 = 680000
westraps2 = which(traps(NemegtADDITIONAL_ch)$x<mineast2)
eastraps2 = which(traps(NemegtADDITIONAL_ch)$x>=mineast2)
NemegtADDwest2_ch = subset(NemegtADDITIONAL_ch,traps=westraps2)
NemegtADDeast2_ch = subset(NemegtADDITIONAL_ch,traps=eastraps2)
xlim.west2 = range(c(traps(NemegtADDwest2_ch)$x,min(NemegtMask$x)))
xlim.east2 = range(traps(NemegtADDeast2_ch)$x)
# Look at the two subsets
quartz(h=8,w=8)
par(mfrow=c(2,1))
plotcovariate(NemegtMask, covariate="GC", contour=FALSE, col=terrain.colors(16), 
              zlim=zlim,asp=1,xlim=xlim.west2)
plot(traps(NemegtADDwest2_ch),add=TRUE)
plot(NemegtADDwest2_ch,tracks=TRUE,add=TRUE)
plotcovariate(NemegtMask, covariate="GC", contour=FALSE, col=terrain.colors(16), 
              zlim=zlim,asp=1,xlim=xlim.east2)
plot(traps(NemegtADDeast2_ch),add=TRUE)
plot(NemegtADDeast2_ch,tracks=TRUE,add=TRUE)

# Fit same model to east and west bits:
NemegtADDwest2<-secr.fit(NemegtADDwest2_ch, detectfn="HHN", mask=NemegtMask,
                        model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                        details = list(userdist = userdfn1),
                        start = list(noneuc = 1))
NemegtADDeast2<-secr.fit(NemegtADDeast2_ch, detectfn="HHN", mask=NemegtMask,
                         model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                         details = list(userdist = userdfn1),
                         start = list(noneuc = 1))
# WTF??! Both of these have NEGATIVE coefficient for D..stdGC and noneuc.stdGC?!!


# Try with simpler least-cost distance functions
# ----------------------------------------------
mytransfun = function(x) exp(mean(x))
myLCdistfun <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') # required by secr
  require(gdistance) # to load transition and geoCorrection functions
  Sraster <- raster(mask, 'noneuc') # Make raster from mesh with covariate 'noneuc'
  # Calculate all the conductances, using mytransfun
  trans <- transition(Sraster, transitionFunction = mytransfun,  directions = 16)
  # Adjust for distance difference between square and diagonal neighbours
  trans <- geoCorrection(trans)
  # calculate the least-cost distance
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}
# Fit same model to east and west bits, using simpler least-cost distance function:
NemegtADDwest2a<-secr.fit(NemegtADDwest2_ch, detectfn="HHN", mask=NemegtMask,
                         model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                         details = list(userdist = myLCdistfun),
                         start = list(noneuc = 1),
                         link = list(noneuc="identity"))
NemegtADDeast2a<-secr.fit(NemegtADDeast2_ch, detectfn="HHN", mask=NemegtMask,
                         model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                         details = list(userdist = myLCdistfun),
                         start = list(noneuc = 1),
                         link = list(noneuc="identity"))
DhatADDwest2a = predictDsurface(NemegtADDwest2a,parameter="D")
NEhatADDwest2a = predictDsurface(NemegtADDwest2a,parameter="noneuc")
DhatADDeast2a = predictDsurface(NemegtADDeast2a,parameter="D")
NEhatADDeast2a = predictDsurface(NemegtADDeast2a,parameter="noneuc")

quartz(h=5,w=10)
par(mfrow=c(2,2))
plotcovariate(DhatADDwest2a,covariate="D.0",asp=1,main="Density",contour=FALSE)
plotcovariate(NEhatADDwest2a,covariate="noneuc.0",asp=1,main="noneuc",contour=FALSE)
plotcovariate(DhatADDeast2a,covariate="D.0",asp=1,main="Density",contour=FALSE)
plotcovariate(NEhatADDeast2a,covariate="noneuc.0",asp=1,main="noneuc",contour=FALSE)

# Try combined east and west, using simpler least-cost distance function:
NemegtADDa<-secr.fit(NemegtADDITIONAL_ch, detectfn="HHN", mask=NemegtMask,
                          model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                          details = list(userdist = myLCdistfun),
                          start = list(noneuc = 1),
                     link = list(noneuc="identity"))
DhatADDa = predictDsurface(NemegtADDa,parameter="D")
NEhatADDa = predictDsurface(NemegtADDa,parameter="noneuc")
quartz(h=8,w=8)
par(mfrow=c(2,1))
plotcovariate(DhatADDa,covariate="D.0",asp=1,main="Density",contour=FALSE)
plotcovariate(NEhatADDa,covariate="noneuc.0",asp=1,main="noneuc",contour=FALSE)
# Hoorayyyyy - that seems to work !!!
# Try without the artificial cameras:
NemegtALLa<-secr.fit(NemegtALL_ch, detectfn="HHN", mask=NemegtMask,
                     model=list(D~stdGC, lambda0~Water, sigma~1, noneuc ~ stdGC -1), 
                     details = list(userdist = myLCdistfun),
                     start = list(noneuc = 1),
                     link = list(noneuc="identity"))
NemegtDhatALLa = predictDsurface(NemegtALLa,parameter="D")
NemegtNEhatALLa = predictDsurface(NemegtALLa,parameter="noneuc")
#quartz(h=8,w=8)
pdf(file="./Analysis4paper/Plots/NemegtLCdistFit.pdf")
par(mfrow=c(2,1))
plotcovariate(NemegtDhatALLa,covariate="D.0",asp=1,main="Nemegt Density",contour=FALSE)
plotcovariate(NemegtNEhatALLa,covariate="noneuc.0",asp=1,main="Nemegt noneuc",contour=FALSE)
dev.off()

# Now try the other two regions using simpler least-cost distance function (fingers crossed!!):
Noyon.NUstdGCa<-secr.fit(Noyon_ch, detectfn="HHN", mask=NoyonMask,
                        model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                        details = list(userdist = myLCdistfun),
                        start = list(noneuc = 1),
                        link = list(noneuc="identity"))
Tost.NUstdGCa<-secr.fit(Tost_ch, detectfn="HHN", mask=TostMask,
                       model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                       details = list(userdist = myLCdistfun),
                       start = list(noneuc = 1),
                       link = list(noneuc="identity"))
# Did not converge, so try again starting where left off
Tost.NUstdGCa<-secr.fit(Tost_ch, detectfn="HHN", mask=TostMask,
                        model=list(D~stdGC, lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                        details = list(userdist = myLCdistfun),
                        start = Tost.NUstdGCa,
                        link = list(noneuc="identity"))

NoyonDhatALLa = predictDsurface(Noyon.NUstdGCa,parameter="D")
NoyonNEhatALLa = predictDsurface(Noyon.NUstdGCa,parameter="noneuc")
#quartz(h=8,w=8)
pdf(file="./Analysis4paper/Plots/NoyonLCdistFit.pdf")
par(mfrow=c(2,1))
plotcovariate(NoyonDhatALLa,covariate="D.0",asp=1,main="Noyon Density",contour=FALSE)
plotcovariate(NoyonNEhatALLa,covariate="noneuc.0",asp=1,main="Noyon noneuc",contour=FALSE)
dev.off()

TostDhatALLa = predictDsurface(Tost.NUstdGCa,parameter="D")
TostNEhatALLa = predictDsurface(Tost.NUstdGCa,parameter="noneuc")
#quartz(h=8,w=8)
pdf(file="./Analysis4paper/Plots/TostLCdistFit.pdf")
par(mfrow=c(2,1))
plotcovariate(TostDhatALLa,covariate="D.0",asp=1,main="Tost Density",contour=FALSE)
plotcovariate(TostNEhatALLa,covariate="noneuc.0",asp=1,main="Tost noneuc",contour=FALSE)
dev.off()
