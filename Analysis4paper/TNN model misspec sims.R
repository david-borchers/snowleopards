library(secr)
library(fields)
library(maptools)
source("scrplotting.r")
library(scrmlebook)
library(parallel)
library(tcltk2) # for progress bar function


if(.Platform$OS.type=="windows") { # this to make the command quartz( ) work on windows machines
  quartz<-function() windows()
}


#-------------------------------------------------------------------------------
# Some functions needed to simulate detection histories when there are 
# trap-specific covariates
#---------------------------------------------------------------------
# Calculate expected encouner rate or detection probability from list of parameters with a0 and sigma
er_from_a0sig = function(parlist,dists,usage=rep(1,length(parlist)),calcprob=FALSE) {
  if(!inherits(parlist,"list")) stop("'parlist' must be a list.")
  npars = length(parlist)
  if(is.vector(dists)) dists = matrix(dists,ncol=1)
  if(!is.matrix(dists)) stop("dists must be a vector or a matrix.")
  if(dim(dists)[1]!=npars) stop("dists must be of same length as npars, or if matrix, have npars rows.")
  ndists = dim(dists)[2]
  er = matrix(rep(NA,npars*ndists),nrow=npars)
  for(i in 1:npars) {
    a0 = parlist[[i]]["a0","estimate"] * 10^4 # conversion from hectaers to sq m
    sigma = parlist[[i]]["sigma","estimate"]
    lambda0 = a0/(2*pi*sigma^2)
    er[i,] = lambda0*exp(-dists[i,]^2/(2*sigma^2)) * usage[i]
  }
  if(calcprob) {
    p = 1 - exp(-er)
    return(p)
  } else {
    return(er)
  }
}



distances <- function (X, Y) {
  ## X and Y are 2-column matrices of coordinates
  onerow <- function (xy) {
    d <- function(xy2) {
      sqrt(sum((xy2 - xy)^2))
    }
    apply(Y, 1, d)
  }
  t(apply(X, 1, onerow))
}

R.N = function(fit) {
  Ns = region.N(fit)
  nNs = length(Ns)
  Nhat = rep(NA,nNs)
  for(i in 1:nNs) Nhat[i] = Ns[[i]]["R.N","estimate"]
  return(Nhat)
}
#-------------------------------------------------------------------------------
# A function to plot a shaded confidence interval
# -----------------------------------------------
#' @title Draws a shaded confidence region for a line
#'
#' @description
#'  Adds a shaded polygon covering the region between lower and upper confidence
#'  bounds, to an existing plot.
#'  
#' @param xyci a data frame or matrix with first element being x-coordinates, 
#' second being the corresponding lower (or upper) confidence bound values, 
#' third being the corresponding upper (or lower) confidence bound values. 
#' @param col an argument of \code{polygon} specifying shading colour.
#' @param border an argument of \code{polygon} specifying border colour
#' @param ... other arguments to \code{polygon}.
#' 
#' @export addCIpoly
#' 
shadedCI = function(xyci,col=adjustcolor("gray", alpha.f=0.4),border=NA,...) {
  if(!inherits(xyci,"data.frame") & !inherits(xyci,"matrix")) 
    stop("xyci must be a matrix or data frame.")
  if(ncol(xyci)!=3) stop("xyci must have 3 columns.")
  ci.x <- c(xyci[,1], rev(xyci[,1]))
  ci.y <- c(xyci[,2], rev(xyci[,3]))
  polygon(x=ci.x, y=ci.y,col=col,border=border,...)
}
#-------------------------------------------------------------------------------

ncores = detectCores()

#----------------------- Load the RData objects made by Make_TNN_RData.r --------------------------
load("./Analysis4paper/TNN_boundaries.RData") # Tostboundary,Noyonboundary,Nemegtboundary
load("./Analysis4paper/TNN_masks.RData") # TostMask,NoyonMask,NemegtMask
load("./Analysis4paper/TNN_caphists.RData") # Tost_ch,Noyon_ch,Nemegt_chTNN_ch

#----------------------- Load fitted models --------- --------------------------
load("./Analysis4paper/TNNnoNEfits.RData")

# best model among those with Euclidian distance
fit = TNNfit.DGrid.a0Waterxsess_topo.sig1
# ... OR
# 2nd best:
fit2 = TNNfit.DGrid.a0Waterxsess.sig1

# best model among those with flat density
fit0 = TNNfit.D1.LamSess.sig1
fitsite = TNNfit.Dsess.Lamsess.sigsess


# Extract data
Tost = dat$Tost
Nemegt = dat$Nemegt
Noyon = dat$Noyon

# predict density 
Dhats.mask = predictDsurface(fit,se.D=TRUE,cl.D=TRUE)
Dhats = list(Tost=NA, Noyon=NA, Nemegt=NA)
mean.Dhat = rep(NA,3)
for(i in 1:3) {
#  covariates(Dhats.mask[[i]])$D.0 = covariates(Dhats.mask[[i]])$D.0 #*10^4 # scale to number per 100 km^2
  Dhats[[i]]=covariates(Dhats.mask[[i]])$D.0
  mean.Dhat[i] = mean(Dhats[[i]])
}

Dhat = Dhat1 = vector(length=3,mode="list")
Dmean = D0 = Dsite = D = D.lcl = D.ucl = D1 = D1.lcl = D1.ucl = stdGC = data.frame(Tost=rep(NA,100),Noyon=rep(NA,100),Nemgt=rep(NA,100))
for(i in 1:3) {
  stdGCrange = range(covariates(fit$mask[[i]])$stdGC)
  stdGC[,i] = seq(stdGCrange[1],stdGCrange[2],length=100)
  newdata = data.frame(
    stdGC=stdGC[,i],
    Topo = rep(covariates(traps(fit$capthist[[1]]))$Topo[1],100),  # [1] is arbitrary value
    Water = rep(covariates(traps(fit$capthist[[1]]))$Water[1],100),  # [1] is arbitrary value
    session = rep(factor(as.character(i),levels=c("1","2","3")),100)
  )
  newdata1 = data.frame(
    stdGC=stdGC[,i],
    Topo = rep(factor("Steppe",levels=c("Canyon","Saddle","Steppe")),100),  # arbitrary value
    Water = rep(factor("Yes",levels=c("Yes","No")),100),  # arbitrary value
    session = rep(factor(as.character(i),levels=c("1","2","3")),100)
  )
  Dhat[[i]] = predict(fit,newdata=newdata)
  Dhat1[[i]] = predict(fit,newdata=newdata1)
}
# See if get same with different Water and Topo:
Dhat0 =  predict(fit0)
Dhatsite =  predict(fitsite)
pclow = rep(NA,3)
for(j in 1:3) {
  pclow[j] = 100*(Dhat0[[j]][1,2]-mean.Dhat[j])/Dhat0[[j]][1,2]
  for(i in 1:100) {
    D[i,j] = Dhat[[j]][[i]][1,2]*10^4
    D.lcl[i,j] = Dhat[[j]][[i]][1,4]*10^4
    D.ucl[i,j] = Dhat[[j]][[i]][1,5]*10^4
    D1[i,j] = Dhat1[[j]][[i]][1,2]*10^4
    D1.lcl[i,j] = Dhat1[[j]][[i]][1,4]*10^4
    D1.ucl[i,j] = Dhat1[[j]][[i]][1,5]*10^4
    D0[,j] = rep(Dhat0[[j]][1,2]*10^4,100)
    Dsite[,j] = rep(Dhatsite[[j]][1,2]*10^4,100)
    Dmean[,j] = rep(mean.Dhat[j]*10^4,100)
  }
}

# Estimate total density over region using best and constant-D models
# ------------------------------------------------------------------
Nhat = region.N(fit)
N2hat = region.N(fit2)
N0hat = region.N(fit0)
Nsitehat = region.N(fitsite)
Dest = D2est = D0est = Dsiteest = data.frame(Dhat=rep(NA,3),lcl=rep(NA,3),ucl=rep(NA,3))
for(i in 1:3) {
  Dest[i,] = Nhat[[i]]["E.N",c(1,3,4)]/(attr(fit$mask[[i]],"area")*dim(fit$mask[[i]])[1])*10^4
  D2est[i,] = N2hat[[i]]["E.N",c(1,3,4)]/(attr(fit2$mask[[i]],"area")*dim(fit2$mask[[i]])[1])*10^4
  Dsiteest[i,] = Nsitehat[[i]]["E.N",c(1,3,4)]/(attr(fitsite$mask[[i]],"area")*dim(fitsite$mask[[i]])[1])*10^4
  D0est[i,] = N0hat[[i]]["E.N",c(1,3,4)]/(attr(fit0$mask[[i]],"area")*dim(fit0$mask[[i]])[1])*10^4
}
Dest
D2est

100*(Dest$Dhat-D2est$Dhat)/Dest$Dhat

Dsiteest
D0est


pcDsite = 100*(Dest$Dhat-Dsiteest$Dhat)/Dsiteest$Dhat
pcD0 = 100*(Dest$Dhat-D0est$Dhat)/D0est$Dhat
pcDsite
pcD0

DCIwidth = Dest$ucl - Dest$lcl
D2CIwidth = D2est$ucl - D2est$lcl
100*(DCIwidth-D2CIwidth)/DCIwidth

DsiteCIwidth = Dsiteest$ucl - Dsiteest$lcl
D0CIwidth = D0est$ucl - D0est$lcl
pcCIwidthsite = 100*(DCIwidth-DsiteCIwidth)/DsiteCIwidth
pcCIwidth0 = 100*(DCIwidth-D0CIwidth)/D0CIwidth
pcCIwidthsite
pcCIwidth0

# differences as % of CI width
100*(Dest$Dhat-Dsiteest$Dhat)/DsiteCIwidth
100*(Dest$Dhat-D0est$Dhat)/D0CIwidth


max(D.ucl)
max(D1.ucl)
max(covariates(Dhats.mask[[1]])$ucl.0*10^4)
max(covariates(Dhats.mask[[2]])$ucl.0*10^4)
max(covariates(Dhats.mask[[3]])$ucl.0*10^4)

pclow # % by which constant density model density is below mean of best model density


# Plot density against stdGC
pdf(file="./Analysis4paper/Plots/DensityVsSdGC.pdf",h=3,w=11)
par(mfrow=c(1,3))
par(mar=c(4,5,2,2))
for(i in 1:3) {
  plot(stdGC[,i],D[,i],type="l",ylim=c(min(D.lcl[,i]),2*max(D[,i])),main=names(D)[i],xlab="stdGC",ylab=expression(hat(D)))
  shadedCI(data.frame(stdGC[,i],D.lcl[,i],D.ucl[,i]))
  lines(stdGC[,i],Dmean[,i],col="darkgray")
  lines(stdGC[,i],Dsite[,i],lty=2,col="darkgray")
  trapcov = covariates(traps(fit$capthist[[i]]))$stdGC
  points(trapcov,rep(min(D.lcl[,i]),length(trapcov)),pch="|")
}
dev.off()
# again, to check with differen Topo and Water values:
quartz(h=3,w=11)
par(mfrow=c(1,3))
par(mar=c(4,5,2,2))
for(i in 1:3) {
  plot(stdGC[,i],D1[,i],type="l",ylim=c(min(D1.lcl[,i]),2*max(D1[,i])),main=names(D1)[i],xlab="stdGC",ylab=expression(hat(D1)))
  shadedCI(data.frame(stdGC[,i],D1.lcl[,i],D1.ucl[,i]))
  lines(stdGC[,i],Dmean[,i],col="darkgray")
  lines(stdGC[,i],Dsite[,i],lty=2,col="darkgray")
  trapcov = covariates(traps(fit$capthist[[i]]))$stdGC
  points(trapcov,rep(min(D1.lcl[,i]),length(trapcov)),pch="|")
}


quartz(h=8,w=3)
par(mfrow=c(3,1))
hist(covariates(fit$mask[[1]])$GC)
hist(covariates(fit$mask[[2]])$GC)
hist(covariates(fit$mask[[3]])$GC)

hist(covariates(fit$mask[[1]])$GRIDCODE)
hist(covariates(fit$mask[[2]])$GRIDCODE)
hist(covariates(fit$mask[[3]])$GRIDCODE)

# Plot stdGC in space
xlim = range(fit$mask[[1]]$x,fit$mask[[2]]$x,fit$mask[[3]]$x)
ylim = range(fit$mask[[1]]$y,fit$mask[[2]]$y,fit$mask[[3]]$y)
quartz(h=4,w=9)
layout(mat=matrix(1:2,nrow=1),widths=c(.9,.1))
splotcovariate(fit$mask[[1]],"stdGC",what="image",zlim=stdGCrange,xlim=xlim,ylim=ylim,col=parula(40),
               main="stdGC",cex.main=0.75)
plot(traps(fit$capthist[[1]]),add=TRUE)
splotcovariate(fit$mask[[2]],"stdGC",what="image",zlim=stdGCrange,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[2]]),add=TRUE)
splotcovariate(fit$mask[[3]],"stdGC",what="image",zlim=stdGCrange,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[3]]),add=TRUE)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
splotcovariate(fit$mask[[1]],"stdGC",what="scale",col=parula(40))


# Plot saddle/canyon/steppe at traps
#quartz(h=4,w=9)
pdf(file="./Analysis4paper/Plots/trapsTopo.pdf",h=4,w=9)
layout(mat=matrix(1:2,nrow=1),widths=c(.9,.1))
splotcovariate(fit$mask[[1]],"stdGC",what="image",zlim=stdGCrange,xlim=xlim,ylim=ylim,col=parula(40),
               main="stdGC",cex.main=0.75)
i=1
cams = traps(fit$capthist)[[i]]
tcol=rep(0,dim(cams)[1])
tcol[covariates(cams)$Topo=="Saddle"] = "red"
tcol[covariates(cams)$Topo=="Canyon"] = "black"
tcol[covariates(cams)$Topo=="Steppe"] = "green"
points(cams$x,cams$y,pch=19,col=tcol)
splotcovariate(fit$mask[[2]],"stdGC",what="image",zlim=stdGCrange,add=TRUE,col=parula(40))
i=2
cams = traps(fit$capthist)[[i]]
tcol=rep(0,dim(cams)[1])
tcol[covariates(cams)$Topo=="Saddle"] = "red"
tcol[covariates(cams)$Topo=="Canyon"] = "black"
tcol[covariates(cams)$Topo=="Steppe"] = "green"
points(cams$x,cams$y,pch=19,col=tcol)
splotcovariate(fit$mask[[3]],"stdGC",what="image",zlim=stdGCrange,add=TRUE,col=parula(40))
legend("topleft",bty="n",pch=rep(19,3),col=c("red","black","green"),legend=c("Saddle","Canyon","Steppe"),cex=0.75)
i=3
cams = traps(fit$capthist)[[i]]
tcol=rep(0,dim(cams)[1])
tcol[covariates(cams)$Topo=="Saddle"] = "red"
tcol[covariates(cams)$Topo=="Canyon"] = "black"
tcol[covariates(cams)$Topo=="Steppe"] = "green"
points(cams$x,cams$y,pch=19,col=tcol)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
splotcovariate(fit$mask[[1]],"stdGC",what="scale",col=parula(40))
dev.off()


# Plot Water at traps
#quartz(h=4,w=9)
pdf(file="./Analysis4paper/Plots/trapsWater.pdf",h=4,w=9)
layout(mat=matrix(1:2,nrow=1),widths=c(.9,.1))
splotcovariate(fit$mask[[1]],"stdGC",what="image",zlim=stdGCrange,xlim=xlim,ylim=ylim,col=parula(40),
               main="stdGC",cex.main=0.75)
i=1
cams = traps(fit$capthist)[[i]]
wcol=rep(0,dim(cams)[1])
wcol[covariates(cams)$Water=="Yes"] = "red"
wcol[covariates(cams)$Water=="No"] = "black"
points(cams$x,cams$y,pch=19,col=wcol)
splotcovariate(fit$mask[[2]],"stdGC",what="image",zlim=stdGCrange,add=TRUE,col=parula(40))
i=2
cams = traps(fit$capthist)[[i]]
wcol=rep(0,dim(cams)[1])
wcol[covariates(cams)$Water=="Yes"] = "red"
wcol[covariates(cams)$Water=="No"] = "black"
points(cams$x,cams$y,pch=19,col=wcol)
splotcovariate(fit$mask[[3]],"stdGC",what="image",zlim=stdGCrange,add=TRUE,col=parula(40))
legend("topleft",bty="n",pch=rep(19,2),col=c("red","black"),legend=c("Yes","No"),cex=0.75)
i=3
cams = traps(fit$capthist)[[i]]
wcol=rep(0,dim(cams)[1])
wcol[covariates(cams)$Water=="Yes"] = "red"
wcol[covariates(cams)$Water=="No"] = "black"
points(cams$x,cams$y,pch=19,col=wcol)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
splotcovariate(fit$mask[[1]],"stdGC",what="scale",col=parula(40))
dev.off()


# Look at the "design" for Topo and Water at traps:
# ------------------------------------------------
tab = vector(length=3,mode="list")
for(i in 1:3) tab[[i]] = table(covariates(traps(fit$capthist[[i]]))[,c("Topo","Water")])
tab



# Plot Density in space
# ------------------------------------------------
D = Dhats.mask
Dmin = DCVmin = Dlclmin = Duclmin = 99999999
Dmax = DCVmax = Dlclmax = Duclmax = -99999999
for(i in 1:3) {
  covariates(D[[i]])$Dhat = covariates(D[[i]])$D.0 *10^4
  covariates(D[[i]])$Dlcl = covariates(D[[i]])$lcl.0 *10^4
  covariates(D[[i]])$Ducl = covariates(D[[i]])$ucl.0 *10^4
  covariates(D[[i]])$DhatCV = 100*covariates(D[[i]])$SE.0/covariates(D[[i]])$D.0
  Dmin = min(Dmin,covariates(D[[i]])$Dhat)
  Dmax = max(Dmax,covariates(D[[i]])$Dhat)
  DCVmin = min(DCVmin,covariates(D[[i]])$DhatCV)
  DCVmax = max(DCVmax,covariates(D[[i]])$DhatCV)
  Dlclmin = min(Dlclmin,covariates(D[[i]])$Dlcl)
  Dlclmax = max(Dlclmax,covariates(D[[i]])$Dlcl)
  Duclmin = min(Duclmin,covariates(D[[i]])$Ducl)
  Duclmax = max(Duclmax,covariates(D[[i]])$Ducl)
}

# Plot density maps
# -----------------
xlim = range(D[[1]]$x,D[[2]]$x,D[[3]]$x)
ylim = range(D[[1]]$y,D[[2]]$y,D[[3]]$y)
Dlim=c(Dmin,Dmax)
#quartz(h=4,w=9)
pdf(file="./Analysis4paper/Plots/Dhat.pdf",h=4,w=9)
layout(mat=matrix(1:2,nrow=1),widths=c(.9,.1))
splotcovariate(D[[1]],"Dhat",what="image",zlim=Dlim,xlim=xlim,ylim=ylim,col=parula(40),
               main="Number per 100 sq km",cex.main=0.75)
plot(traps(fit$capthist[[1]]),add=TRUE)
splotcovariate(D[[2]],"Dhat",what="image",zlim=Dlim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[2]]),add=TRUE)
splotcovariate(D[[3]],"Dhat",what="image",zlim=Dlim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[3]]),add=TRUE)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
splotcovariate(D[[1]],"Dhat",what="scale",col=parula(40))
dev.off()

# Plot density CV maps
# -----------------
xlim = range(D[[1]]$x,D[[2]]$x,D[[3]]$x)
ylim = range(D[[1]]$y,D[[2]]$y,D[[3]]$y)
DCVlim=c(DCVmin,DCVmax)
#quartz(h=4,w=9)
pdf(file="./Analysis4paper/Plots/DpcCV.pdf",h=4,w=9)
layout(mat=matrix(1:2,nrow=1),widths=c(.9,.1))
splotcovariate(D[[1]],"DhatCV",what="image",zlim=DCVlim,xlim=xlim,ylim=ylim,col=parula(40),
               main="% Coefficient of variation",cex.main=0.75)
plot(traps(fit$capthist[[1]]),add=TRUE)
splotcovariate(D[[2]],"DhatCV",what="image",zlim=DCVlim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[2]]),add=TRUE)
splotcovariate(D[[3]],"DhatCV",what="image",zlim=DCVlim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[3]]),add=TRUE)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
splotcovariate(D[[1]],"DhatCV",what="scale",col=parula(40))
dev.off()

# Plot density lcl and ucl maps on same scale
# -------------------------------------------
xlim = range(D[[1]]$x,D[[2]]$x,D[[3]]$x)
ylim = range(D[[1]]$y,D[[2]]$y,D[[3]]$y)
DCIlim=range(Dlclmin,Dlclmax,Duclmin,Duclmax)
#quartz(h=4,w=9)
pdf(file="./Analysis4paper/Plots/DlowerCI.pdf",h=4,w=9)
layout(mat=matrix(1:2,nrow=1),widths=c(.9,.1))
splotcovariate(D[[1]],"Dlcl",what="image",zlim=DCIlim,xlim=xlim,ylim=ylim,col=parula(40),
               main="Lower CI for D",cex.main=0.75)
plot(traps(fit$capthist[[1]]),add=TRUE)
splotcovariate(D[[2]],"Dlcl",what="image",zlim=DCIlim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[2]]),add=TRUE)
splotcovariate(D[[3]],"Dlcl",what="image",zlim=DCIlim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[3]]),add=TRUE)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
splotcovariate(D[[1]],"Dlcl",what="scale",zlim=DCIlim,col=parula(40))
dev.off()

# Plot density ucl maps
# -----------------
xlim = range(D[[1]]$x,D[[2]]$x,D[[3]]$x)
ylim = range(D[[1]]$y,D[[2]]$y,D[[3]]$y)
#quartz(h=4,w=9)
pdf(file="./Analysis4paper/Plots/DupperCI.pdf",h=4,w=9)
layout(mat=matrix(1:2,nrow=1),widths=c(.9,.1))
splotcovariate(D[[1]],"Ducl",what="image",zlim=DCIlim,xlim=xlim,ylim=ylim,col=parula(40),
               main="Upper CI for D",cex.main=0.75)
plot(traps(fit$capthist[[1]]),add=TRUE)
splotcovariate(D[[2]],"Ducl",what="image",zlim=DCIlim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[2]]),add=TRUE)
splotcovariate(D[[3]],"Ducl",what="image",zlim=DCIlim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[3]]),add=TRUE)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
splotcovariate(D[[1]],"Ducl",what="scale",zlim=DCIlim,col=parula(40))
dev.off()

# Calculate and plot coverage probability
# ---------------------------------------
p = vector(length=3,mode="list")
for(i in 1:3) {
  Dsum = sum(covariates(D[[i]])$Dhat)
  cams = traps(fit$capthist[[i]])
  ncams = dim(cams)[1]
  Topo = covariates(cams)$Topo
  Water = covariates(cams)$Water
  # create the session factor for this region (needed for predict() to work, 
  # as seesion is in original fitted model)
  sess = factor(as.character(i),levels=as.character(1:3))
  # predict using a data frame with a row per trap, with associated covariates
  # stdGC has to be in because it is in the fitted model, but we don't use it
  # as we extract only the detection function parameters from predpar
  predpar = predict(fit,newdata=data.frame(stdGC=rep(0,ncams),Topo=Topo,Water=Water,session=rep(sess,ncams)))
  # Calulate distances from all traps to all mask points
  dists = distances(cams,fit$mask[[i]]) # distances from traps to each activity centre
  # Get expected number of encouners at each trap using detection parameters 
  # (a0 and sigma in this case) at each trap. These parameters have been 
  # obtained by predict() above, using the appropriate covariates at each trap.
  en = er_from_a0sig(predpar,dists,usage(cams)) # expected number of encouters
  entot = apply(en,2,sum) # expected number of encounters at each mask point, aggregating en across traps
  p[[i]] = 1 - exp(-entot) # p at each mask point
  covariates(D[[i]])$er = entot
  covariates(D[[i]])$p = p[[i]] # p(see|location)
  covariates(D[[i]])$cover = p[[i]]/dim(D[[i]])[1] # space coverage prob
  covariates(D[[i]])$cover.pop = p[[i]]*covariates(D[[i]])$Dhat/Dsum # pop coverage prob
}

# Plot expected encounter rate:
xlim = range(D[[1]]$x,D[[2]]$x,D[[3]]$x)
ylim = range(D[[1]]$y,D[[2]]$y,D[[3]]$y)
erlim = range(covariates(D[[1]])$er,covariates(D[[2]])$er,covariates(D[[3]])$er)
#quartz(h=4,w=9)
pdf(file="./Analysis4paper/Plots/ExpectedEncounters.pdf",h=4,w=9)
layout(mat=matrix(1:2,nrow=1),widths=c(.9,.1))
splotcovariate(D[[1]],"er",what="image",zlim=erlim,xlim=xlim,ylim=ylim,col=parula(40),
               main="Expected number of encounters",cex.main=0.75)
plot(traps(fit$capthist[[1]]),add=TRUE)
splotcovariate(D[[2]],"er",what="image",zlim=erlim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[2]]),add=TRUE)
splotcovariate(D[[3]],"er",what="image",zlim=erlim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[3]]),add=TRUE)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
splotcovariate(D[[1]],"er",what="scale",zlim=erlim,col=parula(40))
dev.off()


# Plot space coverage probability by points, p:
xlim = range(D[[1]]$x,D[[2]]$x,D[[3]]$x)
ylim = range(D[[1]]$y,D[[2]]$y,D[[3]]$y)
plim = c(0,1)
#quartz(h=4,w=9)
pdf(file="./Analysis4paper/Plots/CoverageProb.pdf",h=4,w=9)
layout(mat=matrix(1:2,nrow=1),widths=c(.9,.1))
splotcovariate(D[[1]],"p",what="image",zlim=plim,xlim=xlim,ylim=ylim,col=parula(40),
               main="Coverage probability",cex.main=0.75)
plot(traps(fit$capthist[[1]]),add=TRUE)
splotcovariate(D[[2]],"p",what="image",zlim=plim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[2]]),add=TRUE)
splotcovariate(D[[3]],"p",what="image",zlim=plim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[3]]),add=TRUE)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
splotcovariate(D[[1]],"p",what="scale",zlim=plim,col=parula(40))
dev.off()



# Plot population coverage probability:
xlim = range(D[[1]]$x,D[[2]]$x,D[[3]]$x)
ylim = range(D[[1]]$y,D[[2]]$y,D[[3]]$y)
plim = range(covariates(D[[1]])$cover.pop,covariates(D[[2]])$cover.pop,covariates(D[[3]])$cover.pop)
#quartz(h=4,w=9)
pdf(file="./Analysis4paper/Plots/PopCoverageProb.pdf",h=4,w=9)
layout(mat=matrix(1:2,nrow=1),widths=c(.9,.1))
splotcovariate(D[[1]],"cover.pop",what="image",zlim=plim,xlim=xlim,ylim=ylim,col=parula(40),
               main="Coverage probability",cex.main=0.75)
plot(traps(fit$capthist[[1]]),add=TRUE)
splotcovariate(D[[2]],"cover.pop",what="image",zlim=plim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[2]]),add=TRUE)
splotcovariate(D[[3]],"cover.pop",what="image",zlim=plim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[3]]),add=TRUE)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
splotcovariate(D[[1]],"cover.pop",what="scale",zlim=plim,col=parula(40))
dev.off()



# Plot spatial coverage probability:
xlim = range(D[[1]]$x,D[[2]]$x,D[[3]]$x)
ylim = range(D[[1]]$y,D[[2]]$y,D[[3]]$y)
plim = range(covariates(D[[1]])$cover,covariates(D[[2]])$cover,covariates(D[[3]])$cover)
#quartz(h=4,w=9)
pdf(file="./Analysis4paper/Plots/SpaceCoverageProb.pdf",h=4,w=9)
layout(mat=matrix(1:2,nrow=1),widths=c(.9,.1))
splotcovariate(D[[1]],"cover",what="image",zlim=plim,xlim=xlim,ylim=ylim,col=parula(40),
               main="Coverage probability",cex.main=0.75)
plot(traps(fit$capthist[[1]]),add=TRUE)
splotcovariate(D[[2]],"cover",what="image",zlim=plim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[2]]),add=TRUE)
splotcovariate(D[[3]],"cover",what="image",zlim=plim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[3]]),add=TRUE)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
splotcovariate(D[[1]],"cover",what="scale",zlim=plim,col=parula(40))
dev.off()


# Look at total abundance estimates, and CIs and CVs for best and const-site-D model:
# -------------------------------------------------------------------------------
Nbest = region.N(fit)
Nsite = region.N(fitsite)
Nhats = array(rep(NA,5*3*2),dim=c(3,6,2),
              dimnames=list(region=c("Tost","Noyon","Nemegt"),
                            est=c("point","SE","lcl","ucl","n","CV"),
                            model=c("best","const-site-D")))
for(i in 1:3) {
  Nhats[i,,1] = c(as.numeric(Nbest[[i]]["E.N",]),Nbest[[i]]["E.N",2]/Nbest[[i]]["E.N",1])
  Nhats[i,,2] = c(as.numeric(Nsite[[i]]["E.N",]),Nsite[[i]]["E.N",2]/Nsite[[i]]["E.N",1])
}
Nhats
pcCVdiff = 100*(Nhats[,"CV",1]-Nhats[,"CV",2])/Nhats[,"CV",1]
pcCVdiff

# Now look at estimates, CIs and CVs in high-stdGC region for best and const-site-D model
# -------------------------------------------------------------------------------
quartz(h=6,w=3)
par(mfrow=c(3,1))
breaks = seq(-1.5,5,length=31)
hist(covariates(fit$mask[[1]])$stdGC,breaks=breaks)
hist(covariates(fit$mask[[2]])$stdGC,breaks=breaks)
hist(covariates(fit$mask[[3]])$stdGC,breaks=breaks)

stdcut = 2 # stdGC cutpoint for inclusion below
maskeep = list(Tost=which(covariates(fit$mask[[1]])$stdGC>=stdcut),
               Noyon=which(covariates(fit$mask[[2]])$stdGC>=stdcut),
               Nemegt=which(covariates(fit$mask[[3]])$stdGC>=stdcut))
camkeep = list(Tost=which(covariates(traps(fit$capthist[[1]]))$stdGC>=stdcut),
               Noyon=which(covariates(traps(fit$capthist[[2]]))$stdGC>=stdcut),
               Nemegt=which(covariates(traps(fit$capthist[[3]]))$stdGC>=stdcut))
# uncomment below if want regions BELOW stdcut
#maskeep = list(Tost=which(covariates(fit$mask[[1]])$stdGC<stdcut),
#            Noyon=which(covariates(fit$mask[[2]])$stdGC<stdcut),
#            Nemegt=which(covariates(fit$mask[[3]])$stdGC<stdcut))
subcams = predmask = vector(length=3,mode="list")
nsubcams = array(rep(NA,3*3),dim=c(3,3),dimnames=list(Region=c("Tost","Noyon","Nemegt"),Set=c("All","subset","%kept")))
for(i in 1:3) {
  predmask[[i]] = subset(fit$mask[[i]],maskeep[[i]])
  subcams[[i]] = subset(traps(fit$capthist[[i]]),camkeep[[i]])
  nsubcams[i,1] = dim(traps(fit$capthist[[i]]))[1]
  nsubcams[i,2] = dim(subcams[[i]])[1]
}
nsubcams[,3] = round(100*nsubcams[,2]/nsubcams[,1],1)
nsubcams

sub.p = vector(length=3,mode="list")
for(i in 1:3) {
  Dsum = sum(covariates(D[[i]])$Dhat)
  cams = subcams[[i]]
  ncams = dim(cams)[1]
  Topo = covariates(cams)$Topo
  Water = covariates(cams)$Water
  # create the session factor for this region (needed for predict() to work, 
  # as seesion is in original fitted model)
  sess = factor(as.character(i),levels=as.character(1:3))
  # predict using a data frame with a row per trap, with associated covariates
  # stdGC has to be in because it is in the fitted model, but we don't use it
  # as we extract only the detection function parameters from predpar
  predpar = predict(fit,newdata=data.frame(stdGC=rep(0,ncams),Topo=Topo,Water=Water,session=rep(sess,ncams)))
  # Calulate distances from all traps to all mask points
  dists = distances(cams,fit$mask[[i]]) # distances from traps to each activity centre
  # Get expected number of encouners at each trap using detection parameters 
  # (a0 and sigma in this case) at each trap. These parameters have been 
  # obtained by predict() above, using the appropriate covariates at each trap.
  en = er_from_a0sig(predpar,dists,usage(cams)) # expected number of encouters
  entot = apply(en,2,sum) # expected number of encounters at each mask point, aggregating en across traps
  sub.p[[i]] = 1 - exp(-entot) # p at each mask point
  covariates(D[[i]])$sub.er = entot
  covariates(D[[i]])$sub.p = sub.p[[i]] # p(see|location)
  covariates(D[[i]])$sub.cover = sub.p[[i]]/dim(D[[i]])[1] # space coverage prob
  covariates(D[[i]])$sub.cover.pop = sub.p[[i]]*covariates(D[[i]])$Dhat/Dsum # pop coverage prob
}





xlim = range(D[[1]]$x,D[[2]]$x,D[[3]]$x)
ylim = range(D[[1]]$y,D[[2]]$y,D[[3]]$y)
plim = range(covariates(D[[1]])$sub.p,covariates(D[[2]])$sub.p,covariates(D[[3]])$sub.p)
#quartz(h=4,w=9)
pdf(file="./Analysis4paper/Plots/CoverageProbSubsetCams.pdf",h=4,w=9)
layout(mat=matrix(1:2,nrow=1),widths=c(.9,.1))
splotcovariate(D[[1]],"sub.p",what="image",zlim=plim,xlim=xlim,ylim=ylim,col=parula(40),
               main=paste("Coverage probability (camera stdGC>",stdcut,")",sep=""),cex.main=0.75)
plot(subcams[[1]],add=TRUE)
splotcovariate(D[[2]],"sub.p",what="image",zlim=plim,add=TRUE,col=parula(40))
plot(subcams[[2]],add=TRUE)
splotcovariate(D[[3]],"sub.p",what="image",zlim=plim,add=TRUE,col=parula(40))
plot(subcams[[3]],add=TRUE)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
splotcovariate(D[[1]],"sub.p",what="scale",zlim=plim,col=parula(40))
dev.off()


# Plot population coverage probability:
xlim = range(D[[1]]$x,D[[2]]$x,D[[3]]$x)
ylim = range(D[[1]]$y,D[[2]]$y,D[[3]]$y)
plim = range(covariates(D[[1]])$sub.cover.pop,covariates(D[[2]])$sub.cover.pop,covariates(D[[3]])$sub.cover.pop)
#quartz(h=4,w=9)
pdf(file="./Analysis4paper/Plots/PopCoverageProbSubsetCams.pdf",h=4,w=9)
layout(mat=matrix(1:2,nrow=1),widths=c(.9,.1))
splotcovariate(D[[1]],"sub.cover.pop",what="image",zlim=plim,xlim=xlim,ylim=ylim,col=parula(40),
               main=paste("Population coverage probability (camera stdGC>",stdcut,")",sep=""),cex.main=0.75)
plot(traps(fit$capthist[[1]]),add=TRUE)
splotcovariate(D[[2]],"sub.cover.pop",what="image",zlim=plim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[2]]),add=TRUE)
splotcovariate(D[[3]],"sub.cover.pop",what="image",zlim=plim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[3]]),add=TRUE)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
splotcovariate(D[[1]],"sub.cover.pop",what="scale",zlim=plim,col=parula(40))
dev.off()



# Plot spatial coverage probability:
xlim = range(D[[1]]$x,D[[2]]$x,D[[3]]$x)
ylim = range(D[[1]]$y,D[[2]]$y,D[[3]]$y)
plim = range(covariates(D[[1]])$sub.cover,covariates(D[[2]])$sub.cover,covariates(D[[3]])$sub.cover)
#quartz(h=4,w=9)
pdf(file="./Analysis4paper/Plots/SpaceCoverageProbSubsetCams.pdf",h=4,w=9)
layout(mat=matrix(1:2,nrow=1),widths=c(.9,.1))
splotcovariate(D[[1]],"sub.cover",what="image",zlim=plim,xlim=xlim,ylim=ylim,col=parula(40),
               main=paste("Space coverage probability (camera stdGC>",stdcut,")",sep=""),cex.main=0.75)
plot(traps(fit$capthist[[1]]),add=TRUE)
splotcovariate(D[[2]],"sub.cover",what="image",zlim=plim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[2]]),add=TRUE)
splotcovariate(D[[3]],"sub.cover",what="image",zlim=plim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[3]]),add=TRUE)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
splotcovariate(D[[1]],"sub.cover",what="scale",zlim=plim,col=parula(40))
dev.off()


# Plot subsetted stdGC in space (need stdGCrange from abaove)
xlim = range(fit$mask[[1]]$x,fit$mask[[2]]$x,fit$mask[[3]]$x)
ylim = range(fit$mask[[1]]$y,fit$mask[[2]]$y,fit$mask[[3]]$y)
quartz(h=4,w=9)
layout(mat=matrix(1:2,nrow=1),widths=c(.9,.1))
splotcovariate(predmask[[1]],"stdGC",what="image",zlim=stdGCrange,xlim=xlim,ylim=ylim,col=parula(40),
               main="stdGC",cex.main=0.75)
plot(subcams[[1]],add=TRUE)
splotcovariate(predmask[[2]],"stdGC",what="image",zlim=stdGCrange,add=TRUE,col=parula(40))
plot(subcams[[2]],add=TRUE)
splotcovariate(predmask[[3]],"stdGC",what="image",zlim=stdGCrange,add=TRUE,col=parula(40))
plot(subcams[[3]],add=TRUE)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
splotcovariate(predmask[[1]],"stdGC",zlim=stdGCrange,what="scale",col=parula(40))

# Look at total abundance estimates, and CIs and CVs for best and const-site-D model:
# IN selected subset of stdGC REGIONS ONLY
# -------------------------------------------------------------------------------
Nbest = region.N(fit,region=predmask)
Nsite = region.N(fitsite,region=predmask)
Nhats = array(rep(NA,5*3*2),dim=c(3,8,2),
                        dimnames=list(region=c("Tost","Noyon","Nemegt"),
                                      est=c("point","SE","lcl","ucl","n","CV","empBiasPC","RMSE"),
                                      model=c("best","const-site-D")))
for(i in 1:3) {
  # add CV column and "bias" column
  Nhats[i,,1] = c(as.numeric(Nbest[[i]]["E.N",]),Nbest[[i]]["E.N",2]/Nbest[[i]]["E.N",1],0,
                  Nbest[[i]]["E.N",2]^2)
  Nhats[i,,2] = c(as.numeric(Nsite[[i]]["E.N",]),Nsite[[i]]["E.N",2]/Nsite[[i]]["E.N",1],
                  100*(Nbest[[i]]["E.N",1]-Nsite[[i]]["E.N",1])/Nbest[[i]]["E.N",1],
                  Nsite[[i]]["E.N",2]^2+((Nbest[[i]]["E.N",1]-Nsite[[i]]["E.N",1]))^2)
}
Nhats















#-------------------------------------------------------------------------------
# Simulations with trap-specific covariates
#-------------------------------------------------------------------------------
Dhats.core = predictDsurface(fit) # Densities to use in simulating population
Nsim = 2 # number of simulaions (these take quite a LONG time)
# object to hold estimated densities:
#Nest = array(rep(NA,Nsim*3*3),dim=c(Nsim,3,3),
#             dimnames=list(sim=1:Nsim,region=c("Tost","Noyon","Nemegt"),method=c("null","DstdGC","full")))
#Nest = array(rep(NA,Nsim*3*2),dim=c(Nsim,3,2),
#             dimnames=list(sim=1:Nsim,region=c("Tost","Noyon","Nemegt"),method=c("cosnt-site-D","full")))
Nest = n = array(rep(NA,Nsim*3),dim=c(Nsim,3),
             dimnames=list(sim=1:Nsim,region=c("Tost","Noyon","Nemegt")))
subsetcams = FALSE 
subsetcams = TRUE # Set to FALSE if want to use all cameras, not selected subset
subsetmult = 100/nsubcams[,"%kept"]
scaleup = rep(NA,3)
# Do the simulations
# Set up progress bar, before looping:
pb <- tkProgressBar(title=paste("Simulation Progress (Nsim=",Nsim,")",sep=""), min=0, max=Nsim, width=400)
for(sim in 1:Nsim) {
  chdf=data.frame(session=NULL,ID=NULL,occasion=NULL,trap=NULL) # data frame for capture histories
  for(reg in 1:3) { # we have 3 regions, need to simulate separately in each
    if(subsetcams) {
      cams = subcams[[reg]]
      # scale up density in inverse proportion to % cameras kept
      if(sim==1) {
        scaleup[reg] = sum(covariates(D[[reg]])$cover.pop)/sum(covariates(D[[reg]])$sub.cover.pop)
        covariates(Dhats.core[[reg]])$D.0 = covariates(Dhats.core[[reg]])$D.0 * scaleup[reg]
      }
    }
    else cams = traps(fit$capthist[[reg]])
    ntraps = dim(cams)[1]
    # get predictor variables for the traps in this region
    Topo = covariates(cams)$Topo
    Water = covariates(cams)$Water
    # create the session factor for this region (needed for predict() to work, 
    # as seesion is in original fitted model)
    sess = factor(as.character(reg),levels=as.character(1:3))
    # predict using a data frame with a row per trap, with associated covariates
    # stdGC has to be in because it is in the fitted model, but we don't use it
    # as we extract only the detection function parameters from predpar
    predpar = predict(fit,newdata=data.frame(stdGC=rep(0,ntraps),Topo=Topo,Water=Water,session=rep(sess,ntraps)))
    # simulate a population with required density surface
    pop = sim.popn("D.0",core=Dhats.core[[reg]],model2D="IHP",Ndist="fixed")
    # Now need to generate capture histories, which depend on distances and on
    # trap covariates, so need distances
    dists = distances(cams,pop) # distances from traps to each activity centre
    nindiv = dim(pop)[1] # number of activity centres
    # Get expected number of encouners at each trap using detection parameters 
    # (a0 and sigma in this case) at each trap. These parameters have been 
    # obtained by predict() above, using the appropriate covariates at each trap.
    en = er_from_a0sig(predpar,dists,usage(cams)) # expected number of encouters
    trapname = row.names(cams) # preserve trap names
    # now for each individual, generate (Poisson) number encounters at each trap,
    # using en calculated above, and add these to the capture history data frame, chdf
    for(i in 1:nindiv) {
      cnt = rpois(ntraps,en[,i])
      for(j in 1:ntraps) {
        if(cnt[j]>0) for(k in 1:cnt[j]) 
          chdf = rbind(chdf,data.frame(session=reg,ID=i,occasion=1,trap=trapname[j]))
      }
    }
  }
  simch = make.capthist(chdf,traps(fit$capthist)) # turn chdf into a capthist object
  for(reg in 1:3) n[sim,reg] = dim(simch[[reg]])[1]
  
  # Now fit whatever models you want, to this new capture history:
#  testfit.DGrid.a0Waterxsess_topo.sig1<-secr.fit(simch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
#                                                 model=list(D~stdGC, a0~Topo+Water*session, sigma~1), ncores=ncores)
#  testfit.DGrid<-secr.fit(simch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
#                          model=list(D~stdGC, a0~1, sigma~1), ncores=ncores)
#  testfit0<-secr.fit(simch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
#                     model=list(D~1, lambda0~session, sigma~1), ncores=ncores)
  testfit.sess<-secr.fit(simch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                     model=list(D~session, lambda0~session, sigma~session), ncores=ncores)
  # store the realised abundance estimate for each model
#  Nest[sim,,1] = R.N(testfit0)
#  Nest[sim,,1] = R.N(testfit.sess)
  Nest[sim,] = R.N(testfit.sess)
  #  Nest[sim,,2] = R.N(testfit.DGrid)
#  Nest[sim,,3] = R.N(testfit.DGrid.a0Waterxsess_topo.sig1)
#  Nest[sim,,2] = R.N(testfit.DGrid.a0Waterxsess_topo.sig1)
  
  # Progress bar update
  setTkProgressBar(pb, sim, label=paste( round(sim/Nsim*100, 0),"% done"))
  
}
# Close progress bar
close(pb)

Nest

n

true.N = A = rep(NA,3)
for(i in 1:3) {
  A[i] = attr(fit$mask[[i]],"area")*dim(fit$mask[[i]])[1]
  true.N[i] = mean(covariates(Dhats.core[[i]])$D.0*A[i])
}
mean.Nest = apply(Nest,2,mean)
100*(mean.Nest - true.N)/true.N



