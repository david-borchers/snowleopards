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
er_from_a0sig = function(parlist,dists,usage=rep(1,length(parlist))) {
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
  return(er)
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
#fit = TNNfit.DGrid.a0Waterxsess.sig1

# best model among those with flat density
fit0 = TNNfit.D1.LamSess.sig1

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

Dhat = vector(length=3,mode="list")
Dmean = D0 = D = D.lcl = D.ucl = stdGC = data.frame(Tost=rep(NA,100),Noyon=rep(NA,100),Nemgt=rep(NA,100))
for(i in 1:3) {
  stdGCrange = range(covariates(fit$mask[[i]])$stdGC)
  stdGC[,i] = seq(stdGCrange[1],stdGCrange[2],length=100)
  newdata = data.frame(
    stdGC=stdGC[,i],
    Topo = rep(covariates(traps(fit$capthist[[1]]))$Topo[1],100),  # [1] is arbitrary value
    Water = rep(covariates(traps(fit$capthist[[1]]))$Water[1],100),  # [1] is arbitrary value
    session = rep(factor(as.character(i),levels=c("1","2","3")),100)
  )
  Dhat[[i]] = predict(fit,newdata=newdata)
}
Dhat0 =  predict(TNNfit.D1.LamSess.sig1)
pclow = rep(NA,3)
for(j in 1:3) {
  pclow[j] = 100*(Dhat0[[j]][1,2]-mean.Dhat[j])/Dhat0[[j]][1,2]
  for(i in 1:100) {
    D[i,j] = Dhat[[j]][[i]][1,2]*10^5
    D.lcl[i,j] = Dhat[[j]][[i]][1,4]*10^5
    D.ucl[i,j] = Dhat[[j]][[i]][1,5]*10^5
    D0[,j] = rep(Dhat0[[j]][1,2]*10^5,100)
    Dmean[,j] = rep(mean.Dhat[j]*10^5,100)
  }
}

pclow # % by which constant density model density is below mean of best model density


# Plot density against stdGC
quartz(h=3,w=11)
par(mfrow=c(1,3))
par(mar=c(4,5,2,2))
for(i in 1:3) {
  plot(stdGC[,i],D[,i],type="l",ylim=c(min(D.lcl[,i]),2*max(D[,i])),main=names(D)[i],xlab="stdGC",ylab=expression(hat(D)))
  shadedCI(data.frame(stdGC[,i],D.lcl[,i],D.ucl[,i]))
  lines(stdGC[,i],Dmean[,i],col="darkgray")
  lines(stdGC[,i],D0[,i],lty=2,col="darkgray")
  trapcov = covariates(traps(fit$capthist[[i]]))$stdGC
  points(trapcov,rep(1,length(trapcov)),pch="|")
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
legend("topleft",bty="n",lty=rep(1,3),col=c("red","black","green"),legend=c("Saddle","Canyon","Steppe"),cex=0.75)
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
legend("topleft",bty="n",lty=rep(1,2),col=c("red","black"),legend=c("Yes","No"),cex=0.75)
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

tab = vector(length=3,mode="list")
for(i in 1:3) tab[[i]] = table(covariates(traps(fit$capthist[[i]]))[,c("Topo","Water")])
tab

# Plot Density in space
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


#-------------------------------------------------------------------------------
# Simulations with trap-specific covariates
#-------------------------------------------------------------------------------
Dhats.core = predictDsurface(fit) # Densities to use in simulating population
Nsim = 30 # number of simulaions (these take quite a LONG time)
# object to hold estimated densities:
Nest = array(rep(NA,Nsim*3*3),dim=c(Nsim,3,3),
             dimnames=list(sim=1:Nsim,region=c("Tost","Noyon","Nemegt"),method=c("null","DstdGC","full")))
# Do the simulations
# Set up progress bar, before looping:
pb <- tkProgressBar(title=paste("Simulation Progress (Nsim=",Nsim,")",sep=""), min=0, max=Nsim, width=400)
for(sim in 1:Nsim) {
  chdf=data.frame(session=NULL,ID=NULL,occasion=NULL,trap=NULL) # data frame for capture histories
  for(reg in 1:3) { # we have 3 regions, need to simulate separately in each
    cams = traps(fit$capthist[[reg]])
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
    pop = sim.popn("D.0",core=Dhats.core[[reg]],model2D="IHP")
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
  
  # Now fit whatever models you want, to this new capture history:
  testfit.DGrid.a0Waterxsess_topo.sig1<-secr.fit(simch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                                 model=list(D~stdGC, a0~Topo+Water*session, sigma~1), ncores=ncores)
  testfit.DGrid<-secr.fit(simch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                          model=list(D~stdGC, a0~1, sigma~1), ncores=ncores)
  testfit0<-secr.fit(simch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                     model=list(D~1, a0~1, sigma~1), ncores=ncores)
  # store the realised abundance estimate for each model
  Nest[sim,,1] = R.N(testfit0)
  Nest[sim,,2] = R.N(testfit.DGrid)
  Nest[sim,,3] = R.N(testfit.DGrid.a0Waterxsess_topo.sig1)
  
  # Progress bar update
  setTkProgressBar(pb, i, label=paste( round(i/Nsim*100, 0),"% done"))
  
}
# Close progress bar
close(pb)

apply(Nest,c(2,3),mean) # calculate mean realised abundances for each method



