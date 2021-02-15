load("Analysis4paper/nonEuc-mdls.RData")

# best model:
fit = TNNfit.DGrid.lam0Waterxsess_topo.sig3
# second best model:
constRfit = TNNfit.DGrid.lam0Waterxsess_topo.sig0
Dhatconst = exp(coefficients(TNNfit.DGrid.lam0Waterxsess_topo.sig0)[1,1])
  
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
Dmean = D = D.lcl = D.ucl = stdGC = data.frame(Tost=rep(NA,100),Noyon=rep(NA,100),Nemgt=rep(NA,100))
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

# Plot stdGC in space
xlim = range(fit$mask[[1]]$x,fit$mask[[2]]$x,fit$mask[[3]]$x)
ylim = range(fit$mask[[1]]$y,fit$mask[[2]]$y,fit$mask[[3]]$y)
pdf(file="./Analysis4paper/Plots/trapsTopo.pdf",h=4,w=9)
layout(mat=matrix(1:2,nrow=1),widths=c(.9,.1))
splotcovariate(fit$mask[[1]],"stdGC",what="image",zlim=stdGCrange,xlim=xlim,ylim=ylim,col=terrain.colors(40),
               main="",cex.main=0.75)
i=1
cams = traps(fit$capthist)[[i]]
tcol=rep(0,dim(cams)[1])
tcol[covariates(cams)$Topo=="Saddle"] = "red"
tcol[covariates(cams)$Topo=="Canyon"] = "black"
tcol[covariates(cams)$Topo=="Steppe"] = "green"
wshp=rep(0,dim(cams)[1])
wshp[covariates(cams)$Water=="Yes"] = 20
wshp[covariates(cams)$Water=="No"] = 18
points(cams$x,cams$y,col=tcol,pch=wshp)
splotcovariate(fit$mask[[2]],"stdGC",what="image",zlim=stdGCrange,add=TRUE,col=terrain.colors(40))
i=2
cams = traps(fit$capthist)[[i]]
tcol=rep(0,dim(cams)[1])
tcol[covariates(cams)$Topo=="Saddle"] = "red"
tcol[covariates(cams)$Topo=="Canyon"] = "black"
tcol[covariates(cams)$Topo=="Steppe"] = "green"
wshp=rep(0,dim(cams)[1])
wshp[covariates(cams)$Water=="Yes"] = 20
wshp[covariates(cams)$Water=="No"] = 18
points(cams$x,cams$y,col=tcol,pch=wshp)
splotcovariate(fit$mask[[3]],"stdGC",what="image",zlim=stdGCrange,add=TRUE,col=terrain.colors(40))
legend("topleft",bty="n",pch=rep(19,3),col=c("red","black","green"),legend=c("Saddle","Canyon","Steppe"),cex=0.75)
i=3
cams = traps(fit$capthist)[[i]]
tcol=rep(0,dim(cams)[1])
tcol[covariates(cams)$Topo=="Saddle"] = "red"
tcol[covariates(cams)$Topo=="Canyon"] = "black"
tcol[covariates(cams)$Topo=="Steppe"] = "green"
wshp=rep(0,dim(cams)[1])
wshp[covariates(cams)$Water=="Yes"] = 20
wshp[covariates(cams)$Water=="No"] = 18
points(cams$x,cams$y,col=tcol,pch=wshp)
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
splotcovariate(fit$mask[[1]],"stdGC",what="scale",col=terrain.colors(40))
dev.off()



# Plot density against stdGC
# --------------------------
for(j in 1:3) {
  for(i in 1:100) {
    D[i,j] = Dhat[[j]][[i]][1,2]*10^4
    D.lcl[i,j] = Dhat[[j]][[i]][1,4]*10^4
    D.ucl[i,j] = Dhat[[j]][[i]][1,5]*10^4
    Dmean[,j] = rep(mean.Dhat[j]*10^4,100)
  }
}
Dhatconst = Dhatconst*10^4
Dmax = 0
for(i in 1:3) Dmax=max(Dmax,c(min(D.lcl[,i]),2*max(D[,i])))
pdf(file="./Analysis4paper/Plots/DensityVsSdGC.pdf",h=3,w=11)
par(mfrow=c(1,3))
par(mar=c(4,5,2,2))
for(i in 1:3) {
  plot(stdGC[,i],D[,i],type="l",ylim=c(0,Dmax),main=names(D)[i],xlab="stdGC",ylab=expression(hat(D)))
  shadedCI(data.frame(stdGC[,i],D.lcl[,i],D.ucl[,i]))
  lines(stdGC[,i],rep(Dhatconst,length(stdGC[,i])),col="black",lty=2)
  trapcov = covariates(traps(fit$capthist[[i]]))$stdGC
  points(trapcov,rep(min(D.lcl[,i]),length(trapcov)),pch="|")
}
dev.off()


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

# Density maps
# -----------------
xlim = range(D[[1]]$x,D[[2]]$x,D[[3]]$x)
ylim = range(D[[1]]$y,D[[2]]$y,D[[3]]$y)
Dlim=c(Dmin,Dmax)
#quartz(h=4,w=9)
pdf(file="./Analysis4paper/Plots/Dhat.pdf",h=4,w=9)
layout(mat=matrix(1:2,nrow=1),widths=c(.9,.1))
splotcovariate(D[[1]],"Dhat",what="image",zlim=Dlim,xlim=xlim,ylim=ylim,col=parula(40),
               main="",cex.main=0.75)
plot(traps(fit$capthist[[1]]),add=TRUE,detpar=list(col="darkgray",cex=0.5))
splotcovariate(D[[2]],"Dhat",what="image",zlim=Dlim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[2]]),add=TRUE,detpar=list(col="darkgray",cex=0.5))
splotcovariate(D[[3]],"Dhat",what="image",zlim=Dlim,add=TRUE,col=parula(40))
plot(traps(fit$capthist[[3]]),add=TRUE,detpar=list(col="darkgray",cex=0.5))
# add region labels
text(680000,4825000,labels="Nemegt")
text(730000,4805000,labels="Noyon")
text(600000,4798000,labels="Tost")
splotcovariate(D[[1]],"Dhat",what="scale",col=parula(40))
dev.off()

# Density CV maps
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

# Density lcl and ucl maps on same scale
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

# Density ucl maps
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
# Need function to calculate expected encouner rate or detection probability 
# from list of parameters with a0 and sigma
er_from_a0lamsig = function(parlist,dists,usage=rep(1,length(parlist)),calcprob=FALSE) {
  if(!inherits(parlist,"list")) stop("'parlist' must be a list.")
  if(is.element("a0",row.names(parlist[[1]]))) partype = "a0"
  else if(is.element("lambda0",row.names(parlist[[1]]))) partype = "lambda0"
  else stop("parlist elements must have row names that includ either 'a0' or 'lambda0'.")
  npars = length(parlist)
  if(is.vector(dists)) dists = matrix(dists,ncol=1)
  if(!is.matrix(dists)) stop("dists must be a vector or a matrix.")
  if(dim(dists)[1]!=npars) stop("dists must be of same length as npars, or if matrix, have npars rows.")
  ndists = dim(dists)[2]
  er = matrix(rep(NA,npars*ndists),nrow=npars)
  for(i in 1:npars) {
    sigma = parlist[[i]]["sigma","estimate"]
    if(partype=="a0") {
      a0 = parlist[[i]]["a0","estimate"] * 10^4 # conversion from hectaers to sq m
      lambda0 = a0/(2*pi*sigma^2)
    } else {
      lambda0 = parlist[[i]]["lambda0","estimate"]
    }
    er[i,] = lambda0*exp(-dists[i,]^2/(2*sigma^2)) * usage[i]
  }
  if(calcprob) {
    p = 1 - exp(-er)
    return(p)
  } else {
    return(er)
  }
}


#noneucPred = predictDsurface(fit,parameter="noneuc")
p = vector(length=3,mode="list")
for(i in 1:3) {
  cams = traps(fit$capthist[[i]]) # extract cameras for this site
  ncams = dim(cams)[1]
  # make predictor variables for the site, and put into a data frame for prediction
  Topo = covariates(cams)$Topo
  Water = covariates(cams)$Water
#  noneuc = covariates(noneucPred[[i]])$noneuc.0
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
  en = er_from_a0lamsig(predpar,dists,usage(cams)) # expected number of encouters
  entot = apply(en,2,sum) # expected number of encounters at each mask point, aggregating en across traps
  p[[i]] = 1 - exp(-entot) # p at each mask point
  covariates(D[[i]])$er = entot
  covariates(D[[i]])$p = p[[i]] # p(see|location)
  covariates(D[[i]])$cover = p[[i]]/dim(D[[i]])[1] # space coverage prob
  Dsum = sum(covariates(D[[i]])$Dhat)
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
               main="",cex.main=0.75)
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


