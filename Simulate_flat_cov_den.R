# Kruger Park leopards density and camera frequencey plots
# --------------------------------------------------------
# fitted model is in erhr.fit,  mask in mask.utm
#Do I use the model from TNanalysisNoNe.r or Tost estimates obtained independently?
predD0 = predict(er.hr.fit)["D",2:5]*100^2 # multiply by 100^2 to convert animals per hectare to per sq km
nhs = 200
hrange = range(covariates(mask.utm)$habitat.cov)
hs = read.traps(data=data.frame(x=1:nhs,y=1:nhs))
habitat.cov = seq(hrange[1],hrange[2],length=nhs)
covariates(hs)$habitat.cov = habitat.cov
hsmooth = covariates(predictDsurface(fit.hab.hhn,hs,se.D=TRUE,cl.D=TRUE))*100^2
ylim=range(hsmooth$lcl.0,hsmooth$ucl.0)

plot(habitat.cov,hsmooth$D.0,ylim=ylim,type="l",xlab="Habitat suitability covariate",ylab="Density")
lines(habitat.cov,hsmooth$lcl.0,lty=2)
lines(habitat.cov,hsmooth$ucl.0,lty=2)
lines(c(min(habitat.cov),max(habitat.cov)),rep(predD0$estimate,2),col="gray")

# Look at distribution of traps across true density
kruger.cams = addCovariates(kruger.cams,mask.utm,column="D",replace=TRUE)
Drange = range(covariates(kruger.cams)$D)

h =hist(covariates(kruger.cams)$D,breaks=seq(Drange[1],Drange[2],length=21),plot=FALSE)
barplot(c(rep(0,4),h$counts,rep(0,17)),space=0,horiz=TRUE,xlab="Number of cameras",ylab="Density")


predict(TNNfit.DGrid.a0Topo.sig1)["D"]
predD0 = predict(TNNfit.DGrid.a0Topo.sig1)["D",2:5]*100*100^2 # multiply by 100^2 to convert animals per hectare to per sq km

# Tost-TosonBumba-Nemegt snow leopards simulation code
# ------------------------------------
# Simulate with true density to check if flat density model is biased because of trap placement:
simcount = function(sim,marks=c(5,10,50)) {
  if(sim==1) cat("\nCounter: \n")
  if(sim %% marks[3] == 0) cat(sim,"\n")
  else if(sim %% marks[2] == 0) cat(sim)
  else if(sim %% marks[1] == 0) cat("+")
  else cat("-")
}
dpar = detectpar(Tost.hhnx)

### Do we need to call predict with the top model to get D before this step?
#trueD = mean(D)

nsim = 100
cDest = data.frame(Dhat=rep(NA,nsim),lambda0=rep(NA,nsim),sigma=rep(NA,nsim))
startt = date()
dbad = cbad = 0
for(i in 1:nsim) {
  simpop = sim.popn(D,TostMask, buffer=0, model2D="IHP")
  simch = sim.capthist(cams,popn=simpop,detectfn="HHN",detectpar=list(labmda0=dpar$lambda0,sigma=dpar$sigma),noccasions=1)
  cDfit = secr.fit(simch, mask=TostMask,detectfn="HHN", model = list(D~session)) #start=list(D=trueD,labmda0=dpar$lambda0,sigma=dpar$sigma),trace=FALSE)
  if(!any(is.na(hDfit$beta.vcv))) {
    cDest[i,] = exp(coefficients(cDfit)$beta)
  } else {
    cbad = cbad+1
  }
  simcount(i)
}
endt = date()

100*(mean(na.omit(cDest[,1]))-trueD)/trueD # Estimated bias of constant-D estimator
Dlim = range(,na.omit(cDest[,1]))
#use predict to compare constant density model with whatever the truth is!