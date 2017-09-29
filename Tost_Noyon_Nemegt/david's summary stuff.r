load("./Tost/Tost-nonEuc-fitsx.RData")

AICTost=AIC(Tost.hhnx, Tost.hhn.detTopo10x, Tost.hhn.detWaterx, Tost.hhn.DHabx, Tost.hhn.DHab.nonUx, 
            Tost.hhn.DHab.nonU.Topo10x,Tost.hhn.DHab.nonU.Wx, Tost.hhn.DHab.nonU.T01Wx, 
            Tost.hhn.D.nonUx, Tost.hhn.DHab.nonU.GBx, Tost.hhn.DHab.nonU.GBGCx)

AICTost
coefficients(Tost.hhn.DHab.nonU.GBx)
region.N(Tost.hhn.DHab.nonU.GBx)
FXTost<-fx.total(Tost.hhnx)
FXTost<-fx.total(Tost.hhn.DHab.nonU.GBx)
# Plot using plot.Dsurface from secr 
# (Note scrplotting.r also has a function called plot.Dsurface, but it has diffent arguments.)
secr:::plot.Dsurface(FXTost, covariate = 'D.sum', breaks = seq(0,10e-5,1e-5), poly = FALSE)
plot(Tost.hhn.DHab.nonU.GBx$capthist,add=TRUE,tracks=TRUE)
plot(traps(Tost.hhn.DHab.nonU.GBx$capthist),add=TRUE)

# Plot using plotcovariate from secrplotting.r
plotcovariate(FXTost, covariate = 'D.sum',xaxt="n",yaxt="n",xlab="",ylab="",
              bty="n",contour=FALSE,asp=1)
plotMaskEdge(Tost.hhn.DHab.nonU.GBx$mask,add=TRUE)
plot(Tost.hhn.DHab.nonU.GBx$capthist,add=TRUE,tracks=TRUE)
plot(traps(Tost.hhn.DHab.nonU.GBx$capthist),add=TRUE)
# Compare to density or log(density) plot:
TostDhat = predictDsurface(Tost.hhn.DHab.nonU.GBx)
covariates(TostDhat)$logD.0 = log(covariates(TostDhat)$D.0)
plotcovariate(TostDhat, covariate = 'D.0',xaxt="n",yaxt="n",xlab="",ylab="",
              bty="n",contour=FALSE,asp=1)
plotcovariate(TostDhat, covariate = 'logD.0',xaxt="n",yaxt="n",xlab="",ylab="",
              bty="n",contour=FALSE,asp=1)
plotMaskEdge(Tost.hhn.DHab.nonU.GBx$mask,add=TRUE)
plot(Tost.hhn.DHab.nonU.GBx$capthist,add=TRUE,tracks=TRUE)
plot(traps(Tost.hhn.DHab.nonU.GBx$capthist),add=TRUE)



Tost_top<-region.N(Tost.hhn.DHab.nonU.GBx)
Tost_null<-region.N(Tost.hhnx)
Tost_Diff=(Tost_null[2,1]-Tost_top[2,1])*100/Tost_null[2,1]
Tost_Diff
summary(TostMask1)

D_Tost<-(Tost_top*100/2116.50)

# Mean density Tost:
Tostarea = dim(Tost.hhn.DHab.nonU.GBx$mask)[1]* summary(Tost.hhn.DHab.nonU.GBx$mask)$cellarea
TostmeanD = Tost_top["R.N",]/Tostarea

load("./Noyon2013/Noyon-nonEuc-fitsx.RData") #Final round analysis (with water and topography)
AICNoyon=AIC(Noyon.hhnx, Noyon.hhn.detrgdx, Noyon.hhn.DHabx, Noyon.hhn.DHab.DetRgd01x, Noyon.hhn.detWaterx, 
             Noyon.hhn.detTopo10x, Noyon.hhn.detTopoWaterx, Noyon.hhn.DHab.nonUx, Noyon.hhn.D.nonUx, 
             Noyon.hhn.DHab.Topo10.nonUx, Noyon.hhn.DHab.DetW.nonUx, Noyon.hhn.DHab.Topo10W.nonUx,
             Noyon.hhn.DHab.nonU.GBx, Noyon.hhn.DHab.nonU.GBGCx, Noyon.hhn.DHab.DetRgd10x)
AICNoyon
coefficients(Noyon.hhn.DHab.nonU.GBx)
region.N(Noyon.hhn.DHab.nonU.GBx)

Noyon_Null<-region.N(Noyon.hhnx)
Noyon_top<-region.N(Noyon.hhn.DHab.nonU.GBx)
Noyon_Diff=(Noyon_Null[2,1]-Noyon_top[2,1])*100/Noyon_Null[2,1]
Noyon_Diff

summary(NoyonMask1)

D_Noyon<-(Noyon_top*100/2481.75)

# Mean density Noyon:
Noyarea = dim(Noyon.hhn.DHab.nonU.GBx$mask)[1]* summary(Noyon.hhn.DHab.nonU.GBx$mask)$cellarea
NoyonmeanD = Noyon_top["R.N",]/Noyarea


load("./Nemegt/Nemegt-nonEuc-fit2xR.RData")
NemegtAIC2xR=AIC(Nemegt.hhn2xR, Nemegt.hhn.detrgd2xR, Nemegt.hhn.DHab2xR, Nemegt.hhn.DHab.detrgd102xR, 
                 Nemegt.hhn.DHab.detrgd012xR, Nemegt.hhn.DHab.nonU2xR, Nemegt.hhn.D.nonU2xR,  
                 Nemegt.hhn.DHab.nonU.GB2xR, Nemegt.hhn.DHab.nonU.LamTopoR, Nemegt.hhn.DHab.nonU.LamW2xR,
                 Nemegt.hhn.DHab.nonU.LamTopoW2R)
NemegtAIC2xR
coefficients(Nemegt.hhn.DHab.nonU.LamW2xR)
Nemegt_Top<-region.N(Nemegt.hhn.DHab.nonU.LamW2xR)
Nemegt_Null<-region.N(Nemegt.hhn2xR)
Nemegt_Diff=(Nemegt_Null[2,1]-Nemegt_Top[2,1])*100/Nemegt_Null[2,1]
Nemegt_Diff
summary(NemegtMask1)

D_Nemegt<-(Nemegt_Top*100/2145.50)


# Mean density Nemegt:
Nemarea = dim(Nemegt.hhn.DHab.nonU.LamW2xR$mask)[1]* summary(Nemegt.hhn.DHab.nonU.LamW2xR$mask)$cellarea
NemegtmeanD = Nemegt_Top["R.N",]/Nemarea

# Average densities:
meanD = rbind(TostmeanD,NoyonmeanD,NemegtmeanD)
row.names(meanD) = c("Tost","Noyon","Nemegt")
meanD[,1:4] = meanD[,1:4]*10^4
meanD
ylim=range(meanD[,c("estimate","lcl","ucl")])
plot(as.factor(row.names(meanD)),meanD$estimate,ylim=ylim,ylab=expression(hat(D)))
arrows(3:1,meanD$lcl,3:1,meanD$ucl,angle=90,code=3)
plot(as.factor(row.names(meanD)),log(meanD$estimate),ylim=log(ylim),ylab=expression(log(hat(D))))
arrows(3:1,log(meanD$lcl),3:1,log(meanD$ucl),angle=90,code=3)

pdf("Dhats.pdf",w=8,h=3.5)
par(mar=c(3, 5, 2, 2) + 0.1)
plot(as.factor(row.names(meanD)),log(meanD$estimate),ylim=log(ylim),ylab=expression(log(hat(D))))
arrows(3:1,log(meanD$lcl),3:1,log(meanD$ucl),angle=90,code=3)
dev.off()


D_Tost
D_Noyon
D_Nemegt

#This is only where non Euc goes marginally negative!

load("./Tost_Noyon_Nemegt/TNN-NonEuc-fits2xR.RData")

TNNAIC2xR<- AIC(TNN.hhn.DRgd.sessR, TNN.hhn.DRgd.sess.DetWR, TNN.hhn.DRgd.DetTopo10WR, TNN.hhn.DHab.nonUR, 
                TNN.hhn.DHab.DetTopo10.nonUR, TNN.hhn.DHab.DetToposess.nonUR, 
                TNN.hhn.DHab.S.DetTopo10.nonUR, TNN.hhn.DHab_S.DetTopo10.nonUR, 
                TNN.hhn.DHab.LamTopoWat.nonUR, TNN.hhn.nonUR,TNN.hhn.DGC.DetTopo10.nonUGBR,
                TNN.hhn.DGB.DetTopo10.nonUR)
TNNAIC2xR
coefficients(TNN.hhn.DGC.DetTopo10.nonUGBR)
coefficients(TNN.hhn.DHab.DetTopo10.nonUR)
coefficients(TNN.hhn.DGB.DetTopo10.nonUR)

TNN.hhn.DHabS3.nonU
Tost.hhn.DHab.nonU.GBGCx
Noyon.hhn.DHab.nonUx
Nemegt.hhn.DHab.nonU.LamW2x

AIC(TNN.hhn.DHabS3.nonU)
AICTostx = read.csv("AICTostx.csv")
AICTostx
AICNoyonx = read.csv("AICNoyonx.csv")
AICNoyonx
AICNemegtx = read.csv("AICNemegtx.csv")
AICNemegtx

# Best separate AICs:
aics = c(465.603,436.753,286.066)
aicsum = sum(aics);aicsum # Combined AIC for models fitted separately to each stratum
AIC(TNN.hhn.DHabS3.nonU)

# number of parameters:
q1 = 6
q2 = 5
q3 = 6
q. = 6

n. = round(2*q.*(q.+1)/(1015.83-1013.496)+q.+1)
qs = q1+q2+q3; qs
AICc.all = aicsum + 2*qs*(qs+1)/(n.-qs-1)
AICc.all; AIC(TNN.hhn.DHabS3.nonU)

D0 =-8.7322280
D1 = -1.4867679
D2 = 1.5543137

range1 = range(covariates(TostMask1)$stdGC)
range2 = range(covariates(NoyonMask1)$stdGC)
range3 = range(covariates(NemegtMask1)$stdGC)
stdGC1 = seq(range1[1],range2[2],length=100)
stdGC2 = seq(range1[1],range2[2],length=100)
stdGC3 = seq(range1[1],range2[2],length=100)
stdGC = seq(min(range1,range2,range3),max(range1,range2,range3),length=100)



# To check effect of non-Euc:
# fit with fixed sigma:
test<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                               model=list(D~s(stdGC,k=3), lambda0~1, noneuc ~ stdGC -1), 
                               details = list(userdist = userdfn1), fixed=list(sigma=1)) 

# Fit with inverse noneuc fn:
userdfn1i <- function (xy1, xy2, mask) {
  if (missing(xy1)) return('noneuc') #When function is called, more of a jargon. Tells that it is a non-euclidean function
  require(gdistance) #to load transition and geoCorrection functions
  Sraster <- raster(mask, 'noneuc') #Creates a raster from a set of coordinates and attributes and turn that into a raster. noneuc needs to be made in advance in the mask that is being used in the analysis
  ## conductance is inverse of friction
#  trans <- transition(Sraster, transitionFunction = function(x) 1/mean(x),  directions = 16)
  trans <- transition(Sraster, transitionFunction = function(x) mean(x),  directions = 16)
  trans <- geoCorrection(trans) #takes care of earth's curvature and also the distance differences between square and diagonally neighbouring cells
  costDistance(trans, as.matrix(xy1), as.matrix(xy2))
}

TNN.hhn.DHabS3.nonUi<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                              model=list(D~s(stdGC,k=3), lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                              details = list(userdist = userdfn1i),
                              start = list(noneuc = 1)) #-1 gets rid of the intercept
# Fit again, using estimates from above as starting values:
TNN.hhn.DHabS3a.nonUi<-secr.fit(all.data.TNN, detectfn="HHN", mask=list(TostMask1, NoyonMask1,NemegtMask1), 
                               model=list(D~s(stdGC,k=3), lambda0~1, sigma~1, noneuc ~ stdGC -1), 
                               details = list(userdist = userdfn1i),
                               start = TNN.hhn.DHabS3.nonUi) 
TNN.hhn.DHabS3.nonUi = TNN.hhn.DHabS3a.nonUi

coefficients(TNN.hhn.DHabS3.nonU)
coefficients(TNN.hhn.DHabS3.nonUi)
