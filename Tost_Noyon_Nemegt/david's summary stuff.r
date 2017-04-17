load("./Tost/Tost-nonEuc-fitsx.RData")

AICTost=AIC(Tost.hhnx, Tost.hhn.detTopo10x, Tost.hhn.detWaterx, Tost.hhn.DHabx, Tost.hhn.DHab.nonUx, 
            Tost.hhn.DHab.nonU.Topo10x,Tost.hhn.DHab.nonU.Wx, Tost.hhn.DHab.nonU.T01Wx, 
            Tost.hhn.D.nonUx, Tost.hhn.DHab.nonU.GBx, Tost.hhn.DHab.nonU.GBGCx)

AICTost
coefficients(Tost.hhn.DHab.nonU.GBx)

load("./Noyon2013/Noyon-nonEuc-fitsx.RData") #Final round analysis (with water and topography)
AICNoyon=AIC(Noyon.hhnx, Noyon.hhn.detrgdx, Noyon.hhn.DHabx, Noyon.hhn.DHab.DetRgd01x, Noyon.hhn.detWaterx, 
             Noyon.hhn.detTopo10x, Noyon.hhn.detTopoWaterx, Noyon.hhn.DHab.nonUx, Noyon.hhn.D.nonUx, 
             Noyon.hhn.DHab.Topo10.nonUx, Noyon.hhn.DHab.DetW.nonUx, Noyon.hhn.DHab.Topo10W.nonUx,
             Noyon.hhn.DHab.nonU.GBx, Noyon.hhn.DHab.nonU.GBGCx, Noyon.hhn.DHab.DetRgd10x)
AICNoyon
coefficients(Noyon.hhn.DHab.nonU.GBx)

load("./Nemegt/Nemegt-nonEuc-fit2xR.RData")
NemegtAIC2xR=AIC(Nemegt.hhn2xR, Nemegt.hhn.detrgd2xR, Nemegt.hhn.DHab2xR, Nemegt.hhn.DHab.detrgd102xR, 
                 Nemegt.hhn.DHab.detrgd012xR, Nemegt.hhn.DHab.nonU2xR, Nemegt.hhn.D.nonU2xR,  
                 Nemegt.hhn.DHab.nonU.GB2xR, Nemegt.hhn.DHab.nonU.LamTopoR, Nemegt.hhn.DHab.nonU.LamW2xR,
                 Nemegt.hhn.DHab.nonU.LamTopoW2R)
NemegtAIC2xR
coefficients(Nemegt.hhn.DHab.nonU.LamW2xR)
#This is only where non Euc goes marginally negative!

load("./Tost_Noyon_Nemegt/TNN-NonEuc-fits2xR.RData")

TNNAIC2xR<- AIC(TNN.hhn.DRgd.sessR, TNN.hhn.DRgd.sess.DetWR, TNN.hhn.DRgd.DetTopo10WR, TNN.hhn.DHab.nonUR, 
                TNN.hhn.DHab.DetTopo10.nonUR, TNN.hhn.DHab.DetToposess.nonUR, 
                TNN.hhn.DHab.S.DetTopo10.nonUR, TNN.hhn.DHab_S.DetTopo10.nonUR, 
                TNN.hhn.DHab.LamTopoWat.nonUR, TNN.hhn.nonUR,TNN.hhn.DGC.DetTopo10.nonUGBR,
                TNN.hhn.DGB.DetTopo10.nonUR)
TNNAIC2xR
coefficients(TNN.hhn.DGC.DetTopo10.nonUGBR)

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


