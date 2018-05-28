# ================== MODELS WITHOUT NONEUC ===================================

library(secr)
library(fields)
library(maptools)
source("scrplotting.r")
library(parallel)

detectCores()
setwd("E:/public/users/koustubh")

#----------------------- Load the RData objects made by Make_TNN_RData.r --------------------------
load("./Analysis4paper/TNN_boundaries.RData") # Tostboundary,Noyonboundary,Nemegtboundary
load("./Analysis4paper/TNN_masks.RData") # TostMask,NoyonMask,NemegtMask
#load("./Analysis4paper/TN_caphists.RData") # Tost_ch,Noyon_ch,TN_ch
load("./Analysis4paper/TNN_caphists.RData") # Tost_ch,Noyon_ch,Nemegt_sh,TNN_ch)
#----------------------- ----------------------------------------------- --------------------------
summary(TNN_ch)

# ===================== Tost, Noyon and Nemegt ====================

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

#8a Density varies with ruggedness, is constan across 3 areas. a0 is a function of Topography OR session, and sigma is constant
TNNfit.DGrid.a0Topo_sess.sig1<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                   model=list(D~stdGC, a0~Topo+session, sigma~1), ncores=16)

#8b Density varies with ruggedness OR session, is constan across 3 areas. a0 is a function of Topography , and sigma is constant
TNNfit.DGridxSess.a0Topo.sig1<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                        model=list(D~stdGC*session, a0~Topo, sigma~1), ncores=16)

#8c Density varies with ruggedness OR session, is constan across 3 areas. a0 is a function of Topography , and sigma is constant
TNNfit.DGrid_Sess.a0Topo.sig1<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                        model=list(D~stdGC+session, a0~Topo, sigma~1), ncores=16)

# 8d Density varies with ruggedness, constant for three session. a0 varies with water at trap site
TNNfit.DGrid.a0Water.sig1<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                   model=list(D~stdGC, a0~Water, sigma~1), ncores=16)

# 8e Density varies with ruggedness, constant for three session. a0 varies with water at trap site and session
TNNfit.DGrid.a0Waterxsess.sig1<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                    model=list(D~stdGC, a0~Water*session, sigma~1), ncores=16)

# 8f Density varies with ruggedness, constant for three session. a0 varies with water at trap site and session
TNNfit.DGrid.a0Water_topo.sig1<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                         model=list(D~stdGC, a0~Water+Topo, sigma~1), ncores=16)

# 8g Density varies with ruggedness, constant for three session. a0 varies with water at trap site and session
TNNfit.DGrid.a0Waterxsess_topo.sig1<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                         model=list(D~stdGC, a0~Topo+Water*session, sigma~1), ncores=16)

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
TNNfit.D1.a0TopoSess.sig1<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                    model=list(D~1, a0~Topo*session, sigma~1), ncores=16)

#13 Density a function of ruggedness and varies differently in 3 areas. a0 varies with topography differently for 3 areas. Sigma is constant
TNNfit.DGridSess.a0TopoSess.sig1<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                           model=list(D~rmeanGC + rmeanGCdev: session + session, a0~Topo*session, sigma~1), ncores=16)

#14 Density a function of ruggedness and varies differently between 3 areas. a0 varies with topography differently for 3 areas.
#   Sigma different for each area
TNNfit.DGridSess.a0TopoSess.sigsess<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                              model=list(D~rmeanGC + rmeanGCdev: session + session, a0~Topo*session, sigma~session), ncores=16)

############################################################################
#12 listed
TNNfit.D1.a0TopoSess.sig1000<-list(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                    model=list(D~1, a0~Topo*session, sigma~1), ncores=3)

#13 listed
TNNfit.DGridSess.a0TopoSess.sig1000<-list(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                           model=list(D~rmeanGC + rmeanGCdev: session + session, a0~Topo*session, sigma~1), ncores=3)

#14 listed
TNNfit.DGridSess.a0TopoSess.sigsess000<-list(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                              model=list(D~rmeanGC + rmeanGCdev: session + session, a0~Topo*session, sigma~session), ncores=3)

fits12_13_14_000<-par.secr.fit (c('TNNfit.D1.a0TopoSess.sig1000', 'TNNfit.DGridSess.a0TopoSess.sig1000', 'TNNfit.DGridSess.a0TopoSess.sigsess000'), ncores = 32)

coefficients(fits12_13_14_000)
TNNfit.D1.a0TopoSess.sig1<-fits12_13_14_000[[1]]
TNNfit.DGridSess.a0TopoSess.sig1<-fits12_13_14_000[[2]]
TNNfit.DGridSess.a0TopoSess.sigsess<-fits12_13_14_000[[3]]


############################################################################

#15 Density a function of ruggedness and varies differently between 3 areas. a0 different for three areas, sigma different for 3 areas.
TNNfit.DGridSess.a0Sess.sigsess<-secr.fit(TNN_ch,detectfn="HHN", list(TostMask, NoyonMask, NemegtMask),
                                          model=list(D~rmeanGC + rmeanGCdev: session + session, a0~session, sigma~session), ncores=16)
TNNfit.DGridSess.a0Sess.sigsess

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

##### Ran 16, 17 18, 19 and 8, 9, 10 because these three were not using stdGC but a summation model! Also 12, 13 14.
##### Will need to merge all model outputs as a single rdata file once complete. Will have to ignore model runs for 8, 9, 10 on thinkpad
##### 

save(TNNfit.DGridSess.lamSess.SigSess, TNN.DGrid.LamTopo.Sig1, TNN.DGridXsess.LamTopoXsess, TNN.DGridXsess.DetTopoXsess, TNNfit.D1.LamSess.sig1, 
    TNNfit.D1.Lam1s.sig1, TNNfit.Dsess.Lamsess.sigsess, TNNfit.DGrid.a0Topo.sig1, TNNfit.DGrid.a0Topo_sess.sig1, TNNfit.DGridxSess.a0Topo.sig1, 
    TNNfit.DGrid_Sess.a0Topo.sig1, TNNfit.DGrid.a0Water.sig1, TNNfit.DGrid.a0Waterxsess.sig1, TNNfit.DGrid.a0Water_topo.sig1, 
    TNNfit.DGrid.a0Waterxsess_topo.sig1, TNNfit.DGrid.a0TopoSess.sig1, TNNfit.DGrid.a0TopoSess.sigsess, TNNfit.DSess.a0TopoSess.sigSess, 
    TNNfit.D1.a0TopoSess.sig1, TNNfit.DGridSess.a0TopoSess.sig1, TNNfit.DGridSess.a0TopoSess.sigsess, TNNfit.DGridSess.a0Sess.sigsess, 
    TNNfit.DGridXSess.a0Sess.sigsess, TNNfit.DGrid_Sess.a0Sess.sigsess, TNNfit.DGrid_Sess.a0TopoXSess.sigsess, TNNfit.DGrid_Sess.a0Topo.Sess.sigsess,
    file = "E:/public/users/koustubh/Analysis4paper/TNNnoNEfits.RData")

AllAIC<-AIC(TNNfit.DGridSess.lamSess.SigSess, TNN.DGrid.LamTopo.Sig1, TNN.DGridXsess.LamTopoXsess, TNN.DGridXsess.DetTopoXsess, TNNfit.D1.LamSess.sig1, 
     TNNfit.D1.Lam1s.sig1, TNNfit.Dsess.Lamsess.sigsess, TNNfit.DGrid.a0Topo.sig1, TNNfit.DGrid.a0Topo_sess.sig1, TNNfit.DGridxSess.a0Topo.sig1, 
     TNNfit.DGrid_Sess.a0Topo.sig1, TNNfit.DGrid.a0Water.sig1, TNNfit.DGrid.a0Waterxsess.sig1, TNNfit.DGrid.a0Water_topo.sig1, 
     TNNfit.DGrid.a0Waterxsess_topo.sig1, TNNfit.DGrid.a0TopoSess.sig1, TNNfit.DGrid.a0TopoSess.sigsess, TNNfit.DSess.a0TopoSess.sigSess, 
     TNNfit.D1.a0TopoSess.sig1, TNNfit.DGridSess.a0TopoSess.sig1, TNNfit.DGridSess.a0TopoSess.sigsess, TNNfit.DGridSess.a0Sess.sigsess, 
     TNNfit.DGridXSess.a0Sess.sigsess, TNNfit.DGrid_Sess.a0Sess.sigsess, TNNfit.DGrid_Sess.a0TopoXSess.sigsess, TNNfit.DGrid_Sess.a0Topo.Sess.sigsess)

save(TNNfit.DGridXSess.a0Sess.sigsess, TNNfit.DGrid_Sess.a0Sess.sigsess, TNNfit.DGrid_Sess.a0TopoXSess.sigsess, TNNfit.DGrid_Sess.a0Topo.Sess.sigsess, 
     TNNfit.DGrid.a0Topo.sig1, TNNfit.DGrid.a0TopoSess.sig1, TNNfit.DGrid.a0TopoSess.sigsess, TNNfit.DGridSess.lamSess.SigSess, TNN.DGrid.LamTopo.Sig1, 
     TNN.DGridXsess.LamTopoXsess, TNN.DGridXsess.DetTopoXsess, TNNfit.D1.LamSess.sig1, TNNfit.D1.Lam1s.sig1, TNNfit.Dsess.Lamsess.sigsess, 
     TNNfit.DSess.a0TopoSess.sigSess, TNNfit.DGrid_Sess.a0Topo.sig1, TNNfit.DGridxSess.a0Topo.sig1, TNNfit.DGrid.a0Topo_sess.sig1, 
     file= "./Analysis4paper/TNNnoNEfits_all.RData")

getwd()
# Load models run in thinkpad
load("./Analysis4paper/TNNnoNEfits_Thinkpad.RData")
load("E:/public/users/koustubh/Analysis4paper/TNNnoNEfits_all.RData")


TNN_AICc<-AIC(TNNfit.DGridXSess.a0Sess.sigsess, TNNfit.DGrid_Sess.a0Sess.sigsess, TNNfit.DGrid_Sess.a0TopoXSess.sigsess, TNNfit.DGrid_Sess.a0Topo.Sess.sigsess, 
    TNNfit.DGrid.a0Topo.sig1, TNNfit.DGrid.a0TopoSess.sig1, TNNfit.DGrid.a0TopoSess.sigsess, TNNfit.DGridSess.lamSess.SigSess, TNN.DGrid.LamTopo.Sig1, 
    TNN.DGridXsess.LamTopoXsess, TNN.DGridXsess.DetTopoXsess, TNNfit.D1.LamSess.sig1, TNNfit.D1.Lam1s.sig1, TNNfit.Dsess.Lamsess.sigsess, 
    TNNfit.DSess.a0TopoSess.sigSess, TNNfit.DGrid_Sess.a0Topo.sig1, TNNfit.DGridxSess.a0Topo.sig1, TNNfit.DGrid.a0Topo_sess.sig1)

write.csv(TNN_AICc, "./Analysis4paper/TNN_AICC.csv")

coefficients(TNNfit.DGrid.a0Topo.sig1)
