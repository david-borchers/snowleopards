# Looking at all models, including Ian's non-Euc models

# Models from Ian fitting (plus the lambda0 version of best model):
load("Analysis4paper/nonEuc-mdls.RData")
AIC(TNNfit.DGrid.a0Waterxsess_topo.sig0,TNNfit.DGrid.a0Waterxsess_topo.sig1,
    TNNfit.DGrid.a0Waterxsess_topo.sig2,TNNfit.DGrid.a0Waterxsess_topo.sig3,
    TNN.GCmean_dev.NE,TNN.GCmean_dev.WW)


# Old models (these use lambda0, not a0)
load("./Analysis4paper/TNNFits1.RData")
AIC(TNN.GCmean_dev.WW.session.sigma,
     TNN.GCmean_dev.WW.session,
     TNN.GCmean_dev.WW)

#  Old models created by Koustubh (*_Thinkpad.RData) and David (*_all.RData) (some of these use lambda0, not a0)
load("./Analysis4paper/TNNnoNEfits_Thinkpad.RData")
load("./Analysis4paper/TNNnoNEfits_all.RData")
TNN_AICc<-AIC(TNNfit.DGridXSess.a0Sess.sigsess, TNNfit.DGrid_Sess.a0Sess.sigsess, TNNfit.DGrid_Sess.a0TopoXSess.sigsess, TNNfit.DGrid_Sess.a0Topo.Sess.sigsess, 
              TNNfit.DGrid.a0Topo.sig1, TNNfit.DGrid.a0TopoSess.sig1, TNNfit.DGrid.a0TopoSess.sigsess, TNNfit.DGridSess.lamSess.SigSess, TNN.DGrid.LamTopo.Sig1, 
              TNN.DGridXsess.LamTopoXsess, TNN.DGridXsess.DetTopoXsess, TNNfit.D1.LamSess.sig1, TNNfit.D1.Lam1s.sig1, TNNfit.Dsess.Lamsess.sigsess, 
              TNNfit.DSess.a0TopoSess.sigSess, TNNfit.DGrid_Sess.a0Topo.sig1, TNNfit.DGridxSess.a0Topo.sig1, TNNfit.DGrid.a0Topo_sess.sig1)
#              TNNfit.D2Grid.a0Waterxsess_topo.sig1) # removed this because it did not exist (looks like might not have converged, with quadratic on density)


# Look at AICc for new and old models together, and check lambda0 and a0 for first 4 models below are same:
AIC(TNNfit.DGrid.a0Waterxsess_topo.sig0,TNNfit.DGrid.a0Waterxsess_topo.sig1,
    TNNfit.DGrid.a0Waterxsess_topo.sig2,TNNfit.DGrid.a0Waterxsess_topo.sig3,
    TNNfit.DGrid.lam0Waterxsess_topo.sig0,TNNfit.DGrid.lam0Waterxsess_topo.sig1, # these, and 
    TNNfit.DGrid.lam0Waterxsess_topo.sig2,TNNfit.DGrid.lam0Waterxsess_topo.sig3, # these, are same as above, but with lambda0, not a0
    TNN.GCmean_dev.NE,TNN.GCmean_dev.WW,
    TNN.GCmean_dev.WW.session.sigma,
    TNN.GCmean_dev.WW.session,
    TNNfit.DGridXSess.a0Sess.sigsess, TNNfit.DGrid_Sess.a0Sess.sigsess, TNNfit.DGrid_Sess.a0TopoXSess.sigsess, TNNfit.DGrid_Sess.a0Topo.Sess.sigsess, 
    TNNfit.DGrid.a0Topo.sig1, TNNfit.DGrid.lam0Topo.sig1, # 5th best model
    TNNfit.DGrid.a0TopoSess.sig1, TNNfit.DGrid.a0TopoSess.sigsess, TNNfit.DGridSess.lamSess.SigSess, TNN.DGrid.LamTopo.Sig1, 
    TNN.DGridXsess.LamTopoXsess, TNN.DGridXsess.DetTopoXsess, TNNfit.D1.LamSess.sig1, TNNfit.D1.Lam1s.sig1, TNNfit.Dsess.Lamsess.sigsess, 
    TNNfit.DSess.a0TopoSess.sigSess, TNNfit.DGrid_Sess.a0Topo.sig1, TNNfit.DGridxSess.a0Topo.sig1, TNNfit.DGrid.a0Topo_sess.sig1)


# Look at AICc only for models without convergence problems: a0 instead of lambda0 parameterization for best models:
AIC(TNNfit.DGrid.a0Waterxsess_topo.sig0,TNNfit.DGrid.a0Waterxsess_topo.sig1,
    TNNfit.DGrid.a0Waterxsess_topo.sig2,TNNfit.DGrid.a0Waterxsess_topo.sig3,
    TNNfit.DGridXSess.a0Sess.sigsess, TNNfit.DGrid_Sess.a0Sess.sigsess, TNNfit.DGrid_Sess.a0TopoXSess.sigsess, TNNfit.DGrid_Sess.a0Topo.Sess.sigsess, 
    TNNfit.DGrid.a0Topo.sig1, TNNfit.DGrid.a0TopoSess.sig1, TNNfit.DGrid.a0TopoSess.sigsess, TNNfit.DGridSess.lamSess.SigSess, TNN.DGrid.LamTopo.Sig1, 
    TNN.DGridXsess.LamTopoXsess, TNN.DGridXsess.DetTopoXsess, TNNfit.D1.LamSess.sig1, TNNfit.D1.Lam1s.sig1, TNNfit.Dsess.Lamsess.sigsess, 
    TNNfit.DSess.a0TopoSess.sigSess, TNNfit.DGrid_Sess.a0Topo.sig1, TNNfit.DGridxSess.a0Topo.sig1, TNNfit.DGrid.a0Topo_sess.sig1)

# Look at AICc only for models without convergence problems: lambda0 parameterizations for best 5 models:
AIC(TNNfit.DGrid.lam0Waterxsess_topo.sig0,TNNfit.DGrid.lam0Waterxsess_topo.sig1, # these, and 
    TNNfit.DGrid.lam0Waterxsess_topo.sig2,TNNfit.DGrid.lam0Waterxsess_topo.sig3, # these, are the 4 bst models
    TNNfit.DGridXSess.a0Sess.sigsess, TNNfit.DGrid_Sess.a0Sess.sigsess, TNNfit.DGrid_Sess.a0TopoXSess.sigsess, TNNfit.DGrid_Sess.a0Topo.Sess.sigsess, 
    TNNfit.DGrid.lam0Topo.sig1, # 5th best model - replaced a0 with lambda0
    TNNfit.DGrid.a0TopoSess.sig1, TNNfit.DGrid.a0TopoSess.sigsess, TNNfit.DGridSess.lamSess.SigSess, TNN.DGrid.LamTopo.Sig1, 
    TNN.DGridXsess.LamTopoXsess, TNN.DGridXsess.DetTopoXsess, TNNfit.D1.LamSess.sig1, TNNfit.D1.Lam1s.sig1, TNNfit.Dsess.Lamsess.sigsess, 
    TNNfit.DSess.a0TopoSess.sigSess, TNNfit.DGrid_Sess.a0Topo.sig1, TNNfit.DGridxSess.a0Topo.sig1, TNNfit.DGrid.a0Topo_sess.sig1)

# get coefficients for best model, for tabulation in paper:
coefficients(TNNfit.DGrid.lam0Waterxsess_topo.sig3)

# get abundances for top 2 models, for use in paper:
rank1.N = region.N(TNNfit.DGrid.lam0Waterxsess_topo.sig3)
rank2.N = region.N(TNNfit.DGrid.lam0Waterxsess_topo.sig0)
Nhats = matrix(
  c(rank1.N[[1]][1,1:4],rank1.N[[2]][1,1:4],rank1.N[[3]][1,1:4],
    rank2.N[[1]][1,1:4],rank2.N[[2]][1,1:4],rank2.N[[3]][1,1:4]),
    byrow=TRUE,nrow=2)
row.names(Nhats) = c("Best","Second")
colnames(Nhats) = c(paste("Tost.",colnames(rank1.N[[1]][1:4]),sep=""),
                     paste("Noyon.",colnames(rank1.N[[1]][1:4]),sep=""),
                     paste("Nemegt.",colnames(rank1.N[[1]][1:4]),sep=""))
for(i in 1:2) for(j in 1:12) Nhats[i,j] = round(as.numeric(Nhats[i,j]),2)
Nhats

# Shit, realised N SEs are NA! WTF!
# Bootstrap instead:
library(mvtnorm)
library(tcltk2) # for progress bar

# best model
fit1 = TNNfit.DGrid.lam0Waterxsess_topo.sig3
vcv1 = fit1$beta.vcv
mvnmean1 = fit1$fit$par
Nhat1 = region.N(fit1)
# 2nd best model
fit2 = TNNfit.DGrid.lam0Waterxsess_topo.sig0
vcv2 = fit2$beta.vcv
mvnmean2 = fit2$fit$par
Nhat2 = region.N(fit2)

B = 1000 # number resamples
b.par1 = rmvnorm(n=B, mean=mvnmean1, sigma=vcv1) # get the resampled parameter estimates
b.par2 = rmvnorm(n=B, mean=mvnmean2, sigma=vcv2) # get the resampled parameter estimates
clone1 = fit1
clone2 = fit2
b.EN1 = b.RN1 = matrix(rep(NA,B*length(Nhat1)),nrow=B)
b.EN2 = b.RN2 = matrix(rep(NA,B*length(Nhat2)),nrow=B)
set.seed(1234)
# Progress bar setup
pb <- tkProgressBar(title=paste("Bootstrap Progress (B=",B,")",sep=""), min=0, max=B, width=400)
for(i in 1:B) {
  # best model
  clone1$fit$par = clone1$fit$estimate = b.par1[i,]
  clone1.Nhat = region.N(clone1)
  b.EN1[i,] = c(clone1.Nhat[[1]][1,1],clone1.Nhat[[2]][1,1],clone1.Nhat[[3]][1,1])
  b.RN1[i,] = c(clone1.Nhat[[1]][2,1],clone1.Nhat[[2]][2,1],clone1.Nhat[[3]][2,1])
  # 2nd best model
  clone2$fit$par = clone2$fit$estimate = b.par2[i,]
  clone2.Nhat = region.N(clone2)
  b.EN2[i,] = c(clone2.Nhat[[1]][1,1],clone2.Nhat[[2]][1,1],clone2.Nhat[[3]][1,1])
  b.RN2[i,] = c(clone2.Nhat[[1]][2,1],clone2.Nhat[[2]][2,1],clone2.Nhat[[3]][2,1])
  
  # update progress bar
  setTkProgressBar(pb, i, label=paste( round(i/B*100, 0),"% done"))
}
close(pb) # close progress bar

save(b.RN1,b.RN2,b.EN1,b.EN2,file="./Analysis4paper/TNNbest2fitsBootN.RData")

n = c(Nhat1[[1]]$n[1],Nhat1[[2]]$n[1],Nhat1[[3]]$n[1])
b.mean1 = apply(b.RN1,2,mean)
b.mean2 = apply(b.RN2,2,mean)
b.se1 = apply(b.RN1,2,sd)
b.se2 = apply(b.RN2,2,sd)
n
b.mean1;b.mean2
b.mean2/b.mean1
b.se1;b.se2

