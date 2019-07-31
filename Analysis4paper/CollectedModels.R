# Looking at all models, including Ian's non-Euc models

# Models from Ian fitting:
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


# Look at AICc for new and old models together:
AIC(TNNfit.DGrid.a0Waterxsess_topo.sig0,TNNfit.DGrid.a0Waterxsess_topo.sig1,
    TNNfit.DGrid.a0Waterxsess_topo.sig2,TNNfit.DGrid.a0Waterxsess_topo.sig3,
    TNN.GCmean_dev.NE,TNN.GCmean_dev.WW,
    TNN.GCmean_dev.WW.session.sigma,
    TNN.GCmean_dev.WW.session,
    TNNfit.DGridXSess.a0Sess.sigsess, TNNfit.DGrid_Sess.a0Sess.sigsess, TNNfit.DGrid_Sess.a0TopoXSess.sigsess, TNNfit.DGrid_Sess.a0Topo.Sess.sigsess, 
    TNNfit.DGrid.a0Topo.sig1, TNNfit.DGrid.a0TopoSess.sig1, TNNfit.DGrid.a0TopoSess.sigsess, TNNfit.DGridSess.lamSess.SigSess, TNN.DGrid.LamTopo.Sig1, 
    TNN.DGridXsess.LamTopoXsess, TNN.DGridXsess.DetTopoXsess, TNNfit.D1.LamSess.sig1, TNNfit.D1.Lam1s.sig1, TNNfit.Dsess.Lamsess.sigsess, 
    TNNfit.DSess.a0TopoSess.sigSess, TNNfit.DGrid_Sess.a0Topo.sig1, TNNfit.DGridxSess.a0Topo.sig1, TNNfit.DGrid.a0Topo_sess.sig1)
