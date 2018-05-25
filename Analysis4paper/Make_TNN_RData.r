library(secr)
library(fields)
library(maptools)
source("scrplotting.r")

# This file makes the .RData objects that are used by TNN.r
# =========================================================

# -------------------- Start Capture Histories --------------------------
# Combined regions
TNN.trapfiles = c(
  "./Analysis4paper/Data/Tost_cams.csv",
  "./Analysis4paper/Data/Noyon_cams.csv",
  "./Analysis4paper/Data/Nemegt_cams.csv"
)
TNN_ch<-read.capthist(captfile = "./Analysis4paper/Data/TNN_capthist.csv", 
                      binary.usage = FALSE, trapfile = TNN.trapfiles, 
                      detector="count", fmt = "trapID", 
                      trapcovnames = c("Rgd","Topo", "Water", "Winter"))

all.data.Tost<-read.capthist(captfile = "./Tost/Tost_capthist2012.csv", binary.usage = FALSE,  
                             trapfile = "./Tost/Tost_cams_rugged2012.csv", detector="count", 
                             fmt = "trapID", trapcovnames = c("Rgd", "Topo",	"Altidute",	"Water"))

summary(TNN_ch)
summary(covariates(traps(TNN_ch)))
names(covariates(traps(TNN_ch)))

# Tost and Noyon
TN.trapfiles = c(
  "./Analysis4paper/Data/Tost_cams.csv",
  "./Analysis4paper/Data/Noyon_cams.csv"
)
TN_ch<-read.capthist(captfile = "./Analysis4paper/Data/TN_capthist.csv", 
                      binary.usage = FALSE, trapfile = TN.trapfiles, 
                      detector="count", fmt = "trapID", 
                      trapcovnames = c("Rgd","Topo", "Water", "Winter"))
summary(TN_ch)

# Individual regions
Tost_ch<-read.capthist(captfile = "./Analysis4paper/Data/Tost_capthist.csv",   
                       trapfile = TNN.trapfiles[1], 
                       detector="count", binary.usage = FALSE, fmt = "trapID", 
                       trapcovnames = c("Rgd", "Topo", "Water", "Winter"))
Noyon_ch<-read.capthist(captfile = "./Analysis4paper/Data/Noyon_capthist.csv", 
                        trapfile = TNN.trapfiles[2], 
                        detector="count", binary.usage=FALSE, fmt = "trapID", 
                        trapcovnames = c("Rgd", "Topo", "Water", "Winter"))
Nemegt_ch<-read.capthist(captfile = "./Analysis4paper/Data/Nemegt_capthist.csv", 
                         trapfile = TNN.trapfiles[3], 
                         detector="count", binary.usage=FALSE, fmt = "trapID", 
                         trapcovnames = c("Rgd", "Topo", "Water", "Winter"))

# # Standarize Rgd on traps across all regions
# # ====
# # # NOTE: This is now redundant - replaced by code for covariates after dealing with mask covariates
# ====
# # ------------------------------------------
# TNN.cams=traps(TNN_ch)
# # look at covariates for each session:
# lapply(covariates(TNN.cams),summary)
# # To standardise in same way over all sessions, need to combine, standarise and then separate:
# Rgds = c(covariates(TNN.cams[[1]])$Rgd,
#          covariates(TNN.cams[[2]])$Rgd,
#          covariates(TNN.cams[[3]])$Rgd)
# stdRgd = scale(Rgds)
# n1 = dim(TNN.cams[[1]])[1]
# n2 = dim(TNN.cams[[2]])[1]
# n3 = dim(TNN.cams[[3]])[1]
# covariates(TNN.cams[[1]])$stdRgd = stdRgd[1:n1]
# covariates(TNN.cams[[2]])$stdRgd = stdRgd[(n1+1):(n1+n2)]
# covariates(TNN.cams[[3]])$stdRgd = stdRgd[(n1+n2+1):(n1+n2+n3)]
# 
# # Put trap covariates stdRgd into individual capture history files
# covariates(traps(Tost_ch))$stdRgd = covariates(TNN.cams[[1]])$stdRgd
# covariates(traps(Noyon_ch))$stdRgd = covariates(TNN.cams[[2]])$stdRgd
# covariates(traps(Nemegt_ch))$stdRgd = covariates(TNN.cams[[3]])$stdRgd
# 
# # look at covariates for each session again:
# lapply(covariates(TNN.cams),summary)
# 
# # put traps back into capthist (if don't put back by list elements, all.data.TNN becomes a 'traps' object!)
# traps(TNN_ch[[1]]) = TNN.cams[[1]]
# traps(TNN_ch[[2]]) = TNN.cams[[2]]
# traps(TNN_ch[[3]]) = TNN.cams[[3]]
# 
# summary(covariates(traps(TNN_ch))[[1]])
# summary(covariates(traps(TNN_ch))[[2]])
# summary(covariates(traps(TNN_ch))[[3]])
# 
# 
# # Standarize Rgd on traps across Tost and Noyon
# # ====
# # NOTE: This is now redundant - replaced by code for covariates after dealing with mask covariates
# # ====
# # ---------------------------------------------
# TN.cams=traps(TN_ch)
# # look at covariates for each session:
# lapply(covariates(TN.cams),summary)
# # To standardise in same way over all sessions, need to combine, standarise and then separate:
# Rgds = c(covariates(TN.cams[[1]])$Rgd,
#          covariates(TN.cams[[2]])$Rgd)
# stdRgd = scale(Rgds)
# n1 = dim(TN.cams[[1]])[1]
# n2 = dim(TN.cams[[2]])[1]
# covariates(TN.cams[[1]])$stdRgd = stdRgd[1:n1]
# covariates(TN.cams[[2]])$stdRgd = stdRgd[(n1+1):(n1+n2)]
# 
# # Put trap covariates stdRgd into individual capture history files
# covariates(traps(Tost_ch))$stdRgd = covariates(TN.cams[[1]])$stdRgd
# covariates(traps(Noyon_ch))$stdRgd = covariates(TN.cams[[2]])$stdRgd
# 
# # look at covariates for each session again:
# lapply(covariates(TN.cams),summary)
# 
# # put traps back into capthist (if don't put back by list elements, all.data.TN becomes a 'traps' object!)
# traps(TN_ch[[1]]) = TN.cams[[1]]
# traps(TN_ch[[2]]) = TN.cams[[2]]
# 
# summary(covariates(traps(TN_ch))[[1]])
# summary(covariates(traps(TN_ch))[[2]])
# # -------------------- End Capture Histories --------------------------


# -------------------- Start Boundary Files  --------------------------
boundaryNemegt=readShapeSpatial("./Nemegt//Habitat/Nemegt_StudyArea.shp")
boundaryNoyon=readShapeSpatial("./Noyon2013/Habitat/NoyonStudy_Area.shp")
boundaryTost=readShapeSpatial("./Tost//Habitat/TostStudy_Area.shp")
# -------------------- End Boundary Files  --------------------------


# -------------------- Start Masks --------------------------
NemegtMask=make.mask(traps(Nemegt_ch), spacing=500, buffer = 25000, type="trapbuffer", poly=boundaryNemegt)
NoyonMask=make.mask(traps(Noyon_ch), spacing=500, buffer = 25000, type="trapbuffer", poly=boundaryNoyon)
TostMask=make.mask(traps(Tost_ch), spacing=500, buffer = 25000, type="trapbuffer", poly=boundaryTost)

# Read ruggedness covariate and put it into mask covariate GRIDCODE
SLCost.Tost<-readShapePoly("./Analysis4paper/Data/Habitat/tost_rgd500m.shp")  #ruggedness pixels averaged over 500m radius
SLCost.Noyon<-readShapePoly("./Analysis4paper/Data/Habitat/noyon_rgd500m.shp")  #ruggedness pixels averaged over 500m radius
SLCost.Nemegt<-readShapePoly("./Analysis4paper/Data/Habitat/nemegt_rgd500m.shp")  #ruggedness pixels averaged over 500m radius

# Read logistic binary SL habitat created using telemetry data
SLCostBINARY<-readShapePoly("./Analysis4paper/Data/Habitat/tost_sl.shp")  

# Add covariates to masks
TostMask<-addCovariates(TostMask, SLCost.Tost)
NoyonMask<-addCovariates(NoyonMask, SLCost.Noyon)
NemegtMask<-addCovariates(NemegtMask, SLCost.Nemegt)
names(covariates(TostMask))
plot(TostMask,covariate="GRIDCODE")
names(covariates(NoyonMask))
plot(NoyonMask,covariate="GRIDCODE")
names(covariates(NemegtMask))
plot(NemegtMask,covariate="GRIDCODE")

TostMask<-addCovariates(TostMask, SLCostBINARY)
names(covariates(TostMask))[3:4] = c("binaryID","BINCODE") #Rename headers
covariates(TostMask)$BINCODE[is.na(covariates(TostMask)$BINCODE)] = 0 # make NAs in BINCODE zeros
summary(covariates(TostMask))
NoyonMask<-addCovariates(NoyonMask, SLCostBINARY)
names(covariates(NoyonMask))[3:4] = c("binaryID","BINCODE") #Rename headers
covariates(NoyonMask)$BINCODE[is.na(covariates(NoyonMask)$BINCODE)] = 0 # make NAs in BINCODE zeros
summary(covariates(NoyonMask))
NemegtMask<-addCovariates(NemegtMask, SLCostBINARY)
names(covariates(NemegtMask))[3:4] = c("binaryID","BINCODE") #Rename headers
covariates(NemegtMask)$BINCODE[is.na(covariates(NemegtMask)$BINCODE)] = 0 # make NAs in BINCODE zeros
summary(covariates(NemegtMask))

# Standarize GRIDCODE (in stdGC) ACROSS ALL REGIONS
# (BINCODE could be factor but then getting region mean is awkward.)
# ------------------------------------------------------------------------
NemGC = covariates(NemegtMask)$GRIDCODE
NemBC = covariates(NemegtMask)$BINCODE
NoyGC = covariates(NoyonMask)$GRIDCODE
NoyBC = covariates(NoyonMask)$BINCODE
TosGC = covariates(TostMask)$GRIDCODE
TosBC = covariates(TostMask)$BINCODE
nNem = length(NemGC)
nNoy = length(NoyGC)
nTos = length(TosGC)
# put all GC together 
GC = c(TosGC,NoyGC,NemGC)
# mean and var of GC; to use later
meanGC = mean(GC)
sdGC = sqrt(var(GC))
# Standardize GC
stdGC = scale(GC)
# put all BC together 
BC = c(TosGC,NoyGC,NemGC)
# mean and var of BC; to use later
meanBC = mean(BC)
sdBC = sqrt(var(BC))
# Standardize BC
stdBC = scale(BC)

# Create region, mean GC and BC, and total GC and BC variables
in.Tos = 1:nTos
in.Noy = (nTos+1):(nTos+nNoy)
in.Nem = (nTos+nNoy+1):(nTos+nNoy+nNem)
#region = c(rep("Nemegt",nNem),rep("Noyon",nNoy),rep("Tost",nTos))
rsumGC = rsumBC = rmeanGCdev = rmeanBCdev = rmeanGC = rmeanBC = rep(NA,(nTos+nNoy+nNem))
# create region-specific mean covariates
rmeanGC[in.Tos] = mean(stdGC[in.Tos])
rmeanBC[in.Tos] = mean(stdBC[in.Tos])
rmeanGC[in.Noy] = mean(stdGC[in.Noy])
rmeanBC[in.Noy] = mean(stdBC[in.Noy])
rmeanGC[in.Nem] = mean(stdGC[in.Nem])
rmeanBC[in.Nem] = mean(stdBC[in.Nem])
# create region-specific deviation-from-mean covariates
rmeanGCdev[in.Tos] = stdGC[in.Tos] - rmeanGC[in.Tos]
rmeanBCdev[in.Tos] = stdGC[in.Tos] - rmeanBC[in.Tos]
rmeanGCdev[in.Noy] = stdGC[in.Noy] - rmeanGC[in.Noy]
rmeanBCdev[in.Noy] = stdGC[in.Noy] - rmeanBC[in.Noy]
rmeanGCdev[in.Nem] = stdGC[in.Nem] - rmeanGC[in.Nem]
rmeanBCdev[in.Nem] = stdGC[in.Nem] - rmeanBC[in.Nem]
# create region-specific sum covariates
rsumGC[in.Tos] = sum(stdGC[in.Tos])
rsumBC[in.Tos] = sum(stdBC[in.Tos])
rsumGC[in.Noy] = sum(stdGC[in.Noy])
rsumBC[in.Noy] = sum(stdBC[in.Noy])
rsumGC[in.Nem] = sum(stdGC[in.Nem])
rsumBC[in.Nem] = sum(stdBC[in.Nem])
# scale the three sums:
scaledsumGC = scale(c(rsumGC[in.Tos][1],rsumGC[in.Noy][1],rsumGC[in.Nem][1]))
scaledsumBC = scale(c(rsumBC[in.Tos][1],rsumBC[in.Noy][1],rsumBC[in.Nem][1]))
rsumGC[in.Tos]=scaledsumBC[1]
rsumGC[in.Noy]=scaledsumBC[2]
rsumGC[in.Nem]=scaledsumBC[3]
rsumBC[in.Tos]=scaledsumBC[1]
rsumBC[in.Noy]=scaledsumBC[2]
rsumBC[in.Nem]=scaledsumBC[3]
# Put variables back in masks
covariates(NemegtMask)$stdGC = stdGC[in.Nem]
covariates(NemegtMask)$stdBC = stdBC[in.Nem]
covariates(NemegtMask)$GC = GC[in.Nem]
covariates(NemegtMask)$BC = BC[in.Nem]
covariates(NemegtMask)$rmeanGC = rmeanGC[in.Nem]
covariates(NemegtMask)$rmeanBC = rmeanBC[in.Nem]
covariates(NemegtMask)$rmeanGCdev = rmeanGCdev[in.Nem]
covariates(NemegtMask)$rmeanBCdev = rmeanBCdev[in.Nem]
covariates(NemegtMask)$rsumGC = rsumGC[in.Nem]
covariates(NemegtMask)$rsumBC = rsumBC[in.Nem]

covariates(NoyonMask)$stdGC = stdGC[in.Noy]
covariates(NoyonMask)$stdBC = stdBC[in.Noy]
covariates(NoyonMask)$GC = GC[in.Noy]
covariates(NoyonMask)$BC = BC[in.Noy]
covariates(NoyonMask)$rmeanGC = rmeanGC[in.Noy]
covariates(NoyonMask)$rmeanBC = rmeanBC[in.Noy]
covariates(NoyonMask)$rmeanGCdev = rmeanGCdev[in.Noy]
covariates(NoyonMask)$rmeanBCdev = rmeanBCdev[in.Noy]
covariates(NoyonMask)$rsumGC = rsumGC[in.Noy]
covariates(NoyonMask)$rsumBC = rsumBC[in.Noy]

covariates(TostMask)$stdGC = stdGC[in.Tos]
covariates(TostMask)$stdBC = stdBC[in.Tos]
covariates(TostMask)$GC = GC[in.Tos]
covariates(TostMask)$BC = BC[in.Tos]
covariates(TostMask)$rmeanGC = rmeanGC[in.Tos]
covariates(TostMask)$rmeanBC = rmeanBC[in.Tos]
covariates(TostMask)$rmeanGCdev = rmeanGCdev[in.Tos]
covariates(TostMask)$rmeanBCdev = rmeanBCdev[in.Tos]
covariates(TostMask)$rsumGC = rsumGC[in.Tos]
covariates(TostMask)$rsumBC = rsumBC[in.Tos]

# -------------------- End Masks --------------------------
TNN.cams=traps(TNN_ch)
# look at covariates for each session:
lapply(covariates(TNN.cams),summary)

# Put mask covariates into individual trap files
covnames = c("GRIDCODE","BINCODE","stdGC","stdBC","GC","BC","rmeanGC","rmeanBC","rmeanGCdev","rsumGC","rsumBC")
TNN.cams[[1]] = addCovariates(TNN.cams[[1]],TostMask,columns=covnames,replace=TRUE)
TNN.cams[[2]] = addCovariates(TNN.cams[[2]],NoyonMask,columns=covnames,replace=TRUE)
TNN.cams[[3]] = addCovariates(TNN.cams[[3]],NemegtMask,columns=covnames,replace=TRUE)

# look at covariates for each session again:
lapply(covariates(TNN.cams),summary)

# put traps back into capthist (if don't put back by list elements, all.data.TNN becomes a 'traps' object!)
traps(TNN_ch[[1]]) = TNN.cams[[1]]
traps(TNN_ch[[2]]) = TNN.cams[[2]]
traps(TNN_ch[[3]]) = TNN.cams[[3]]
Tost_ch = TNN_ch[[1]]
Noyon_ch = TNN_ch[[2]]
Nemegt_ch = TNN_ch[[3]]
# -------------------- Put all mask covariates into traps --------------------------


# ----------------- END of put all mask covariates into traps ----------------------

# -------------------- Start RData Save --------------------------
Tostboundary = boundaryTost
Noyonboundary = boundaryNoyon
Nemegtboundary = boundaryNemegt
save(Tostboundary,Noyonboundary,Nemegtboundary,file="./Analysis4paper/TNN_boundaries.RData")
save(TostMask,NoyonMask,NemegtMask,file="./Analysis4paper/TNN_masks.RData")
save(Tost_ch,Noyon_ch,Nemegt_ch,TNN_ch,file="./Analysis4paper/TNN_caphists.RData")
#save(Tost_ch,Noyon_ch,TN_ch,file="./Analysis4paper/TN_caphists.RData")
# -------------------- End RData Save --------------------------
