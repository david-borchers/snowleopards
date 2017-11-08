
quartz(h=5,w=10)
plotcovariate(NemegtMask, covariate="GC", contour = FALSE, col = terrain.colors(16), zlim=zlim,asp=1)
NemegtALL_ch<-read.capthist(captfile = "./Analysis4paper/Data/Nemegt_capthist.csv", 
                         trapfile = "./Analysis4paper/Data/Nemegt_ALLcams.csv", 
                         detector="count", binary.usage=FALSE, fmt = "trapID", 
                         trapcovnames = c("Rgd", "Topo", "Water", "Winter"))
plot(traps(NemegtALL_ch),add=TRUE)
plot(NemegtALL_ch,tracks=TRUE,add=TRUE)
newtraps = locator(17)
newtraps = data.frame(
  x=c(730363.1,735499.8,719645.6,712332.4,657476.3,630080.4,621393.8,616609.1,606037.1,612461.0,
      607762.4,623464.7,615800.3,613404.8,685895.0, 727806.0, 715011.3),
  y=c(4834546,4829676,4827372,4832861,4829669,4835907,4843254,4837188,4831607,4829172,4824698,
      4827227, 4822124,4815395,4833884, 4828713, 4820258)
    )
points(newtraps,pch="+")
newtraps$ID = paste("NemegtSLTADD",1:dim(newtraps)[1],sep="")
NEW.traps = read.traps(data=newtraps,detector="proximity",trapID="ID")
NEW.traps = addCovariates(NEW.traps,NemegtMask,columns="GRIDCODE")
names(covariates(NEW.traps)) = "Rgd"
NEW.traps.df = data.frame(ID=row.names(NEW.traps),X=NEW.traps$x,Y=NEW.traps$y,Rgd=covariates(NEW.traps)$Rgd)
write.csv(NEW.traps.df,file="./Analysis4paper/Data/Nemegt_NEW_traps.csv",row.names=FALSE)
