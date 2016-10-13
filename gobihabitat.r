##
## NOTE: ...GIS Files/gobi/Boundaries/Habitat is not on Dropbox so also not on 
##       github. Should put it on at some point.
##

#Creating SL habitat vs non-habitat in Gobi
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}

RasterYN <- read.csv("C:/GIS Files/gobi/Boundaries/Habitat/RasterYN.csv", row.names=1)
(RasterYN$RASTERVALU1=="1")
length(RasterYN$RASTERVALU)
length(RasterYN$RASTERVALU[RasterYN$RASTERVALU=="3"])
RasterYN$RASTERVALU1<-as.numeric.factor(RasterYN$RASTERVALU)
RasterYN$RASTERVALU1[is.na(RasterYN$RASTERVALU1)]<-0
summary(RasterYN)
head(RasterYN)
catYN<-RasterYN[,c("POINTID", "GRID_CODE","RASTERVALU1")]
head(catYN)
max(catYN$RASTERVALU1)

catY<-catYN[rep(1:nrow(catYN), catYN$RASTERVALU1),1:3]
catY$Pres<-1

catN<-catYN[which(catYN$RASTERVALU1==0),1:3]
catN$Pres<-0

tail(catN)
head(catN)
nrow(catY)
names(catN)
names(catY)
catYN2<-rbind(catY,catN[sample(nrow(catN), 25000),])

#catYN2<-append(catY,catN,after=nrow(catY))
Rows1<-nrow(catY)+nrow(catN)
head(catYN2)
tail(catYN2)
nrow(catYN2)
hist(catYN2$Pres)
glm_full<-glm(Pres~GRID_CODE, family =binomial(logit), data=catYN2)
glm_25kSample<-glm(Pres~GRID_CODE, family =binomial(logit), data=catYN2)
summary(glm_full)
summary(glm_25kSample)
coefficients(glm_full)
coefficients(glm_25kSample)

glm.Area=glm(Extinct~Area_z,family=binomial(logit), data=extinct)

