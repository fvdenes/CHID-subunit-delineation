# First attemmpt to fit CAWA national model as a GW regression

##Step one
###Load dataset
library("GWmodel")
library(mefa4)
library(sp)
library(rgdal)
library(opticut)

load("C:/Users/voeroesd/Dropbox/BAM/Critical Habitat/CHID subunit delineation/pack_2016-12-01.Rdata")



DAT$HAB_NALC1 <- DAT$HABTR
DAT$HAB_NALC2 <- DAT$HAB
DAT$YEAR <- DAT$YR+2013

mm <- Mefa(YY, DAT, TAX, "inner")


# Combine CAWA counts, offsets and covariates in a single dataset
count <- YY[,"CAWA"]
offset <- OFF[,"CAWA"]


# thin_cawa_mm <-cawa_mm[sample(nrow(cawa_mm),2000),] # need to later figure out a proper sampling procedure (to make sure it represents spatial samples)



## resampling blocks
DAT$YR5 <- 0
DAT$YR5[DAT$YEAR > 2000] <- 1
DAT$YR5[DAT$YEAR > 2004] <- 2
DAT$YR5[DAT$YEAR > 2009] <- 3
table(DAT$YEAR,DAT$YR5)
table(DAT$YR5)
DAT$bootg <- interaction(DAT$Units, DAT$YR5, drop=TRUE)

bbfun2 <- function(DAT1, B, out=0.1, seed=1234) {
  set.seed(seed)
  DAT1$SS_YR <- interaction(DAT1$SS, DAT1$YEAR, drop=TRUE)
  DAT1$IDMAP <- seq_len(nrow(DAT1))
  ## randomize input
  DAT1 <- DAT1[sample.int(nrow(DAT1)),]
  #    kk <- floor(nrow(DAT1) * (1-out))
  #    DAT1k <- DAT1[1:kk,] # k as in *k*eep
  kkk <- floor(nlevels(DAT1$SS_YR) * (1-out))
  DAT1k <- DAT1[as.integer(DAT1$SS_YR) <= kkk,]
  if (nlevels(droplevels(DAT1k$bootg)) != nlevels(droplevels(DAT1$bootg)))
    stop("bootg problem: pick larger blocks for validation")
  ## one run
  r1fun <- function(DAT1k, replace=FALSE) {
    ## get rid of resamples
    DAT1k <- DAT1k[sample.int(nrow(DAT1k)),]
    DAT1k <- nonDuplicated(DAT1k, SS_YR)
    id2 <- list()
    for (l in levels(DAT1k$bootg)) {
      sset0 <- which(DAT1k$bootg == l)
      id2[[l]] <- if (length(sset0) < 2)
        sset0 else sample(sset0, length(sset0), replace=replace)
    }
    DAT1k$IDMAP[unname(unlist(id2))]
  }
  BB0 <- r1fun(DAT1k, replace=FALSE)
  BB1 <- pbsapply(seq_len(B), function(i) r1fun(DAT1k, replace=TRUE))
  cbind(BB0, BB1)
  #aa <- unique(BB1)
  #table(selected=DAT1$IDMAP %in% aa)
  #table(selected=DAT1$IDMAP %in% aa, revisit=duplicated(DAT1$SS_YR))
}

cawa_mm<-cbind(count, offset,DAT)

B<- 239
BB <- bbfun2(cawa_mm, B, out=0.95)

thin_cawa_mm <-cawa_mm[BB[,1],]

nrow(BB)/nrow(cawa_mm)
aa <- unique(BB)
bb <- table(selected=seq_len(nrow(cawa_mm)) %in% aa)
bb[2]/sum(bb)




ol<- optilevels(y=thin_cawa_mm$count, x=thin_cawa_mm$HAB, dist="poisson", offset=thin_cawa_mm$offset)
thin_cawa_mm$HAB2<-thin_cawa_mm$HAB
levels(thin_cawa_mm$HAB2)<-ol$level[[length(ol$level)]] 
table(thin_cawa_mm$HAB,thin_cawa_mm$HAB2)


mmsp <- SpatialPointsDataFrame(coords=cbind(thin_cawa_mm$X,thin_cawa_mm$Y),proj4string = CRS("+proj=longlat +ellps=WGS84"),data=thin_cawa_mm)

spplot(pres.ggwr1, "classification", do.log = F,
       key.space=list(x=0.2,y=0.3,corner=c(0,1)),
       scales=list(draw=T),
       main="Classification",
       sp.layout=list(basemap))


## Specify distance matrix:

 # view data points
grd <- SpatialGrid(GridTopology(c(-164,40),c(0.5,0.5),c(225,56)),proj4string=CRS("+proj=longlat +ellps=WGS84"))
plot(grd, axes=TRUE)
points(coordinates(mmsp))

# create distance matrix
DM <- gw.dist(dp.locat = coordinates(mmsp), longlat=TRUE)

# a simple models with HAB and ROAD only

# obtain optimum bandwidth:
bw.ggwr.1 <- bw.ggwr(count ~ HAB + ROAD + offset(mmsp$offset), data = mmsp, approach = "AICc", kernel = "bisquare", adaptive = TRUE, family="poisson", longlat = TRUE, dMat=DM) 

load("C:/Users/voeroesd/Dropbox/BAM/Critical Habitat/CHID subunit delineation/subunits.RData")

# fit model
ggwr.1<-ggwr.basic(count ~ HAB + ROAD + offset(mmsp$offset), data = mmsp, bw = bw.ggwr.1, kernel = "bisquare", adaptive = TRUE, longlat=TRUE, family="poisson", dMat = DM)

# results
print(ggwr.1)





# with climate variables
## model with many covariates, seems to need a larger sample size



bw.ggwr.2 <- bw.ggwr(count ~  HAB2 + ROAD + HGT + HGT2 + CTI + CTI2 + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 + CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT + offset(mmsp$offset), data = mmsp, approach = "AICc", kernel = "exponential", adaptive = TRUE, family="poisson", longlat = TRUE, dMat=DM) 
save.image("C:/Users/voeroesd/Dropbox/BAM/Critical Habitat/CHID subunit delineation/subunits.RData")



# fit model - maybe try with bandwidth optimized for simpler model?
ggwr.2<-ggwr.basic(count ~ HAB2 + ROAD + HGT + HGT2 + CTI + CTI2 + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 + CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT + offset(mmsp$offset), data = mmsp, bw = bw.ggwr.2, kernel = "exponential", adaptive = TRUE, longlat=TRUE, family="poisson", dMat=DM)

save("ggwr.2","C:/Users/voeroesd/Dropbox/BAM/Critical Habitat/CHID subunit delineation/output/ggwr_0.Rdata")


print(ggwr.2)




# Run several models based on different samples of the BAM dataset

ldf <- lapply(X=c(1:9),function(x){cawa_mm[BB[,x+1],]})
names(ldf)<-1:9
nameslist<-as.character(1:9)



getGWmod <- function(x){
  library(opticut)
  library(GWmodel)
  
  thindata<-ldf[[x]]
  ol<- optilevels(y=thindata$count, x=thindata$HAB, dist="poisson", offset=thindata$offset)
  thindata$HAB2<-thindata$HAB
  levels(thindata$HAB2)<-ol$level[[length(ol$level)]] 
  
  mmsp <- SpatialPointsDataFrame(coords=cbind(thindata$X,thindata$Y),proj4string = CRS("+proj=longlat +ellps=WGS84"),data=thindata)
  
  # create distance matrix
  DM <- gw.dist(dp.locat = coordinates(mmsp), longlat=TRUE)
  
  bw <- 1675 # from bw_ggwr.2, using the exponentioal kernel on the 1st sample.
  
  ggwr<-ggwr.basic(count ~ HAB2 + ROAD + HGT + HGT2 + CTI + CTI2 + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 + CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT + offset(mmsp$offset), data = mmsp, bw = bw, kernel = "exponential", adaptive = TRUE, longlat=TRUE, family="poisson", dMat=DM)
  
  wd<-"C:/Users/voeroesd/Dropbox/BAM/Critical Habitat/CHID subunit delineation"
  #wd<-getwd()
  outf <- paste0("ggwr_",x)
  save("ggwr", file.path(wd,"output",paste0(outf,".Rdata")))
  
  return(ggwr)
}


models<- lapply(nameslist,getGWmod)




## Basemap

basemap<-readOGR("province_state_lcc.shp")
proj4string((basemap))


#### Cluster analysis
library(mclust)
str(ggwr.1)

mclust1<-Mclust(data=ggwr.1$SDF@data[which(ggwr.1$SDF@data$y>0),1:10],G=1:9) # clustering only the points where the species is present, i.e. y>0)
summary(mclust1)
summary(mclust1$BIC)


plot(mclust1)




pres.ggwr1<-ggwr.1$SDF[which(ggwr.1$SDF@data$y>0),]

pres.ggwr1$classification<-as.factor(mclust1$classification)
pres.ggwr1$uncertainty<-mclust1$uncertainty

pres.ggwr1<-spTransform(pres.ggwr1,CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"))


spplot(pres.ggwr1, "classification", do.log = F,
       key.space=list(x=0.2,y=0.3,corner=c(0,1)),
       scales=list(draw=T),
       main="Classification",
       sp.layout=list(basemap))

spplot(pres.ggwr1, "uncertainty", do.log = F,
       key.space=list(x=0.2,y=0.3,corner=c(0,1)),
       scales=list(draw=T),
       main="Uncertainty",
       colorkey=T,
       sp.layout=list(basemap))

drmclust1<-MclustDR(mclust1,lambda=1)
summary(drmclust1)
plot(drmclust1, what="contour")
plot(drmclust1, what="boundaries", ngrid=500)

