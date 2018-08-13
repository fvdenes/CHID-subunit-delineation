# CAWA national model as a Geographically-weighted Regression for Canada Warbler (North America extent) ####

# Load packages and data

library("GWmodel")
library(mefa4)
library(sp)
library(rgdal)
library(opticut)
library(raster)
library(scales)
library(mclust)
library(RColorBrewer)
library(colorspace)
library(ggplot2)
library(gridExtra)

#load("C:/Users/voeroesd/Dropbox/BAM/Critical Habitat/CHID subunit delineation/pack_2016-12-01.Rdata")
load("D:/CHID subunit delineation/subunits.RData")
basemap<-readOGR("province_state_lcc.shp") ## Basemap for point plots
BCRs<- readOGR("bcrfinallcc.shp")
brandt<-readOGR("BRANDT_diss_type.shp")
brandt<-spTransform(brandt,CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"))

#### Some data preparation ####
DAT$HAB_NALC1 <- DAT$HABTR
DAT$HAB_NALC2 <- DAT$HAB
DAT$YEAR <- DAT$YR+1997

mm <- Mefa(YY, DAT, TAX, "inner")

# Combine CAWA counts, offsets and covariates in a single dataset
count <- YY[,"CAWA"]
offset <- OFF[,"CAWA"]


### Sampling from the dataset (GW model will not run with full dataset) 

# Define year resampling blocks
DAT$YR5 <- 0
DAT$YR5[DAT$YEAR > 2000] <- 1
DAT$YR5[DAT$YEAR > 2004] <- 2
DAT$YR5[DAT$YEAR > 2009] <- 3
table(DAT$YEAR,DAT$YR5)
table(DAT$YR5)

# Define intersections
DAT$bootg <- interaction(DAT$JURS, DAT$YR5, drop=TRUE)

# Join objects
cawa_mm<-cbind(count, offset,DAT)

allpoints <- SpatialPointsDataFrame(coords=cbind(cawa_mm$X,cawa_mm$Y),proj4string =  CRS("+proj=longlat +ellps=WGS84"),data=cawa_mm)
allpoints<-spTransform(allpoints,CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"))



spplot(allpoints,"count",do.log = F,legendEntries=as.character(0:4),
       key.space=list(x=0.02,y=0.3,corner=c(0,1)),
       sp.layout=list(basemap))


spplot(allpoints[which(allpoints$count>0),],"count",do.log = F,cuts=4,legendEntries=as.character(1:4),
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       sp.layout=list(basemap))


### A function to sample points with count>0 in each jurisdiction
sample_jur<-function(data,N=20,seed=123){
  set.seed(seed)
  data2<-data[which(data$count>0),]
  ls<-as.list(rep(NA,length(levels(data2$JURS))))
  for(i in 1:length(levels(data2$JURS))){
    data22<-data2[which(data2$JURS==levels(data2$JURS)[i]),]
    if(table(data2$JURS)[i]>0){
      if(table(data2$JURS)[i]<N){
        ls[[i]]<-row.names(data22)
      }
      else{
        ls[[i]]<-row.names(data22[sample(nrow(data22),size=N),]) 
      }
    }
  }
  ls<-unlist(ls)
  ls<-ls[-which(is.na(ls))]
  out<-match(x=ls,table = row.names(data))
}

#### A function to randomly sample from jurisdiction
sample_max<-function(data,N=1000,seed=123){
  set.seed(seed)
  ls<-as.list(rep(NA,length(levels(data$JURS))))
  for(i in 1:length(levels(data$JURS))){
    data2<-data[which(data$JURS==levels(data$JURS)[i]),]
    if(table(data$JURS)[i]>0){
      if(table(data$JURS)[i]<N){
        ls[[i]]<-row.names(data2)
      }
      else{
        ls[[i]]<-row.names(data2[sample(nrow(data2),size=N),]) 
      }
    }
  }
  ls<-unlist(ls)
  out<-match(x=ls,table = row.names(data))
}

#### A function to randomly sample from longitude intervals
sample_long<-function(data,N=200,intervals=100,seed=123){
  set.seed(seed)
  maxlongitude<-max(data$X)
  minlongitude<-min(data$X)
  Xwidth<-(maxlongitude-minlongitude)/intervals
  Xlimits<-minlongitude+c(0,Xwidth*1:intervals)
  ls<-as.list(rep(NA,intervals))
  for(i in 1:intervals){
    data2<-data[which(data$X>Xlimits[i] & data$X<Xlimits[i+1]),]
    if(nrow(data2)<N){
      ls[[i]]<-row.names(data2)
    }
    else{
      ls[[i]]<-row.names(data2[sample(nrow(data2),size=N),]) 
    }
  }
  ls<-unlist(ls)
  out<-match(x=ls,table = row.names(data))
}

#### A function to randomly sample from quadrants defined by latitude and longitude intervals.
sample_quad<-function(data,N=20,Xintervals=100,Yintervals=60,seed=123,onlypres=F){
  set.seed(seed)
  if(onlypres==F){
    maxlongitude<-max(data$X)
    minlongitude<-min(data$X)
    Xwidth<-(maxlongitude-minlongitude)/Xintervals
    Xlimits<-minlongitude+c(0,Xwidth*1:Xintervals)
    
    maxlatitude<-max(data$Y)
    minlatitude<-min(data$Y)
    Ywidth<-(maxlatitude-minlatitude)/Yintervals
    Ylimits<-minlatitude+c(0,Ywidth*1:Yintervals)
    
    ls<-array(NA,dim=c(Xintervals,Yintervals,N))
    for(i in 1:Xintervals){
      data2<-data[which(data$X>Xlimits[i] & data$X<Xlimits[i+1]),]
      for(j in 1:Yintervals){
        data3<-data2[which(data2$Y>Ylimits[j] & data2$Y<Ylimits[j+1]),]
        if(nrow(data3)>0){
          if(nrow(data3)<N){
            ls[i,j,1:nrow(data3)]<-row.names(data3)
          }
          else{
            ls[i,j,]<-row.names(data3[sample(nrow(data3),size=N),]) 
          }
        }
      }
    } 
  }
  
  else{
    data2<-data[which(data$count>0),]
    maxlongitude<-max(data2$X)
    minlongitude<-min(data2$X)
    Xwidth<-(maxlongitude-minlongitude)/Xintervals
    Xlimits<-minlongitude+c(0,Xwidth*1:Xintervals)
    
    maxlatitude<-max(data2$Y)
    minlatitude<-min(data2$Y)
    Ywidth<-(maxlatitude-minlatitude)/Yintervals
    Ylimits<-minlatitude+c(0,Ywidth*1:Yintervals)
    
    ls<-array(NA,dim=c(Xintervals,Yintervals,N))
    for(i in 1:Xintervals){
      data3<-data2[which(data2$X>Xlimits[i] & data2$X<Xlimits[i+1]),]
      for(j in 1:Yintervals){
        data4<-data3[which(data3$Y>Ylimits[j] & data3$Y<Ylimits[j+1]),]
        if(nrow(data4)>0){
          if(nrow(data4)<N){
            ls[i,j,1:nrow(data4)]<-row.names(data4)
          }
          else{
            ls[i,j,]<-row.names(data4[sample(nrow(data4),size=N),]) 
          }
        }
      }
    } 
  }
 
    ls<-c(ls)
    ls<-ls[-which(is.na(ls))]
  out<-match(x=ls,table = row.names(data))
}


plot_sampled<-function(obj,index,onlypres=T,show="count"){
  thin_cawa_mm<-obj[index,]
  mmsp <- SpatialPointsDataFrame(coords=cbind(thin_cawa_mm$X,thin_cawa_mm$Y),proj4string =  CRS("+proj=longlat +ellps=WGS84"),data=thin_cawa_mm)
  mmsp<-spTransform(mmsp,CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"))
  if(onlypres==T){
    spplot(mmsp[which(mmsp$count>0),],show,do.log = F,cuts=4,legendEntries=as.character(1:4),
           key.space=list(x=0.5,y=0.9,corner=c(0,1)),
           sp.layout=list(basemap))
  }
  else{
    spplot(mmsp,show,do.log = F,legendEntries=as.character(0:4),
           key.space=list(x=0.02,y=0.3,corner=c(0,1)),
           sp.layout=list(basemap))
  }
}

sample<-sample_jur(cawa_mm, N=30,seed=222)
plot_sampled(obj=cawa_mm, index=sample, onlypres = T) 



sample2<-sample_long(cawa_mm)
plot_sampled(obj=cawa_mm, index=sample2, onlypres = F)
plot_sampled(obj=cawa_mm, index=sample2, onlypres = F, show="YEAR")

sample3<-sample_quad(cawa_mm)
sample4<-sample_quad(cawa_mm,onlypres = T, N=4)

plot_sampled(obj=cawa_mm, index=sample3, onlypres = F)
plot_sampled(obj=cawa_mm, index=sample4, onlypres = T)

#### Check spatial distribution of sampled points with count >0 - make sure they are well-distributed ####
union<-union(sample3,sample4)

# function to plot points where CAWA count >0 from samples
plot_sampled(obj=cawa_mm, index=union, onlypres = F) 
plot_sampled(obj=cawa_mm, index=union, onlypres = T) 

table(subset(thin_cawa_mm,JURS=="SK")$YEAR)
thin_cawa_mm <-cawa_mm[union,]
table(thin_cawa_mm$JURS)
table(thin_cawa_mm$YEAR)

# clear some of the workspace
rm(list=ls()[! ls() %in% c("cawa_mm","thin_cawa_mm","basemap","BCRs","brandt","sample3","sample4","plot_sampled")])         
           
# apply landscape covariate optimization algorithm
ol<- optilevels(y=thin_cawa_mm$count, x=thin_cawa_mm$HAB_NALC1, dist="poisson", offset=thin_cawa_mm$offset)
thin_cawa_mm$HAB_NALC2<-thin_cawa_mm$HAB_NALC1
levels(thin_cawa_mm$HAB_NALC2)<-ol$level[[length(ol$level)]] 
table(thin_cawa_mm$HAB_NALC1,thin_cawa_mm$HAB_NALC2)

# create spatial dataframe from sampled dataset
mmsp <- SpatialPointsDataFrame(coords=cbind(thin_cawa_mm$X,thin_cawa_mm$Y),proj4string = CRS("+proj=longlat +ellps=WGS84"),data=thin_cawa_mm)

#### GW modeling and clustering ####
## Specify distance matrix:
DM <- gw.dist(dp.locat = coordinates(mmsp), longlat=TRUE)

# Optimize bandwidth - exponential kernel function
bw.ggwr.exponential_grid <- bw.ggwr(count ~  HAB_NALC2 + ROAD + HGT + HGT2 + CTI + CTI2 + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 + CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT + offset(mmsp$offset), data = mmsp, approach = "AICc", kernel = "exponential", adaptive = TRUE, family="poisson", longlat = TRUE, dMat=DM) 

# Fit model 
ggwr_exponential_grid<-ggwr.basic(count ~ HAB_NALC2 + ROAD + HGT + HGT2 + CTI + CTI2 + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 + CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT + offset(mmsp$offset), data = mmsp, bw = bw.ggwr.exponential_grid, kernel = "exponential", adaptive = TRUE, longlat=TRUE, family="poisson", dMat=DM)

save.image("D:/CHID subunit delineation/subunits.RData")
save("ggwr_exponential_grid",file="D:/CHID subunit delineation//output/ggwr_exponential_grid.Rdata")
# load("ggwr_exponential","C:/Users/voeroesd/Dropbox/BAM/Critical Habitat/CHID subunit delineation/output/ggwr_exponential.Rdata")

print(ggwr_exponential_grid)


summary(ggwr_exponential_grid$SDF)
ggwr_exponential_grid$SDF$y

#### Clustering 1-9 clusters####
mclust_1<- Mclust(data=ggwr_exponential_grid$SDF@data[which(ggwr_exponential_grid$SDF@data$y>0),c(1:15)[-11]],G=1:9) # vegetation type, tree height and topography
summary(mclust_1)

plot(mclust_1, "BIC")
plot(mclust_1,"classification", dimens=c(12,13),symbols=16)


###Plot survey points where count>0 coloured by cluster

# Load basemap layers
#basemap<-readOGR("province_state_lcc.shp") #provinces and states
#BCRs<- readOGR("bcrfinallcc.shp") # BCRs

# Create spatial dataframe for points to be plotted
pres.ggwr_exponential<-ggwr_exponential_grid$SDF[which(ggwr_exponential_grid$SDF@data$y>0),]
pres.ggwr_exponential$classification<-as.factor(mclust_1$classification)
pres.ggwr_exponential$uncertainty<-mclust_1$uncertainty

# Reproject
pres.ggwr_exponential<-spTransform(pres.ggwr_exponential,CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"))

# Create simple plot of points to obtain plot extents, and crop BCR layer to that extent
p<- spplot(pres.ggwr_exponential, "classification", do.log = F)
p.extent<-extent(rbind(p$x.limits,p$y.limits))
out<-crop(BCRs,p.extent)
out$BCR<-as.factor(out$BCR)
levels(out$BCR)<-1:21


out2<-crop(brandt,p.extent)

set.seed(12)
colsmapBCRs<-c("#ffffff",sample(rainbow_hcl(21,alpha=0.8))[-1])[out$BCR]

# with BCRs
spplot(pres.ggwr_exponential, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=brewer.pal(9,"Set1"),
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out,col="gainsboro",fill=colsmapBCRs),basemap)
)

set.seed(432)
brandtcols<-sample(rainbow_hcl(4,alpha=0.8))
# with Brandt boreal
spplot(pres.ggwr_exponential, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=brewer.pal(9,"Set1"),
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out2,col="gainsboro",fill=brandtcols),basemap)
)


spplot(pres.ggwr_exponential, "uncertainty", do.log = F,
       key.space=list(x=0.2,y=0.3,corner=c(0,1)),
       main="Uncertainty",
       colorkey=T,
       sp.layout=list(basemap))


# Summary of coefficients per clusters


clust.coef.means<-as.data.frame(t(aggregate(pres.ggwr_exponential@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential@data$classification),mean)[,-1]))
clust.coef.means<-cbind(Coef=factor(rownames(clust.coef.means),levels=rownames(clust.coef.means)),clust.coef.means)

clust.coef.means$Coef<-as.character(clust.coef.means$Coef)
clust.coef.means$Coef[2:10]<-substring(as.character(clust.coef.means$Coef)[2:10],first=10)
clust.coef.means$Coef<-factor(clust.coef.means$Coef,levels=clust.coef.means$Coef)
colnames(clust.coef.means)[2:10]<-paste0("mean",1:9)

clust.coef.medians<-cbind(Coef=clust.coef.means$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential@data$classification),median)[,-1])))
clust.coef.sd<-cbind(Coef=clust.coef.means$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential@data$classification),sd)[,-1])))
colnames(clust.coef.sd)[2:10]<-paste0("sd",1:9)
colnames(clust.coef.medians)[2:10]<-paste0("median",1:9)


clust.summary<-cbind(clust.coef.means,clust.coef.medians[,-1],clust.coef.sd[,-1])

p1<-ggplot(clust.summary,aes(Coef,mean1))
p1+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(9,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0.1))+theme(axis.text.x = element_text(size=12,angle = 45, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(9,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.15))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(9,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(9,"Set1")[4],size=1,fatten=2,position=position_nudge(x=-0.2))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(9,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.05))+ #orange
  geom_pointrange(aes(Coef,mean6,ymin=mean6-sd6,ymax=mean6+sd6),col=brewer.pal(9,"Set1")[6],size=1,fatten=2,position=position_nudge(x=0.2))+#yellow
  geom_pointrange(aes(Coef,mean7,ymin=mean7-sd7,ymax=mean7+sd7),col=brewer.pal(9,"Set1")[7],size=1,fatten=2,position=position_nudge(x=0.15))+ #brown
  geom_pointrange(aes(Coef,mean8,ymin=mean8-sd8,ymax=mean8+sd8),col=brewer.pal(9,"Set1")[8],size=1,fatten=2,position=position_nudge(x=-0.1))+ #pink
  geom_pointrange(aes(Coef,mean9,ymin=mean9-sd9,ymax=mean9+sd9),col=brewer.pal(9,"Set1")[9],size=1,fatten=2,position=position_nudge(x=-0.05)) #grey

# Draw convex hull polygons from clusters
cluster1<-pres.ggwr_exponential@coords[which(pres.ggwr_exponential$classification==1),]
ch<-chull(cluster1)
coords1<-cluster1[c(ch,ch[1]),]
spplot(cluster1, pch=16)
lines(coords1, col="red")

spplot(pres.ggwr_exponential[which(pres.ggwr_exponential$classification==1),],"classification", sp.layout=basemap)


spplot(pres.ggwr_exponential, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=brewer.pal(9,"Set1"),
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out2,col="gainsboro",fill=brandtcols),basemap)
)





### Limiting number of clusters to 3: ####
mclust_3<- Mclust(data=ggwr_exponential_grid$SDF@data[which(ggwr_exponential_grid$SDF@data$y>0),c(1:15)[-11]],G=1:3) # vegetation type, tree height and topography
summary(mclust_3)

plot(mclust_3, "BIC")

pres.ggwr_exponential3<-ggwr_exponential_grid$SDF[which(ggwr_exponential_grid$SDF@data$y>0),]
pres.ggwr_exponential3$classification<-as.factor(mclust_3$classification)
pres.ggwr_exponential3$uncertainty<-mclust_3$uncertainty
pres.ggwr_exponential3<-spTransform(pres.ggwr_exponential3,CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"))

# with BCRs
spplot(pres.ggwr_exponential3, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=brewer.pal(4,"Set1")[-3],
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out,col="gainsboro",fill=colsmapBCRs),basemap)
)

set.seed(432)
brandtcols<-sample(rainbow_hcl(4,alpha=0.8))
# with Brandt boreal
spplot(pres.ggwr_exponential3, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=brewer.pal(4,"Set1")[-3],
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out2,col="gainsboro",fill=brandtcols),basemap)
)

# Summary of coefficients per clusters
clust.coef.means3<-as.data.frame(t(aggregate(pres.ggwr_exponential3@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential3@data$classification),mean)[,-1]))
clust.coef.means3<-cbind(Coef=factor(rownames(clust.coef.means3),levels=rownames(clust.coef.means3)),clust.coef.means3)

clust.coef.means3$Coef<-as.character(clust.coef.means3$Coef)
clust.coef.means3$Coef[2:10]<-substring(as.character(clust.coef.means3$Coef)[2:10],first=10)
clust.coef.means3$Coef<-factor(clust.coef.means3$Coef,levels=clust.coef.means3$Coef)
colnames(clust.coef.means3)[2:4]<-paste0("mean",1:3)

clust.coef.medians3<-cbind(Coef=clust.coef.means3$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential3@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential3@data$classification),median)[,-1])))
clust.coef.sd3<-cbind(Coef=clust.coef.means3$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential3@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential3@data$classification),sd)[,-1])))
colnames(clust.coef.sd3)[2:4]<-paste0("sd",1:3)
colnames(clust.coef.medians3)[2:4]<-paste0("median",1:3)


clust.summary3<-cbind(clust.coef.means3,clust.coef.medians3[,-1],clust.coef.sd3[,-1])

p3<-ggplot(clust.summary3,aes(Coef,mean1))
p3+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(4,"Set1")[1],size=1,fatten=2)+theme(axis.text.x = element_text(size=15,angle = 45, hjust = 1),axis.text.y = element_text(size=15))+ylab("")+xlab("")+
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(4,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.05))+
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(4,"Set1")[4],size=1,fatten=2,position=position_nudge(x=0.05))





### Limiting number of clusters to 4: ####
mclust_4<- Mclust(data=ggwr_exponential_grid$SDF@data[which(ggwr_exponential_grid$SDF@data$y>0),c(1:15)[-11]],G=1:4) # vegetation type, tree height and topography
summary(mclust_4)

plot(mclust_4, "BIC")
plot(mclust_4,"classification", dimens=c(12,13),symbols=16)


###Plot survey points where count>0 coloured by cluster

# Load basemap layers
#basemap<-readOGR("province_state_lcc.shp") #provinces and states
#BCRs<- readOGR("bcrfinallcc.shp") # BCRs

# Create spatial dataframe for points to be plotted
pres.ggwr_exponential4<-ggwr_exponential_grid$SDF[which(ggwr_exponential_grid$SDF@data$y>0),]
pres.ggwr_exponential4$classification<-as.factor(mclust_4$classification)
pres.ggwr_exponential4$uncertainty<-mclust_4$uncertainty

# Reproject
pres.ggwr_exponential4<-spTransform(pres.ggwr_exponential4,CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"))


# with BCRs
spplot(pres.ggwr_exponential4, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=brewer.pal(4,"Set1"),
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out,col="gainsboro",fill=colsmapBCRs),basemap)
)

set.seed(432)
brandtcols<-sample(rainbow_hcl(4,alpha=0.8))
# with Brandt boreal
spplot(pres.ggwr_exponential4, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=brewer.pal(4,"Set1"),
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out2,col="gainsboro",fill=brandtcols),basemap)
)

# Summary of coefficients per clusters
clust.coef.means4<-as.data.frame(t(aggregate(pres.ggwr_exponential4@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential4@data$classification),mean)[,-1]))
clust.coef.means4<-cbind(Coef=factor(rownames(clust.coef.means4),levels=rownames(clust.coef.means4)),clust.coef.means4)

clust.coef.means4$Coef<-as.character(clust.coef.means4$Coef)
clust.coef.means4$Coef[2:10]<-substring(as.character(clust.coef.means4$Coef)[2:10],first=10)
clust.coef.means4$Coef<-factor(clust.coef.means4$Coef,levels=clust.coef.means4$Coef)
colnames(clust.coef.means4)[2:5]<-paste0("mean",1:4)

clust.coef.medians4<-cbind(Coef=clust.coef.means4$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential4@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential4@data$classification),median)[,-1])))
clust.coef.sd4<-cbind(Coef=clust.coef.means4$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential4@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential4@data$classification),sd)[,-1])))
colnames(clust.coef.sd4)[2:5]<-paste0("sd",1:4)
colnames(clust.coef.medians4)[2:5]<-paste0("median",1:4)


clust.summary4<-cbind(clust.coef.means4,clust.coef.medians4[,-1],clust.coef.sd4[,-1])

p4<-ggplot(clust.summary4,aes(Coef,mean1))
p4+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(9,"Set1")[1],size=1,fatten=2)+theme(axis.text.x = element_text(size=15,angle = 45, hjust = 1),axis.text.y = element_text(size=15))+ylab("")+xlab("")+
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd4,ymax=mean2+sd4),col=brewer.pal(9,"Set1")[2],size=1,fatten=2,position=position_nudge(x=0.05))+
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(9,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0.1))+
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(9,"Set1")[4],size=1,fatten=2,position=position_nudge(x=0.15))



### Limiting number of clusters to 5: ####
mclust_5<- Mclust(data=ggwr_exponential_grid$SDF@data[which(ggwr_exponential_grid$SDF@data$y>0),c(1:15)[-11]],G=1:5) # vegetation type, tree height and topography
summary(mclust_5)

plot(mclust_5, "BIC")

pres.ggwr_exponential5<-ggwr_exponential_grid$SDF[which(ggwr_exponential_grid$SDF@data$y>0),]
pres.ggwr_exponential5$classification<-as.factor(mclust_5$classification)
pres.ggwr_exponential5$uncertainty<-mclust_5$uncertainty
pres.ggwr_exponential5<-spTransform(pres.ggwr_exponential5,CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"))

# with BCRs
spplot(pres.ggwr_exponential5, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=brewer.pal(5,"Set1"),
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out,col="gainsboro",fill=colsmapBCRs),basemap)
)

set.seed(432)
brandtcols<-sample(rainbow_hcl(4,alpha=0.8))
# with Brandt boreal
spplot(pres.ggwr_exponential5, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=brewer.pal(5,"Set1"),
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out2,col="gainsboro",fill=brandtcols),basemap)
)

# Summary of coefficients per clusters
clust.coef.means5<-as.data.frame(t(aggregate(pres.ggwr_exponential5@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential5@data$classification),mean)[,-1]))
clust.coef.means5<-cbind(Coef=factor(rownames(clust.coef.means5),levels=rownames(clust.coef.means5)),clust.coef.means5)

clust.coef.means5$Coef<-as.character(clust.coef.means5$Coef)
clust.coef.means5$Coef[2:10]<-substring(as.character(clust.coef.means5$Coef)[2:10],first=10)
clust.coef.means5$Coef<-factor(clust.coef.means5$Coef,levels=clust.coef.means5$Coef)
colnames(clust.coef.means5)[2:6]<-paste0("mean",1:5)

clust.coef.medians5<-cbind(Coef=clust.coef.means5$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential5@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential5@data$classification),median)[,-1])))
clust.coef.sd5<-cbind(Coef=clust.coef.means5$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential5@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential5@data$classification),sd)[,-1])))
colnames(clust.coef.sd5)[2:6]<-paste0("sd",1:5)
colnames(clust.coef.medians5)[2:6]<-paste0("median",1:5)


clust.summary5<-cbind(clust.coef.means5,clust.coef.medians5[,-1],clust.coef.sd5[,-1])

p5<-ggplot(clust.summary5,aes(Coef,mean1))
p5+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2)+theme(axis.text.x = element_text(size=15,angle = 45, hjust = 1),axis.text.y = element_text(size=15))+ylab("")+xlab("")+
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.05))+
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=-0.1))+
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(5,"Set1")[4],size=1,fatten=2,position=position_nudge(x=0.05))+
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(5,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.1))



### Limiting number of clusters to 6: ####
mclust_6<- Mclust(data=ggwr_exponential_grid$SDF@data[which(ggwr_exponential_grid$SDF@data$y>0),c(1:15)[-11]],G=1:6) # vegetation type, tree height and topography
summary(mclust_6)

plot(mclust_6, "BIC")

pres.ggwr_exponential6<-ggwr_exponential_grid$SDF[which(ggwr_exponential_grid$SDF@data$y>0),]
pres.ggwr_exponential6$classification<-as.factor(mclust_6$classification)
pres.ggwr_exponential6$uncertainty<-mclust_6$uncertainty
pres.ggwr_exponential6<-spTransform(pres.ggwr_exponential6,CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"))

# with BCRs
spplot(pres.ggwr_exponential6, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=brewer.pal(6,"Set1"),
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out,col="gainsboro",fill=colsmapBCRs),basemap)
)

set.seed(432)
brandtcols<-sample(rainbow_hcl(4,alpha=0.8))
# with Brandt boreal
spplot(pres.ggwr_exponential6, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=brewer.pal(6,"Set1"),
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out2,col="gainsboro",fill=brandtcols),basemap)
)

# Summary of coefficients per clusters
clust.coef.means6<-as.data.frame(t(aggregate(pres.ggwr_exponential6@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential6@data$classification),mean)[,-1]))
clust.coef.means6<-cbind(Coef=factor(rownames(clust.coef.means6),levels=rownames(clust.coef.means6)),clust.coef.means6)

clust.coef.means6$Coef<-as.character(clust.coef.means6$Coef)
clust.coef.means6$Coef[2:10]<-substring(as.character(clust.coef.means6$Coef)[2:10],first=10)
clust.coef.means6$Coef<-factor(clust.coef.means6$Coef,levels=clust.coef.means6$Coef)
colnames(clust.coef.means6)[2:7]<-paste0("mean",1:6)

clust.coef.medians6<-cbind(Coef=clust.coef.means6$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential6@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential6@data$classification),median)[,-1])))
clust.coef.sd6<-cbind(Coef=clust.coef.means6$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential6@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential6@data$classification),sd)[,-1])))
colnames(clust.coef.sd6)[2:7]<-paste0("sd",1:6)
colnames(clust.coef.medians6)[2:7]<-paste0("median",1:6)


clust.summary6<-cbind(clust.coef.means6,clust.coef.medians6[,-1],clust.coef.sd6[,-1])

p6<-ggplot(clust.summary6,aes(Coef,mean1))
p6+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(6,"Set1")[1],size=1,fatten=2)+theme(axis.text.x = element_text(size=15,angle = 45, hjust = 1),axis.text.y = element_text(size=15))+ylab("")+xlab("")+
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(6,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.05))+
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(6,"Set1")[3],size=1,fatten=2,position=position_nudge(x=-0.1))+
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(6,"Set1")[4],size=1,fatten=2,position=position_nudge(x=0.05))+
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(6,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.1))+
  geom_pointrange(aes(Coef,mean6,ymin=mean6-sd6,ymax=mean6+sd6),col=brewer.pal(6,"Set1")[6],size=1,fatten=2,position=position_nudge(x=0.15))



### Limiting number of clusters to 7: ####
mclust_7<- Mclust(data=ggwr_exponential_grid$SDF@data[which(ggwr_exponential_grid$SDF@data$y>0),c(1:15)[-11]],G=1:7) # vegetation type, tree height and topography
summary(mclust_7)

plot(mclust_7, "BIC")

pres.ggwr_exponential7<-ggwr_exponential_grid$SDF[which(ggwr_exponential_grid$SDF@data$y>0),]
pres.ggwr_exponential7$classification<-as.factor(mclust_7$classification)
pres.ggwr_exponential7$uncertainty<-mclust_7$uncertainty
pres.ggwr_exponential7<-spTransform(pres.ggwr_exponential7,CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"))

# with BCRs
spplot(pres.ggwr_exponential7, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=brewer.pal(9,"Set1")[-(7:8)],
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out,col="gainsboro",fill=colsmapBCRs),basemap)
)

set.seed(432)
brandtcols<-sample(rainbow_hcl(4,alpha=0.8))
# with Brandt boreal
spplot(pres.ggwr_exponential7, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=brewer.pal(7,"Set1"),
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out2,col="gainsboro",fill=brandtcols),basemap)
)

# Summary of coefficients per clusters
clust.coef.means7<-as.data.frame(t(aggregate(pres.ggwr_exponential7@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential7@data$classification),mean)[,-1]))
clust.coef.means7<-cbind(Coef=factor(rownames(clust.coef.means7),levels=rownames(clust.coef.means7)),clust.coef.means7)

clust.coef.means7$Coef<-as.character(clust.coef.means7$Coef)
clust.coef.means7$Coef[2:10]<-substring(as.character(clust.coef.means7$Coef)[2:10],first=10)
clust.coef.means7$Coef<-factor(clust.coef.means7$Coef,levels=clust.coef.means7$Coef)
colnames(clust.coef.means7)[2:8]<-paste0("mean",1:7)

clust.coef.medians7<-cbind(Coef=clust.coef.means7$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential7@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential7@data$classification),median)[,-1])))
clust.coef.sd7<-cbind(Coef=clust.coef.means7$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential7@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential7@data$classification),sd)[,-1])))
colnames(clust.coef.sd7)[2:8]<-paste0("sd",1:7)
colnames(clust.coef.medians7)[2:8]<-paste0("median",1:7)


clust.summary7<-cbind(clust.coef.means7,clust.coef.medians7[,-1],clust.coef.sd7[,-1])

p7<-ggplot(clust.summary7,aes(Coef,mean1))
p7+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(7,"Set1")[1],size=1,fatten=2)+theme(axis.text.x = element_text(size=15,angle = 45, hjust = 1),axis.text.y = element_text(size=15))+ylab("")+xlab("")+
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(7,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.05))+
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(7,"Set1")[3],size=1,fatten=2,position=position_nudge(x=-0.1))+
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(7,"Set1")[4],size=1,fatten=2,position=position_nudge(x=-0.15))+
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(7,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.1))+
  geom_pointrange(aes(Coef,mean6,ymin=mean6-sd6,ymax=mean6+sd6),col=brewer.pal(7,"Set1")[6],size=1,fatten=2,position=position_nudge(x=0.05))+
  geom_pointrange(aes(Coef,mean6,ymin=mean7-sd7,ymax=mean7+sd7),col=brewer.pal(7,"Set1")[7],size=1,fatten=2,position=position_nudge(x=0.15))




### Limiting number of clusters to 8: ####
mclust_8<- Mclust(data=ggwr_exponential_grid$SDF@data[which(ggwr_exponential_grid$SDF@data$y>0),c(1:15)[-11]],G=1:8) # vegetation type, tree height and topography
summary(mclust_8)

plot(mclust_8, "BIC")

pres.ggwr_exponential8<-ggwr_exponential_grid$SDF[which(ggwr_exponential_grid$SDF@data$y>0),]
pres.ggwr_exponential8$classification<-as.factor(mclust_8$classification)
pres.ggwr_exponential8$uncertainty<-mclust_8$uncertainty
pres.ggwr_exponential8<-spTransform(pres.ggwr_exponential8,CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"))

# with BCRs
spplot(pres.ggwr_exponential8, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=brewer.pal(8,"Set1"),
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out,col="gainsboro",fill=colsmapBCRs),basemap)
)

set.seed(432)
brandtcols<-sample(rainbow_hcl(4,alpha=0.8))
# with Brandt boreal
spplot(pres.ggwr_exponential8, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=brewer.pal(8,"Set1"),
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out2,col="gainsboro",fill=brandtcols),basemap)
)

# Summary of coefficients per clusters
clust.coef.means8<-as.data.frame(t(aggregate(pres.ggwr_exponential8@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential8@data$classification),mean)[,-1]))
clust.coef.means8<-cbind(Coef=factor(rownames(clust.coef.means8),levels=rownames(clust.coef.means8)),clust.coef.means8)

clust.coef.means8$Coef<-as.character(clust.coef.means8$Coef)
clust.coef.means8$Coef[2:10]<-substring(as.character(clust.coef.means8$Coef)[2:10],first=10)
clust.coef.means8$Coef<-factor(clust.coef.means8$Coef,levels=clust.coef.means8$Coef)
colnames(clust.coef.means8)[2:9]<-paste0("mean",1:8)

clust.coef.medians8<-cbind(Coef=clust.coef.means8$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential8@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential8@data$classification),median)[,-1])))
clust.coef.sd8<-cbind(Coef=clust.coef.means8$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential8@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential8@data$classification),sd)[,-1])))
colnames(clust.coef.sd8)[2:9]<-paste0("sd",1:8)
colnames(clust.coef.medians8)[2:9]<-paste0("median",1:8)


clust.summary8<-cbind(clust.coef.means8,clust.coef.medians8[,-1],clust.coef.sd8[,-1])

p8<-ggplot(clust.summary8,aes(Coef,mean1))
p8+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(8,"Set1")[1],size=1,fatten=2)+theme(axis.text.x = element_text(size=15,angle = 45, hjust = 1),axis.text.y = element_text(size=15))+ylab("")+xlab("")+
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(8,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.05))+
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(8,"Set1")[3],size=1,fatten=2,position=position_nudge(x=-0.1))+
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(8,"Set1")[4],size=1,fatten=2,position=position_nudge(x=-0.15))+
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(8,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.1))+
  geom_pointrange(aes(Coef,mean6,ymin=mean6-sd6,ymax=mean6+sd6),col=brewer.pal(8,"Set1")[6],size=1,fatten=2,position=position_nudge(x=0.05))+
  geom_pointrange(aes(Coef,mean7,ymin=mean7-sd7,ymax=mean7+sd7),col=brewer.pal(8,"Set1")[7],size=1,fatten=2,position=position_nudge(x=0.15))+
  geom_pointrange(aes(Coef,mean8,ymin=mean8-sd8,ymax=mean8+sd8),col=brewer.pal(8,"Set1")[8],size=1,fatten=2,position=position_nudge(x=0.15))



### Limiting number of clusters to 10: ####
mclust_10<- Mclust(data=ggwr_exponential_grid$SDF@data[which(ggwr_exponential_grid$SDF@data$y>0),c(1:15)[-11]],G=1:10) # vegetation type, tree height and topography
summary(mclust_10)

plot(mclust_10, "BIC")

pres.ggwr_exponential10<-ggwr_exponential_grid$SDF[which(ggwr_exponential_grid$SDF@data$y>0),]
pres.ggwr_exponential10$classification<-as.factor(mclust_10$classification)
pres.ggwr_exponential10$uncertainty<-mclust_10$uncertainty
pres.ggwr_exponential10<-spTransform(pres.ggwr_exponential10,CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"))

# with BCRs
spplot(pres.ggwr_exponential10, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=cores[c(1:5,7:11)],
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out,col="gainsboro",fill=colsmapBCRs),basemap)
)

set.seed(432)
brandtcols<-sample(rainbow_hcl(4,alpha=0.8))
# with Brandt boreal
spplot(pres.ggwr_exponential10, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=cores[c(1:5,7:11)],
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out2,col="gainsboro",fill=brandtcols),basemap)
)

# Summary of coefficients per clusters
clust.coef.means10<-as.data.frame(t(aggregate(pres.ggwr_exponential10@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential10@data$classification),mean)[,-1]))
clust.coef.means10<-cbind(Coef=factor(rownames(clust.coef.means10),levels=rownames(clust.coef.means10)),clust.coef.means10)

clust.coef.means10$Coef<-as.character(clust.coef.means10$Coef)
clust.coef.means10$Coef[2:10]<-substring(as.character(clust.coef.means10$Coef)[2:10],first=10)
clust.coef.means10$Coef<-factor(clust.coef.means10$Coef,levels=clust.coef.means10$Coef)
colnames(clust.coef.means10)[2:11]<-paste0("mean",1:10)

clust.coef.medians10<-cbind(Coef=clust.coef.means10$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential10@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential10@data$classification),median)[,-1])))
clust.coef.sd10<-cbind(Coef=clust.coef.means10$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential10@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential10@data$classification),sd)[,-1])))
colnames(clust.coef.sd10)[2:11]<-paste0("sd",1:10)
colnames(clust.coef.medians10)[2:11]<-paste0("median",1:10)


clust.summary10<-cbind(clust.coef.means10,clust.coef.medians10[,-1],clust.coef.sd10[,-1])

p10<-ggplot(clust.summary10,aes(Coef,mean1))
p10+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=cores[1],size=1,fatten=2)+theme(axis.text.x = element_text(size=15,angle = 45, hjust = 1),axis.text.y = element_text(size=15))+ylab("")+xlab("")+
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=cores[2],size=1,fatten=2,position=position_nudge(x=0.025))+
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=cores[3],size=1,fatten=2,position=position_nudge(x=0.05))+
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=cores[4],size=1,fatten=2,position=position_nudge(x=0.075))+
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=cores[5],size=1,fatten=2,position=position_nudge(x=0.1))+
  geom_pointrange(aes(Coef,mean6,ymin=mean6-sd6,ymax=mean6+sd6),col=cores[6],size=1,fatten=2,position=position_nudge(x=0.125))+
  geom_pointrange(aes(Coef,mean7,ymin=mean7-sd7,ymax=mean7+sd7),col=cores[7],size=1,fatten=2,position=position_nudge(x=0.15))+
  geom_pointrange(aes(Coef,mean8,ymin=mean8-sd8,ymax=mean8+sd8),col=cores[8],size=1,fatten=2,position=position_nudge(x=0.175))+
  geom_pointrange(aes(Coef,mean9,ymin=mean9-sd9,ymax=mean9+sd9),col=cores[9],size=1,fatten=2,position=position_nudge(x=-0.025))+
  geom_pointrange(aes(Coef,mean10,ymin=mean10-sd10,ymax=mean10+sd10),col=cores[10],size=1,fatten=2,position=position_nudge(x=-0.05))


### Limiting number of clusters to 11: ####
mclust_11<- Mclust(data=ggwr_exponential_grid$SDF@data[which(ggwr_exponential_grid$SDF@data$y>0),c(1:15)[-11]],G=1:11) # vegetation type, tree height and topography
summary(mclust_11)

plot(mclust_11, "BIC")

pres.ggwr_exponential11<-ggwr_exponential_grid$SDF[which(ggwr_exponential_grid$SDF@data$y>0),]
pres.ggwr_exponential11$classification<-as.factor(mclust_11$classification)
pres.ggwr_exponential11$uncertainty<-mclust_11$uncertainty
pres.ggwr_exponential11<-spTransform(pres.ggwr_exponential11,CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"))

# with BCRs
spplot(pres.ggwr_exponential11, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=cores[c(1:5,7:12)],
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out,col="gainsboro",fill=colsmapBCRs),basemap)
)

set.seed(432)
brandtcols<-sample(rainbow_hcl(4,alpha=0.8))
# with Brandt boreal
spplot(pres.ggwr_exponential11, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=cores[c(1:5,7:12)],
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out2,col="gainsboro",fill=brandtcols),basemap)
)

# Summary of coefficients per clusters
clust.coef.means11<-as.data.frame(t(aggregate(pres.ggwr_exponential11@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential11@data$classification),mean)[,-1]))
clust.coef.means11<-cbind(Coef=factor(rownames(clust.coef.means11),levels=rownames(clust.coef.means11)),clust.coef.means11)

clust.coef.means11$Coef<-as.character(clust.coef.means11$Coef)
clust.coef.means11$Coef[2:10]<-substring(as.character(clust.coef.means11$Coef)[2:10],first=10)
clust.coef.means11$Coef<-factor(clust.coef.means11$Coef,levels=clust.coef.means11$Coef)
colnames(clust.coef.means11)[2:12]<-paste0("mean",1:11)

clust.coef.medians11<-cbind(Coef=clust.coef.means11$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential11@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential11@data$classification),median)[,-1])))
clust.coef.sd11<-cbind(Coef=clust.coef.means11$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential11@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential11@data$classification),sd)[,-1])))
colnames(clust.coef.sd11)[2:12]<-paste0("sd",1:11)
colnames(clust.coef.medians11)[2:12]<-paste0("median",1:11)


clust.summary11<-cbind(clust.coef.means11,clust.coef.medians11[,-1],clust.coef.sd11[,-1])

p11<-ggplot(clust.summary11,aes(Coef,mean1))
p11+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=cores[1],size=1,fatten=2)+theme(axis.text.x = element_text(size=15,angle = 45, hjust = 1),axis.text.y = element_text(size=15))+ylab("")+xlab("")+
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=cores[2],size=1,fatten=2,position=position_nudge(x=0.025))+
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=cores[3],size=1,fatten=2,position=position_nudge(x=0.05))+
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=cores[4],size=1,fatten=2,position=position_nudge(x=0.075))+
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=cores[5],size=1,fatten=2,position=position_nudge(x=0.1))+
  geom_pointrange(aes(Coef,mean6,ymin=mean6-sd6,ymax=mean6+sd6),col=cores[6],size=1,fatten=2,position=position_nudge(x=0.125))+
  geom_pointrange(aes(Coef,mean7,ymin=mean7-sd7,ymax=mean7+sd7),col=cores[7],size=1,fatten=2,position=position_nudge(x=0.15))+
  geom_pointrange(aes(Coef,mean8,ymin=mean8-sd8,ymax=mean8+sd8),col=cores[8],size=1,fatten=2,position=position_nudge(x=0.175))+
  geom_pointrange(aes(Coef,mean9,ymin=mean9-sd9,ymax=mean9+sd9),col=cores[9],size=1,fatten=2,position=position_nudge(x=-0.025))+
  geom_pointrange(aes(Coef,mean10,ymin=mean10-sd10,ymax=mean10+sd10),col=cores[10],size=1,fatten=2,position=position_nudge(x=-0.05))+
  geom_pointrange(aes(Coef,mean11,ymin=mean11-sd11,ymax=mean11+sd11),col=cores[11],size=1,fatten=2,position=position_nudge(x=-0.05))

### Limiting number of clusters to 12: ####
mclust_12<- Mclust(data=ggwr_exponential_grid$SDF@data[which(ggwr_exponential_grid$SDF@data$y>0),c(1:15)[-11]],G=1:12) # vegetation type, tree height and topography
summary(mclust_12)

plot(mclust_12, "BIC")

pres.ggwr_exponential12<-ggwr_exponential_grid$SDF[which(ggwr_exponential_grid$SDF@data$y>0),]
pres.ggwr_exponential12$classification<-as.factor(mclust_12$classification)
pres.ggwr_exponential12$uncertainty<-mclust_12$uncertainty
pres.ggwr_exponential12<-spTransform(pres.ggwr_exponential12,CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"))

# with BCRs
spplot(pres.ggwr_exponential12, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=cores[1:12],
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out,col="gainsboro",fill=colsmapBCRs),basemap)
)

set.seed(432)
brandtcols<-sample(rainbow_hcl(4,alpha=0.8))
# with Brandt boreal
spplot(pres.ggwr_exponential12, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=cores[1:12],
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out2,col="gainsboro",fill=brandtcols),basemap)
)

# Summary of coefficients per clusters
clust.coef.means12<-as.data.frame(t(aggregate(pres.ggwr_exponential12@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential12@data$classification),mean)[,-1]))
clust.coef.means12<-cbind(Coef=factor(rownames(clust.coef.means12),levels=rownames(clust.coef.means12)),clust.coef.means12)

clust.coef.means12$Coef<-as.character(clust.coef.means12$Coef)
clust.coef.means12$Coef[2:10]<-substring(as.character(clust.coef.means12$Coef)[2:10],first=10)
clust.coef.means12$Coef<-factor(clust.coef.means12$Coef,levels=clust.coef.means12$Coef)
colnames(clust.coef.means12)[2:13]<-paste0("mean",1:12)

clust.coef.medians12<-cbind(Coef=clust.coef.means12$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential12@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential12@data$classification),median)[,-1])))
clust.coef.sd12<-cbind(Coef=clust.coef.means12$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential12@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential12@data$classification),sd)[,-1])))
colnames(clust.coef.sd12)[2:13]<-paste0("sd",1:12)
colnames(clust.coef.medians12)[2:13]<-paste0("median",1:12)


clust.summary12<-cbind(clust.coef.means12,clust.coef.medians12[,-1],clust.coef.sd12[,-1])

p12<-ggplot(clust.summary12,aes(Coef,mean1))
p12+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=cores[1],size=1,fatten=2)+theme(axis.text.x = element_text(size=15,angle = 45, hjust = 1),axis.text.y = element_text(size=15))+ylab("")+xlab("")+
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=cores[2],size=1,fatten=2,position=position_nudge(x=0.025))+
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=cores[3],size=1,fatten=2,position=position_nudge(x=0.05))+
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=cores[4],size=1,fatten=2,position=position_nudge(x=0.075))+
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=cores[5],size=1,fatten=2,position=position_nudge(x=0.1))+
  geom_pointrange(aes(Coef,mean6,ymin=mean6-sd6,ymax=mean6+sd6),col=cores[6],size=1,fatten=2,position=position_nudge(x=0.125))+
  geom_pointrange(aes(Coef,mean7,ymin=mean7-sd7,ymax=mean7+sd7),col=cores[7],size=1,fatten=2,position=position_nudge(x=0.15))+
  geom_pointrange(aes(Coef,mean8,ymin=mean8-sd8,ymax=mean8+sd8),col=cores[8],size=1,fatten=2,position=position_nudge(x=0.175))+
  geom_pointrange(aes(Coef,mean9,ymin=mean9-sd9,ymax=mean9+sd9),col=cores[9],size=1,fatten=2,position=position_nudge(x=-0.025))+
  geom_pointrange(aes(Coef,mean10,ymin=mean10-sd10,ymax=mean10+sd10),col=cores[10],size=1,fatten=2,position=position_nudge(x=-0.05))+
  geom_pointrange(aes(Coef,mean11,ymin=mean11-sd11,ymax=mean11+sd11),col=cores[11],size=1,fatten=2,position=position_nudge(x=-0.05))+
  geom_pointrange(aes(Coef,mean12,ymin=mean12-sd12,ymax=mean12+sd12),col=cores[12],size=1,fatten=2,position=position_nudge(x=-0.05))

### Limiting number of clusters to 13: ####
mclust_13<- Mclust(data=ggwr_exponential_grid$SDF@data[which(ggwr_exponential_grid$SDF@data$y>0),c(1:15)[-11]],G=1:13) # vegetation type, tree height and topography
summary(mclust_13)

plot(mclust_13, "BIC")

pres.ggwr_exponential13<-ggwr_exponential_grid$SDF[which(ggwr_exponential_grid$SDF@data$y>0),]
pres.ggwr_exponential13$classification<-as.factor(mclust_13$classification)
pres.ggwr_exponential13$uncertainty<-mclust_13$uncertainty
pres.ggwr_exponential13<-spTransform(pres.ggwr_exponential13,CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"))

# with BCRs
spplot(pres.ggwr_exponential13, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=cores[1:13],
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out,col="gainsboro",fill=colsmapBCRs),basemap)
)

set.seed(432)
brandtcols<-sample(rainbow_hcl(4,alpha=0.8))
# with Brandt boreal
spplot(pres.ggwr_exponential13, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=cores[1:13],
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out2,col="gainsboro",fill=brandtcols),basemap)
)

# Summary of coefficients per clusters
clust.coef.means13<-as.data.frame(t(aggregate(pres.ggwr_exponential13@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential13@data$classification),mean)[,-1]))
clust.coef.means13<-cbind(Coef=factor(rownames(clust.coef.means13),levels=rownames(clust.coef.means13)),clust.coef.means13)

clust.coef.means13$Coef<-as.character(clust.coef.means13$Coef)
clust.coef.means13$Coef[2:10]<-substring(as.character(clust.coef.means13$Coef)[2:10],first=10)
clust.coef.means13$Coef<-factor(clust.coef.means13$Coef,levels=clust.coef.means13$Coef)
colnames(clust.coef.means13)[2:14]<-paste0("mean",1:13)

clust.coef.medians13<-cbind(Coef=clust.coef.means13$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential13@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential13@data$classification),median)[,-1])))
clust.coef.sd13<-cbind(Coef=clust.coef.means13$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential13@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential13@data$classification),sd)[,-1])))
colnames(clust.coef.sd13)[2:14]<-paste0("sd",1:13)
colnames(clust.coef.medians13)[2:14]<-paste0("median",1:13)


clust.summary13<-cbind(clust.coef.means13,clust.coef.medians13[,-1],clust.coef.sd13[,-1])

p13<-ggplot(clust.summary13,aes(Coef,mean1))
p13+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=cores[1],size=1,fatten=2)+theme(axis.text.x = element_text(size=15,angle = 45, hjust = 1),axis.text.y = element_text(size=15))+ylab("")+xlab("")+
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=cores[2],size=1,fatten=2,position=position_nudge(x=0.025))+
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=cores[3],size=1,fatten=2,position=position_nudge(x=0.05))+
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=cores[4],size=1,fatten=2,position=position_nudge(x=0.075))+
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=cores[5],size=1,fatten=2,position=position_nudge(x=0.1))+
  geom_pointrange(aes(Coef,mean6,ymin=mean6-sd6,ymax=mean6+sd6),col=cores[6],size=1,fatten=2,position=position_nudge(x=0.125))+
  geom_pointrange(aes(Coef,mean7,ymin=mean7-sd7,ymax=mean7+sd7),col=cores[7],size=1,fatten=2,position=position_nudge(x=0.15))+
  geom_pointrange(aes(Coef,mean8,ymin=mean8-sd8,ymax=mean8+sd8),col=cores[8],size=1,fatten=2,position=position_nudge(x=0.175))+
  geom_pointrange(aes(Coef,mean9,ymin=mean9-sd9,ymax=mean9+sd9),col=cores[9],size=1,fatten=2,position=position_nudge(x=-0.025))+
  geom_pointrange(aes(Coef,mean10,ymin=mean10-sd10,ymax=mean10+sd10),col=cores[10],size=1,fatten=2,position=position_nudge(x=-0.05))+
  geom_pointrange(aes(Coef,mean11,ymin=mean11-sd11,ymax=mean11+sd11),col=cores[11],size=1,fatten=2,position=position_nudge(x=-0.05))+
  geom_pointrange(aes(Coef,mean12,ymin=mean12-sd12,ymax=mean12+sd12),col=cores[12],size=1,fatten=2,position=position_nudge(x=-0.05))+
  geom_pointrange(aes(Coef,mean13,ymin=mean13-sd13,ymax=mean13+sd13),col=cores[13],size=1,fatten=2,position=position_nudge(x=-0.05))


### Limiting number of clusters to 14: ####
mclust_14<- Mclust(data=ggwr_exponential_grid$SDF@data[which(ggwr_exponential_grid$SDF@data$y>0),c(1:15)[-11]],G=1:14) # vegetation type, tree height and topography
summary(mclust_14)

plot(mclust_14, "BIC")

pres.ggwr_exponential14<-ggwr_exponential_grid$SDF[which(ggwr_exponential_grid$SDF@data$y>0),]
pres.ggwr_exponential14$classification<-as.factor(mclust_14$classification)
pres.ggwr_exponential14$uncertainty<-mclust_14$uncertainty
pres.ggwr_exponential14<-spTransform(pres.ggwr_exponential14,CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"))

# with BCRs
spplot(pres.ggwr_exponential14, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=cores[1:14],
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out,col="gainsboro",fill=colsmapBCRs),basemap)
)

set.seed(432)
brandtcols<-sample(rainbow_hcl(4,alpha=0.8))
# with Brandt boreal
spplot(pres.ggwr_exponential14, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="Clustering by vegetation type, tree height, and topography",
       col.regions=cores[1:14],
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out2,col="gainsboro",fill=brandtcols),basemap)
)

# Summary of coefficients per clusters
clust.coef.means14<-as.data.frame(t(aggregate(pres.ggwr_exponential14@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential14@data$classification),mean)[,-1]))
clust.coef.means14<-cbind(Coef=factor(rownames(clust.coef.means14),levels=rownames(clust.coef.means14)),clust.coef.means14)

clust.coef.means14$Coef<-as.character(clust.coef.means14$Coef)
clust.coef.means14$Coef[2:10]<-substring(as.character(clust.coef.means14$Coef)[2:10],first=10)
clust.coef.means14$Coef<-factor(clust.coef.means14$Coef,levels=clust.coef.means14$Coef)
colnames(clust.coef.means14)[2:15]<-paste0("mean",1:143)

clust.coef.medians14<-cbind(Coef=clust.coef.means14$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential14@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential14@data$classification),median)[,-1])))
clust.coef.sd14<-cbind(Coef=clust.coef.means14$Coef,as.data.frame(t(aggregate(pres.ggwr_exponential14@data[,c(1:15)[-11]],list(group=pres.ggwr_exponential14@data$classification),sd)[,-1])))
colnames(clust.coef.sd14)[2:15]<-paste0("sd",1:14)
colnames(clust.coef.medians14)[2:15]<-paste0("median",1:14)


clust.summary14<-cbind(clust.coef.means14,clust.coef.medians14[,-1],clust.coef.sd14[,-1])

p14<-ggplot(clust.summary14,aes(Coef,mean1))
p14+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=cores[1],size=1,fatten=2)+theme(axis.text.x = element_text(size=15,angle = 45, hjust = 1),axis.text.y = element_text(size=15))+ylab("")+xlab("")+
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=cores[2],size=1,fatten=2,position=position_nudge(x=0.025))+
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=cores[3],size=1,fatten=2,position=position_nudge(x=0.05))+
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=cores[4],size=1,fatten=2,position=position_nudge(x=0.075))+
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=cores[5],size=1,fatten=2,position=position_nudge(x=0.1))+
  geom_pointrange(aes(Coef,mean6,ymin=mean6-sd6,ymax=mean6+sd6),col=cores[6],size=1,fatten=2,position=position_nudge(x=0.125))+
  geom_pointrange(aes(Coef,mean7,ymin=mean7-sd7,ymax=mean7+sd7),col=cores[7],size=1,fatten=2,position=position_nudge(x=0.15))+
  geom_pointrange(aes(Coef,mean8,ymin=mean8-sd8,ymax=mean8+sd8),col=cores[8],size=1,fatten=2,position=position_nudge(x=0.175))+
  geom_pointrange(aes(Coef,mean9,ymin=mean9-sd9,ymax=mean9+sd9),col=cores[9],size=1,fatten=2,position=position_nudge(x=-0.025))+
  geom_pointrange(aes(Coef,mean10,ymin=mean10-sd10,ymax=mean10+sd10),col=cores[10],size=1,fatten=2,position=position_nudge(x=-0.05))+
  geom_pointrange(aes(Coef,mean11,ymin=mean11-sd11,ymax=mean11+sd11),col=cores[11],size=1,fatten=2,position=position_nudge(x=-0.05))+
  geom_pointrange(aes(Coef,mean12,ymin=mean12-sd12,ymax=mean12+sd12),col=cores[12],size=1,fatten=2,position=position_nudge(x=-0.05))+
  geom_pointrange(aes(Coef,mean13,ymin=mean13-sd13,ymax=mean13+sd13),col=cores[13],size=1,fatten=2,position=position_nudge(x=-0.05))+
  geom_pointrange(aes(Coef,mean14,ymin=mean14-sd14,ymax=mean14+sd14),col=cores[14],size=1,fatten=2,position=position_nudge(x=-0.05))




#### Re-run model for additional samples from BAM dataset, but only on the group 1 data from 3-cluster model (everything east of mid-Ontario). ####
# clear some of the workspace
rm(list=ls()[! ls() %in% c("cawa_mm","thin_cawa_mm","basemap","BCRs","brandt","plot_sampled","sample_quad","ggwr_exponential_grid","out","out2","pres.ggwr_exponential3")])  


# define longitude that defines the western limit of the eastern cluster
Ylim<-min(pres.ggwr_exponential3@coords[ which(pres.ggwr_exponential3$classification==1),1])

cawa_mm_east<-cawa_mm[which(cawa_mm$Xcl>=Ylim),]

samples_east<-as.list(rep(NA,10))

names(samples_east)<-paste0("sample_east",1:10)

sample_union<-function(seed=123){
  s1<-sample_quad(cawa_mm_east,Xintervals=60,Yintervals=50,seed=seed)
  s2<-sample_quad(cawa_mm_east,Xintervals=60,Yintervals=50,onlypres = T, N=4,seed=seed)
  union<-union(s1,s2)
  union
}


samples_east$sample_east1<- sample_union(seed=1)
samples_east$sample_east2<- sample_union(seed=2)
samples_east$sample_east3<- sample_union(seed=3)
samples_east$sample_east4<- sample_union(seed=4)
samples_east$sample_east5<- sample_union(seed=5)
samples_east$sample_east6<- sample_union(seed=6)
samples_east$sample_east7<- sample_union(seed=7)
samples_east$sample_east8<- sample_union(seed=8)
samples_east$sample_east9<- sample_union(seed=9)
samples_east$sample_east10<- sample_union(seed=10)


# define datasets:
thin_cawa_mm_east <-lapply(samples_east,function(x){cawa_mm_east[x,]})

# sample-specific landcover class optimization
lc_optim<-function(thin_cawa_mm){
  ol<- optilevels(y=thin_cawa_mm$count, x=thin_cawa_mm$HAB_NALC1, dist="poisson", offset=thin_cawa_mm$offset)
  thin_cawa_mm$HAB_NALC2<-thin_cawa_mm$HAB_NALC1
  levels(thin_cawa_mm$HAB_NALC2)<-ol$level[[length(ol$level)]]
  return(thin_cawa_mm)
}

#thin_cawa_mm_east<-lapply(thin_cawa_mm_east,lc_optim)

# landcover class optimization v.2: use union of samples and one optimization for all
union_samples<-union(union(union(union(union(union(union(union(union(samples_east$sample_east1,samples_east$sample2),samples_east$sample3),samples_east$sample4),samples_east$sample5),samples_east$sample6),samples_east$sample7),samples_east$sample8),samples_east$sample9),samples_east$sample10)
cawa_lc_optim<-cawa_mm_east[union_samples,]



lc_optim2<-function(thin_cawa_mm){
  ol<- optilevels(y=cawa_lc_optim$count, x=cawa_lc_optim$HAB_NALC1, dist="poisson", offset=cawa_lc_optim$offset)
  thin_cawa_mm$HAB_NALC2<-thin_cawa_mm$HAB_NALC1
  levels(thin_cawa_mm$HAB_NALC2)<-ol$level[[length(ol$level)]]
  return(thin_cawa_mm)
}
thin_cawa_mm_east<-lapply(thin_cawa_mm_east,lc_optim2)

# centering HGT, HGT2, CTI and CIT2 variables for each sample based on their respective means from the union of all samples

centercovs<-function(sample){
  sample$centerHGT<-sample$HGT-mean(cawa_lc_optim$HGT)
  sample$centerHGT2<-sample$HGT2-mean(cawa_lc_optim$HGT2)
  sample$centerCTI<-sample$CTI-mean(cawa_lc_optim$CTI)
  sample$centerCTI2<-sample$CTI2-mean(cawa_lc_optim$CTI2)
  return(sample)
}

thin_cawa_mm_east<-lapply(thin_cawa_mm_east,centercovs)


# create spatial dataframe from sampled dataset

mmsp_fun<-function(thin_cawa_mm){
  mmsp<- SpatialPointsDataFrame(coords=cbind(thin_cawa_mm$X,thin_cawa_mm$Y),proj4string = CRS("+proj=longlat +ellps=WGS84"),data=thin_cawa_mm)
}


mmsp_east <- lapply(thin_cawa_mm_east,mmsp_fun)

#### GW modeling and clustering - East samples #####


## Specify distance matrix:##

ggwr_east<-function(index,modelname,expandlcc=FALSE){
  DM <- gw.dist(dp.locat = coordinates(mmsp_east[[index]]), longlat=TRUE)
  
  if(expandlcc==TRUE){ # using cell means model parameterization
    bw.ggwr.exponential_grid_east <- bw.ggwr(count ~  (HAB_NALC2-1) + ROAD + centerHGT + centerHGT2 + centerCTI + centerCTI2 + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 + CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT + offset(mmsp_east[[index]]$offset), data = mmsp_east[[index]], approach = "AICc", kernel = "exponential", adaptive = TRUE, family="poisson", longlat = TRUE, dMat=DM) 
    
    # Fit model
    ggwr_exponential_grid_east<-ggwr.basic(count ~ (HAB_NALC2-1) + ROAD + centerHGT + centerHGT2 + centerCTI + centerCTI2 + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 + CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT + offset(mmsp_east[[index]]$offset), data = mmsp_east[[index]], bw = bw.ggwr.exponential_grid_east, kernel = "exponential", adaptive = TRUE, longlat=TRUE, family="poisson", dMat=DM)
  }
  
  else{ # using treatment contrasts parameterization
    bw.ggwr.exponential_grid_east <- bw.ggwr(count ~  HAB_NALC2 + ROAD + centerHGT + centerHGT2 + centerCTI + centerCTI2 + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 + CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT + offset(mmsp_east[[index]]$offset), data = mmsp_east[[index]], approach = "AICc", kernel = "exponential", adaptive = TRUE, family="poisson", longlat = TRUE, dMat=DM) 
    
    # Fit model
    ggwr_exponential_grid_east<-ggwr.basic(count ~ HAB_NALC2 + ROAD + centerHGT + centerHGT2 + centerCTI + centerCTI2 + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 + CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT + offset(mmsp_east[[index]]$offset), data = mmsp_east[[index]], bw = bw.ggwr.exponential_grid_east, kernel = "exponential", adaptive = TRUE, longlat=TRUE, family="poisson", dMat=DM)
  }

  assign(modelname,ggwr_exponential_grid_east)
  
  save(list=c(modelname),file=paste0("D:/CHID subunit delineation/output/",paste0(modelname,".Rdata")))
  
  return(get(modelname))
} 



save.image("D:/CHID subunit delineation/subunits.RData")

ggwr_east_1<-ggwr_east(index=1,modelname="ggwr_east_1",expandlcc = T)
ggwr_east_3<-ggwr_east(index=3,modelname="ggwr_east_3",expandlcc = T)
ggwr_east_5<-ggwr_east(index=5,modelname="ggwr_east_5",expandlcc = T)
ggwr_east_7<-ggwr_east(index=7,modelname="ggwr_east_7",expandlcc = T)
ggwr_east_9<-ggwr_east(index=9,modelname="ggwr_east_9",expandlcc = T)

ggwr_east_2<-ggwr_east(index=2,modelname="ggwr_east_2",expandlcc = T)
ggwr_east_4<-ggwr_east(index=4,modelname="ggwr_east_4",expandlcc = T)
ggwr_east_6<-ggwr_east(index=6,modelname="ggwr_east_6",expandlcc = T)
ggwr_east_8<-ggwr_east(index=8,modelname="ggwr_east_8",expandlcc = T)
ggwr_east_10<-ggwr_east(index=10,modelname="ggwr_east_10",expandlcc = T)

#### Generate east clusters ####
load("D:/CHID subunit delineation/output/ggwr_east_1.Rdata")
load("D:/CHID subunit delineation/output/ggwr_east_2.Rdata")
load("D:/CHID subunit delineation/output/ggwr_east_3.Rdata")
load("D:/CHID subunit delineation/output/ggwr_east_4.Rdata")
load("D:/CHID subunit delineation/output/ggwr_east_5.Rdata")
load("D:/CHID subunit delineation/output/ggwr_east_6.Rdata")
load("D:/CHID subunit delineation/output/ggwr_east_7.Rdata")
load("D:/CHID subunit delineation/output/ggwr_east_8.Rdata")
load("D:/CHID subunit delineation/output/ggwr_east_9.Rdata")
load("D:/CHID subunit delineation/output/ggwr_east_10.Rdata")



# clear some of the workspace
rm(list=ls()[! ls() %in% c("cawa_mm","thin_cawa_mm","basemap","BCRs","brandt","plot_sampled","sample_quad","ggwr_exponential_grid","out","out2","out3","out4","ggwr_east","ggwr_east_1","ggwr_east_10","ggwr_east_2","ggwr_east_3","ggwr_east_4","ggwr_east_5","ggwr_east_6","ggwr_east_7","ggwr_east_8","ggwr_east_9","pres.ggwr_exponential3")])  

# a function for clustering of each model using vegetation type, tree height and topography coefficients
gw_clust<-function(gwmodel,maxclusters=5, out3=NULL,out4=NULL, return.outs=FALSE){
  
  colmax<-which(colnames(gwmodel$SDF@data)=="centerCTI2")
  roadcol<-which(colnames(gwmodel$SDF@data)=="ROAD")
  
  set.seed(123)
  ##Clustering and BIC plot
  mclust<- Mclust(data=gwmodel$SDF@data[which(gwmodel$SDF@data$y>0),c(1:colmax)[-roadcol]],G=1:maxclusters,modelNames="VVV")  
  
  
  ##Plot mapped clusters
  # Create spatial dataframe for points to be plotted
  pres.ggwr<-gwmodel$SDF[which(gwmodel$SDF@data$y>0),]
  pres.ggwr$classification<-as.factor(mclust$classification)
  pres.ggwr$uncertainty<-mclust$uncertainty
  # Reproject
  pres.ggwr<-spTransform(pres.ggwr,CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"))
  # Create simple plot of points to obtain plot extents, and crop BCR layer to that extent
  p<- spplot(pres.ggwr, "classification", do.log = F)
  p.extent<-extent(rbind(p$x.limits,p$y.limits))
  
  if(is.null(out3))  {
    out3<-crop(BCRs,p.extent)
    out3$BCR<-as.factor(out3$BCR)
    levels(out3$BCR)<-1:11
    }
  
  
  
  set.seed(12)
  colsmapBCRs<-c("#ffffff",sample(rainbow_hcl(11,alpha=0.8))[-1])[out3$BCR]

  # with BCRs
  BCRplot<- spplot(pres.ggwr, "classification", do.log = F,
         key.space=list(x=0.8,y=0.3,corner=c(0,1)),
         main="",
         col.regions=brewer.pal(maxclusters,"Set1"),
         #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
         sp.layout=list(list("sp.polygons",out3,col="gainsboro",fill=colsmapBCRs),basemap))
  
  set.seed(432)
  
  if(is.null(out4)){out4<-crop(brandt,p.extent)}
  
  brandtcols<-sample(rainbow_hcl(4,alpha=0.8))
  # with Brandt boreal
  Brandtplot<-spplot(pres.ggwr, "classification", do.log = F,
         key.space=list(x=0.8,y=0.3,corner=c(0,1)),
         main="",
         col.regions=brewer.pal(maxclusters,"Set1"),
         #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
         sp.layout=list(list("sp.polygons",out4,col="gainsboro",fill=brandtcols),basemap))
  
  # Summary of coefficients per clusters
  clust.coef.means<-as.data.frame(t(aggregate(pres.ggwr@data[,c(1:colmax)[-roadcol]],list(group=pres.ggwr@data$classification),mean)[,-1]))
  clust.coef.means<-cbind(Coef=factor(rownames(clust.coef.means),levels=rownames(clust.coef.means)),clust.coef.means)
  
  clust.coef.means$Coef<-as.character(clust.coef.means$Coef)
  clust.coef.means$Coef[1:(roadcol-1)]<-substring(as.character(clust.coef.means$Coef)[1:(roadcol-1)],first=10)
  clust.coef.means$Coef<-factor(clust.coef.means$Coef,levels=clust.coef.means$Coef)
  colnames(clust.coef.means)[-1]<-paste0("mean",1:maxclusters)
  
  clust.coef.sd<-cbind(Coef=clust.coef.means$Coef,as.data.frame(t(aggregate(pres.ggwr@data[,c(1:colmax)[-11]],list(group=pres.ggwr@data$classification),sd)[,-1])))
  colnames(clust.coef.sd)[-1]<-paste0("sd",1:maxclusters)
  
  
  clust.summary<-cbind(clust.coef.means,clust.coef.sd[,-1])
  
  # if maxclusters is different from 5, need to adjust the plot below
  coefplot<-ggplot(clust.summary,aes(Coef,mean1))
  
  if(return.outs==FALSE){out<-list(mclust=mclust,pres.ggwr=pres.ggwr,BCRplot=BCRplot,Brandtplot=Brandtplot,coefplot=coefplot)}
  else {out<-list(mclust=mclust,pres.ggwr=pres.ggwr,BCRplot=BCRplot,Brandtplot=Brandtplot,coefplot=coefplot,out3=out3,out4=out4)}
  
  return(out)
}

### East sample1 - 3 clusters ----
result.3clust.sample1<-gw_clust(ggwr_east_1,maxclusters=3,return.outs=T)
#plot(result.3clust.sample1$mclust,"BIC")
#result.3clust.sample1$BCRplot
#result.3clust.sample1$Brandtplot
p.3.1<-result.3clust.sample1$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0.1))+theme(axis.text.x = element_text(size=9,angle = 45, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.1))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0)) #green


### East sample2 - 3 clusters ----
result.3clust.sample2<-gw_clust(ggwr_east_2,maxclusters=3,out3=result.3clust.sample1$out3,out4=result.3clust.sample1$out4)
#plot(result.3clust.sample2$mclust,"BIC")
#result.3clust.sample2$BCRplot
#result.3clust.sample2$Brandtplot
p.3.2<-result.3clust.sample2$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0))+theme(axis.text.x = element_text(size=9,angle = 45, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=0.1))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=-0.1)) #green

### East sample3 - 3 clusters ----
result.3clust.sample3<-gw_clust(ggwr_east_3,maxclusters=3,out3=result.3clust.sample1$out3,out4=result.3clust.sample1$out4)
#plot(result.3clust.sample3$mclust,"BIC")
#result.3clust.sample3$BCRplot
#result.3clust.sample3$Brandtplot
p.3.3<-result.3clust.sample3$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0.1))+theme(axis.text.x = element_text(size=9,angle = 60, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.1))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0)) #green

### East sample4 - 3 clusters ----
result.3clust.sample4<-gw_clust(ggwr_east_4,maxclusters=3,out3=result.3clust.sample1$out3,out4=result.3clust.sample1$out4)
#plot(result.3clust.sample4$mclust,"BIC")
#result.3clust.sample4$BCRplot
#result.3clust.sample4$Brandtplot
p.3.4<-result.3clust.sample4$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0.1))+theme(axis.text.x = element_text(size=9,angle = 60, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.1))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0)) #green

### East sample5 - 3 clusters ----
result.3clust.sample5<-gw_clust(ggwr_east_5,maxclusters=3,out3=result.3clust.sample1$out3,out4=result.3clust.sample1$out4)
#plot(result.3clust.sample5$mclust,"BIC")
#result.3clust.sample5$BCRplot
#result.3clust.sample5$Brandtplot
p.3.5<-result.3clust.sample5$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0.1))+theme(axis.text.x = element_text(size=9,angle = 60, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=0))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=-0.1)) #green

### East sample6 - 3 clusters ----
result.3clust.sample6<-gw_clust(ggwr_east_6,maxclusters=3,out3=result.3clust.sample1$out3,out4=result.3clust.sample1$out4)
#plot(result.3clust.sample6$mclust,"BIC")
#result.3clust.sample6$BCRplot
#result.3clust.sample6$Brandtplot
p.3.6<-result.3clust.sample6$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0.1))+theme(axis.text.x = element_text(size=9,angle = 60, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=0))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=-0.1)) #green

### East sample7 - 3 clusters ----
result.3clust.sample7<-gw_clust(ggwr_east_7,maxclusters=3,out3=result.3clust.sample1$out3,out4=result.3clust.sample1$out4)
#plot(result.3clust.sample7$mclust,"BIC")
#result.3clust.sample7$BCRplot
#result.3clust.sample7$Brandtplot
p.3.7<-result.3clust.sample7$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0))+theme(axis.text.x = element_text(size=9,angle = 60, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=0.1))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=-0.1)) #green


### East sample8 - 3 clusters ----
result.3clust.sample8<-gw_clust(ggwr_east_8,maxclusters=3,out3=result.3clust.sample1$out3,out4=result.3clust.sample1$out4)
#plot(result.3clust.sample8$mclust,"BIC")
#result.3clust.sample8$BCRplot
#result.3clust.sample8$Brandtplot
p.3.8<-result.3clust.sample8$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0))+theme(axis.text.x = element_text(size=9,angle = 60, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.1))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0.1)) #green

### East sample9 - 3 clusters ----
result.3clust.sample9<-gw_clust(ggwr_east_9,maxclusters=3,out3=result.3clust.sample1$out3,out4=result.3clust.sample1$out4)
#plot(result.3clust.sample9$mclust,"BIC")
#result.3clust.sample9$BCRplot
#result.3clust.sample9$Brandtplot
p.3.9<-result.3clust.sample9$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2,position=position_nudge(x=-0.1))+theme(axis.text.x = element_text(size=9,angle = 60, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=0.1))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0)) #green

### East sample10 - 3 clusters ----
result.3clust.sample10<-gw_clust(ggwr_east_10,maxclusters=3,out3=result.3clust.sample1$out3,out4=result.3clust.sample1$out4)
#plot(result.3clust.sample10$mclust,"BIC")
#result.3clust.sample10$BCRplot
#result.3clust.sample10$Brandtplot
p.3.10<-result.3clust.sample10$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0))+theme(axis.text.x = element_text(size=9,angle = 60, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.1))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0.1)) #green



### East sample1 - 5 clusters ----
result.5clust.sample1<-gw_clust(ggwr_east_1,out3=result.3clust.sample1$out3,out4=result.3clust.sample1$out4,return.outs=T)
#plot(result.5clust.sample1$mclust,"BIC")
#result.5clust.sample1$BCRplot
#result.5clust.sample1$Brandtplot
p.5.1<-result.5clust.sample1$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0.2))+theme(axis.text.x = element_text(size=9,angle = 45, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.1))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(5,"Set1")[4],size=1,fatten=2,position=position_nudge(x=-0.2))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(5,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.1)) #orange

### East sample2 - 5 clusters ----
result.5clust.sample2<-gw_clust(ggwr_east_2,out3=result.5clust.sample1$out3,out4=result.5clust.sample1$out4)
#plot(result.5clust.sample2$mclust,"BIC")
#result.5clust.sample2$BCRplot
#result.5clust.sample2$Brandtplot
p.5.2<-result.5clust.sample2$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2,position=position_nudge(x=-0.1))+theme(axis.text.x = element_text(size=9,angle = 45, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=0))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=-0.2))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(5,"Set1")[4],size=1,fatten=2,position=position_nudge(x=0.1))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(5,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.2)) #orange

### East sample3 - 5 clusters ----
result.5clust.sample3<-gw_clust(ggwr_east_3,out3=result.5clust.sample1$out3,out4=result.5clust.sample1$out4)
#plot(result.5clust.sample3$mclust,"BIC")
#result.5clust.sample3$BCRplot
#result.5clust.sample3$Brandtplot
p.5.3<-result.5clust.sample3$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2,position=position_nudge(x=-0.1))+theme(axis.text.x = element_text(size=9,angle = 60, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=0))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=-0.2))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(5,"Set1")[4],size=1,fatten=2,position=position_nudge(x=0.1))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(5,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.2)) #orange

### East sample4 - 5 clusters ----
result.5clust.sample4<-gw_clust(ggwr_east_4,out3=result.5clust.sample1$out3,out4=result.5clust.sample1$out4)
#plot(result.5clust.sample4$mclust,"BIC")
#result.5clust.sample4$BCRplot
#result.5clust.sample4$Brandtplot
p.5.4<-result.5clust.sample4$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2,position=position_nudge(x=-0.1))+theme(axis.text.x = element_text(size=9,angle = 60, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=0.1))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=-0.2))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(5,"Set1")[4],size=1,fatten=2,position=position_nudge(x=0))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(5,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.2)) #orange

### East sample5 - 5 clusters ----
result.5clust.sample5<-gw_clust(ggwr_east_5,out3=result.5clust.sample1$out3,out4=result.5clust.sample1$out4)
#plot(result.5clust.sample5$mclust,"BIC")
#result.5clust.sample5$BCRplot
#result.5clust.sample5$Brandtplot
p.5.5<-result.5clust.sample5$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0.1))+theme(axis.text.x = element_text(size=9,angle = 60, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.1))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=-0.2))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(5,"Set1")[4],size=1,fatten=2,position=position_nudge(x=0))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(5,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.2)) #orange

### East sample6 - 5 clusters ----
result.5clust.sample6<-gw_clust(ggwr_east_6,out3=result.5clust.sample1$out3,out4=result.5clust.sample1$out4)
#plot(result.5clust.sample6$mclust,"BIC")
#result.5clust.sample6$BCRplot
#result.5clust.sample6$Brandtplot
p.5.6<-result.5clust.sample6$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0.2))+theme(axis.text.x = element_text(size=9,angle = 60, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.1))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(5,"Set1")[4],size=1,fatten=2,position=position_nudge(x=-0.1))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(5,"Set1")[5],size=1,fatten=2,position=position_nudge(x=-0.2)) #orange

### East sample7 - 5 clusters ----
result.5clust.sample7<-gw_clust(ggwr_east_7,out3=result.5clust.sample1$out3,out4=result.5clust.sample1$out4)
#plot(result.5clust.sample7$mclust,"BIC")
#result.5clust.sample7$BCRplot
#result.5clust.sample7$Brandtplot
p.5.7<-result.5clust.sample7$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2,position=position_nudge(x=-0.1))+theme(axis.text.x = element_text(size=9,angle = 60, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=0.2))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=-0.2))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(5,"Set1")[4],size=1,fatten=2,position=position_nudge(x=0))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(5,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.1)) #orange


### East sample8 - 5 clusters ----
result.5clust.sample8<-gw_clust(ggwr_east_8,out3=result.5clust.sample1$out3,out4=result.5clust.sample1$out4)
#plot(result.5clust.sample8$mclust,"BIC")
#result.5clust.sample8$BCRplot
#result.5clust.sample8$Brandtplot
p.5.8<-result.5clust.sample8$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0))+theme(axis.text.x = element_text(size=9,angle = 60, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=0.2))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0.1))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(5,"Set1")[4],size=1,fatten=2,position=position_nudge(x=-0.2))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(5,"Set1")[5],size=1,fatten=2,position=position_nudge(x=-0.1)) #orange

### East sample9 - 5 clusters ----
result.5clust.sample9<-gw_clust(ggwr_east_9,out3=result.5clust.sample1$out3,out4=result.5clust.sample1$out4)
#plot(result.5clust.sample9$mclust,"BIC")
#result.5clust.sample9$BCRplot
#result.5clust.sample9$Brandtplot
p.5.9<-result.5clust.sample9$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2,position=position_nudge(x=-0.1))+theme(axis.text.x = element_text(size=9,angle = 60, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=0.1))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0.2))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(5,"Set1")[4],size=1,fatten=2,position=position_nudge(x=-0.2))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(5,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0)) #orange

### East sample10 - 5 clusters ----
result.5clust.sample10<-gw_clust(ggwr_east_10,out3=result.5clust.sample1$out3,out4=result.5clust.sample1$out4)
#plot(result.5clust.sample10$mclust,"BIC")
#result.5clust.sample10$BCRplot
#result.5clust.sample10$Brandtplot
p.5.10<-result.5clust.sample10$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(5,"Set1")[1],size=1,fatten=2,position=position_nudge(x=-0.1))+theme(axis.text.x = element_text(size=9,angle = 60, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(5,"Set1")[2],size=1,fatten=2,position=position_nudge(x=0.1))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(5,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(5,"Set1")[4],size=1,fatten=2,position=position_nudge(x=-0.2))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(5,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.2)) #orange


### East sample1 - 7 clusters ----
result.7clust.sample1<-gw_clust(ggwr_east_1,maxclusters=7,out3=result.3clust.sample1$out3,out4=result.3clust.sample1$out4,return.outs=T)
#plot(result.7clust.sample1$mclust,"BIC")
#result.7clust.sample1$BCRplot
#result.7clust.sample1$Brandtplot
p.7.1<-result.7clust.sample1$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(7,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0.15))+theme(axis.text.x = element_text(size=9,angle = 45, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(7,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.05))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(7,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(7,"Set1")[4],size=1,fatten=2,position=position_nudge(x=-0.15))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(7,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.05))+#orange
  geom_pointrange(aes(Coef,mean6,ymin=mean6-sd6,ymax=mean6+sd6),col=brewer.pal(7,"Set1")[6],size=1,fatten=2,position=position_nudge(x=-0.1))+ #yellow
  geom_pointrange(aes(Coef,mean7,ymin=mean7-sd7,ymax=mean7+sd7),col=brewer.pal(7,"Set1")[7],size=1,fatten=2,position=position_nudge(x=0.1))#brown

### East sample2 - 7 clusters ----
result.7clust.sample2<-gw_clust(ggwr_east_2,maxclusters=7,out3=result.7clust.sample1$out3,out4=result.7clust.sample1$out4)
#plot(result.7clust.sample2$mclust,"BIC")
#result.7clust.sample2$BCRplot
#result.7clust.sample2$Brandtplot
p.7.2<-result.7clust.sample2$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(7,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0.15))+theme(axis.text.x = element_text(size=9,angle = 45, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(7,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.05))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(7,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(7,"Set1")[4],size=1,fatten=2,position=position_nudge(x=-0.15))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(7,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.05))+#orange
  geom_pointrange(aes(Coef,mean6,ymin=mean6-sd6,ymax=mean6+sd6),col=brewer.pal(7,"Set1")[6],size=1,fatten=2,position=position_nudge(x=-0.1))+ #yellow
  geom_pointrange(aes(Coef,mean7,ymin=mean7-sd7,ymax=mean7+sd7),col=brewer.pal(7,"Set1")[7],size=1,fatten=2,position=position_nudge(x=0.1))#brown

### East sample3 - 7 clusters ----
result.7clust.sample3<-gw_clust(ggwr_east_3,maxclusters=7,out3=result.7clust.sample1$out3,out4=result.7clust.sample1$out4)
#plot(result.7clust.sample3$mclust,"BIC")
#result.7clust.sample3$BCRplot
#result.7clust.sample3$Brandtplot
p.7.3<-result.7clust.sample3$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(7,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0.15))+theme(axis.text.x = element_text(size=9,angle = 45, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(7,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.05))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(7,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(7,"Set1")[4],size=1,fatten=2,position=position_nudge(x=-0.15))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(7,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.05))+#orange
  geom_pointrange(aes(Coef,mean6,ymin=mean6-sd6,ymax=mean6+sd6),col=brewer.pal(7,"Set1")[6],size=1,fatten=2,position=position_nudge(x=-0.1))+ #yellow
  geom_pointrange(aes(Coef,mean7,ymin=mean7-sd7,ymax=mean7+sd7),col=brewer.pal(7,"Set1")[7],size=1,fatten=2,position=position_nudge(x=0.1))#brown

### East sample4 - 7 clusters ----
result.7clust.sample4<-gw_clust(ggwr_east_4,maxclusters=7,out3=result.7clust.sample1$out3,out4=result.7clust.sample1$out4)
#plot(result.7clust.sample4$mclust,"BIC")
#result.7clust.sample4$BCRplot
#result.7clust.sample4$Brandtplot
p.7.4<-result.7clust.sample4$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(7,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0.15))+theme(axis.text.x = element_text(size=9,angle = 45, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(7,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.05))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(7,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(7,"Set1")[4],size=1,fatten=2,position=position_nudge(x=-0.15))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(7,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.05))+#orange
  geom_pointrange(aes(Coef,mean6,ymin=mean6-sd6,ymax=mean6+sd6),col=brewer.pal(7,"Set1")[6],size=1,fatten=2,position=position_nudge(x=-0.1))+ #yellow
  geom_pointrange(aes(Coef,mean7,ymin=mean7-sd7,ymax=mean7+sd7),col=brewer.pal(7,"Set1")[7],size=1,fatten=2,position=position_nudge(x=0.1))#brown


### East sample5 - 7 clusters ----
result.7clust.sample5<-gw_clust(ggwr_east_5,maxclusters=7,out3=result.7clust.sample1$out3,out4=result.7clust.sample1$out4)
#plot(result.7clust.sample5$mclust,"BIC")
#result.7clust.sample5$BCRplot
#result.7clust.sample5$Brandtplot
p.7.5<-result.7clust.sample5$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(7,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0.15))+theme(axis.text.x = element_text(size=9,angle = 45, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(7,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.05))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(7,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(7,"Set1")[4],size=1,fatten=2,position=position_nudge(x=-0.15))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(7,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.05))+#orange
  geom_pointrange(aes(Coef,mean6,ymin=mean6-sd6,ymax=mean6+sd6),col=brewer.pal(7,"Set1")[6],size=1,fatten=2,position=position_nudge(x=-0.1))+ #yellow
  geom_pointrange(aes(Coef,mean7,ymin=mean7-sd7,ymax=mean7+sd7),col=brewer.pal(7,"Set1")[7],size=1,fatten=2,position=position_nudge(x=0.1))#brown


### East sample6 - 7 clusters ----
result.7clust.sample6<-gw_clust(ggwr_east_6,maxclusters=7,out3=result.7clust.sample1$out3,out4=result.7clust.sample1$out4)
#plot(result.7clust.sample6$mclust,"BIC")
#result.7clust.sample6$BCRplot
#result.7clust.sample6$Brandtplot
p.7.6<-result.7clust.sample6$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(7,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0.15))+theme(axis.text.x = element_text(size=9,angle = 45, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(7,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.05))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(7,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(7,"Set1")[4],size=1,fatten=2,position=position_nudge(x=-0.15))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(7,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.05))+#orange
  geom_pointrange(aes(Coef,mean6,ymin=mean6-sd6,ymax=mean6+sd6),col=brewer.pal(7,"Set1")[6],size=1,fatten=2,position=position_nudge(x=-0.1))+ #yellow
  geom_pointrange(aes(Coef,mean7,ymin=mean7-sd7,ymax=mean7+sd7),col=brewer.pal(7,"Set1")[7],size=1,fatten=2,position=position_nudge(x=0.1))#brown


### East sample7 - 7 clusters ----
result.7clust.sample7<-gw_clust(ggwr_east_7,maxclusters=7,out3=result.7clust.sample1$out3,out4=result.7clust.sample1$out4)
#plot(result.7clust.sample7$mclust,"BIC")
#result.7clust.sample7$BCRplot
#result.7clust.sample7$Brandtplot
p.7.7<-result.7clust.sample7$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(7,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0.15))+theme(axis.text.x = element_text(size=9,angle = 45, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(7,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.05))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(7,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(7,"Set1")[4],size=1,fatten=2,position=position_nudge(x=-0.15))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(7,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.05))+#orange
  geom_pointrange(aes(Coef,mean6,ymin=mean6-sd6,ymax=mean6+sd6),col=brewer.pal(7,"Set1")[6],size=1,fatten=2,position=position_nudge(x=-0.1))+ #yellow
  geom_pointrange(aes(Coef,mean7,ymin=mean7-sd7,ymax=mean7+sd7),col=brewer.pal(7,"Set1")[7],size=1,fatten=2,position=position_nudge(x=0.1))#brown



### East sample8 - 7 clusters ----
result.7clust.sample8<-gw_clust(ggwr_east_8,maxclusters=7,out3=result.7clust.sample1$out3,out4=result.7clust.sample1$out4)
#plot(result.7clust.sample8$mclust,"BIC")
#result.7clust.sample8$BCRplot
#result.7clust.sample8$Brandtplot
p.7.8<-result.7clust.sample8$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(7,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0.15))+theme(axis.text.x = element_text(size=9,angle = 45, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(7,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.05))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(7,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(7,"Set1")[4],size=1,fatten=2,position=position_nudge(x=-0.15))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(7,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.05))+#orange
  geom_pointrange(aes(Coef,mean6,ymin=mean6-sd6,ymax=mean6+sd6),col=brewer.pal(7,"Set1")[6],size=1,fatten=2,position=position_nudge(x=-0.1))+ #yellow
  geom_pointrange(aes(Coef,mean7,ymin=mean7-sd7,ymax=mean7+sd7),col=brewer.pal(7,"Set1")[7],size=1,fatten=2,position=position_nudge(x=0.1))#brown


### East sample9 - 7 clusters ----
result.7clust.sample9<-gw_clust(ggwr_east_9,maxclusters=7,out3=result.7clust.sample1$out3,out4=result.7clust.sample1$out4)
#plot(result.7clust.sample9$mclust,"BIC")
#result.7clust.sample9$BCRplot
#result.7clust.sample9$Brandtplot
p.7.9<-result.7clust.sample9$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(7,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0.15))+theme(axis.text.x = element_text(size=9,angle = 45, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(7,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.05))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(7,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(7,"Set1")[4],size=1,fatten=2,position=position_nudge(x=-0.15))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(7,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.05))+#orange
  geom_pointrange(aes(Coef,mean6,ymin=mean6-sd6,ymax=mean6+sd6),col=brewer.pal(7,"Set1")[6],size=1,fatten=2,position=position_nudge(x=-0.1))+ #yellow
  geom_pointrange(aes(Coef,mean7,ymin=mean7-sd7,ymax=mean7+sd7),col=brewer.pal(7,"Set1")[7],size=1,fatten=2,position=position_nudge(x=0.1))#brown


### East sample10 - 7 clusters ----
result.7clust.sample10<-gw_clust(ggwr_east_10,maxclusters=7,out3=result.7clust.sample1$out3,out4=result.7clust.sample1$out4)
#plot(result.7clust.sample10$mclust,"BIC")
#result.7clust.sample10$BCRplot
#result.7clust.sample10$Brandtplot
p.7.10<-result.7clust.sample10$coefplot+geom_pointrange(aes(ymin=mean1-sd1,ymax=mean1+sd1),col=brewer.pal(7,"Set1")[1],size=1,fatten=2,position=position_nudge(x=0.15))+theme(axis.text.x = element_text(size=9,angle = 45, hjust = 1),axis.text.y = element_text(size=12))+ylab("")+xlab("")+ # red
  geom_pointrange(aes(Coef,mean2,ymin=mean2-sd2,ymax=mean2+sd2),col=brewer.pal(7,"Set1")[2],size=1,fatten=2,position=position_nudge(x=-0.05))+ #blue
  geom_pointrange(aes(Coef,mean3,ymin=mean3-sd3,ymax=mean3+sd3),col=brewer.pal(7,"Set1")[3],size=1,fatten=2,position=position_nudge(x=0))+ #green
  geom_pointrange(aes(Coef,mean4,ymin=mean4-sd4,ymax=mean4+sd4),col=brewer.pal(7,"Set1")[4],size=1,fatten=2,position=position_nudge(x=-0.15))+ #purple
  geom_pointrange(aes(Coef,mean5,ymin=mean5-sd5,ymax=mean5+sd5),col=brewer.pal(7,"Set1")[5],size=1,fatten=2,position=position_nudge(x=0.05))+#orange
  geom_pointrange(aes(Coef,mean6,ymin=mean6-sd6,ymax=mean6+sd6),col=brewer.pal(7,"Set1")[6],size=1,fatten=2,position=position_nudge(x=-0.1))+ #yellow
  geom_pointrange(aes(Coef,mean7,ymin=mean7-sd7,ymax=mean7+sd7),col=brewer.pal(7,"Set1")[7],size=1,fatten=2,position=position_nudge(x=0.1))#brown



# composite plots ----


# with BCRs
p1<-spplot(pres.ggwr_exponential3, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="",
       col.regions=brewer.pal(4,"Set1")[-3],
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out,col="gainsboro",fill=colsmapBCRs),basemap)
)
p2<-
# with Brandt boreal
spplot(pres.ggwr_exponential3, "classification", do.log = F,
       key.space=list(x=0.5,y=0.9,corner=c(0,1)),
       main="",
       col.regions=brewer.pal(4,"Set1")[-3],
       #sp.layout=list(list("sp.polygons",out,col="lightgrey",fill=alpha(c("#ffffff",rainbow(18)[-1]),0.3)[out$BCR]),basemap)
       sp.layout=list(list("sp.polygons",out2,col="gainsboro",fill=brandtcols),basemap)
)


jpeg("3clustCanada.jpg",width=20,height=5,units="in",quality = 100, res=72)
grid.arrange(p1,p2,nrow=1)
dev.off()


png("3clustCoefs.jpg",width=20,height=25,units="in",res=120)
grid.arrange(p.3.1,p.3.2,p.3.3,p.3.4,p.3.5,p.3.6,p.3.7,p.3.8,p.3.9,p.3.10,nrow=5)
dev.off()

png("3clustBCR.jpg",width=20,height=25,units="in", res=72)
grid.arrange(result.3clust.sample1$BCRplot,
             result.3clust.sample2$BCRplot,
             result.3clust.sample3$BCRplot,
             result.3clust.sample4$BCRplot,
             result.3clust.sample5$BCRplot,
             result.3clust.sample6$BCRplot,
             result.3clust.sample7$BCRplot,
             result.3clust.sample8$BCRplot,
             result.3clust.sample9$BCRplot,
             result.3clust.sample10$BCRplot,
             nrow=5)
dev.off()


png("3clustBrandt.jpg",width=20,height=25,units="in", res=72)
grid.arrange(result.3clust.sample1$Brandtplot,
             result.3clust.sample2$Brandtplot,
             result.3clust.sample3$Brandtplot,
             result.3clust.sample4$Brandtplot,
             result.3clust.sample5$Brandtplot,
             result.3clust.sample6$Brandtplot,
             result.3clust.sample7$Brandtplot,
             result.3clust.sample8$Brandtplot,
             result.3clust.sample9$Brandtplot,
             result.3clust.sample10$Brandtplot,
             nrow=5)
dev.off()

png("3BICplots.jpg",width=12,height=15,units="in", res=150)
par(mfrow=c(5,2))
plot(result.3clust.sample1$mclust,"BIC")
plot(result.3clust.sample2$mclust,"BIC")
plot(result.3clust.sample3$mclust,"BIC")
plot(result.3clust.sample4$mclust,"BIC")
plot(result.3clust.sample5$mclust,"BIC")
plot(result.3clust.sample6$mclust,"BIC")
plot(result.3clust.sample7$mclust,"BIC")
plot(result.3clust.sample8$mclust,"BIC")
plot(result.3clust.sample9$mclust,"BIC")
plot(result.3clust.sample10$mclust,"BIC")
dev.off()

png("5clustBCR.jpg",width=20,height=25,units="in", res=72)
grid.arrange(result.5clust.sample1$BCRplot,
             result.5clust.sample2$BCRplot,
             result.5clust.sample3$BCRplot,
             result.5clust.sample4$BCRplot,
             result.5clust.sample5$BCRplot,
             result.5clust.sample6$BCRplot,
             result.5clust.sample7$BCRplot,
             result.5clust.sample8$BCRplot,
             result.5clust.sample9$BCRplot,
             result.5clust.sample10$BCRplot,
             nrow=5)
dev.off()

png("5clustBrandt.jpg",width=20,height=25,units="in", res=72)
grid.arrange(result.5clust.sample1$Brandtplot,
             result.5clust.sample2$Brandtplot,
             result.5clust.sample3$Brandtplot,
             result.5clust.sample4$Brandtplot,
             result.5clust.sample5$Brandtplot,
             result.5clust.sample6$Brandtplot,
             result.5clust.sample7$Brandtplot,
             result.5clust.sample8$Brandtplot,
             result.5clust.sample9$Brandtplot,
             result.5clust.sample10$Brandtplot,
             nrow=5)
dev.off()

png("5clustCoefs.jpg",width=20,height=25,units="in", res=120)
grid.arrange(p.5.1,p.5.2,p.5.3,p.5.4,p.5.5,p.5.6,p.5.7,p.5.8,p.5.9,p.5.10,nrow=5)
dev.off()


png("5BICplots.jpg",width=12,height=15,units="in", res=150)
par(mfrow=c(5,2))
plot(result.5clust.sample1$mclust,"BIC")
plot(result.5clust.sample2$mclust,"BIC")
plot(result.5clust.sample3$mclust,"BIC")
plot(result.5clust.sample4$mclust,"BIC")
plot(result.5clust.sample5$mclust,"BIC")
plot(result.5clust.sample6$mclust,"BIC")
plot(result.5clust.sample7$mclust,"BIC")
plot(result.5clust.sample8$mclust,"BIC")
plot(result.5clust.sample9$mclust,"BIC")
plot(result.5clust.sample10$mclust,"BIC")
dev.off()

png("7clustBCR.jpg",width=20,height=25,units="in", res=72)
grid.arrange(result.7clust.sample1$BCRplot,
             result.7clust.sample2$BCRplot,
             result.7clust.sample3$BCRplot,
             result.7clust.sample4$BCRplot,
             result.7clust.sample5$BCRplot,
             result.7clust.sample6$BCRplot,
             result.7clust.sample7$BCRplot,
             result.7clust.sample8$BCRplot,
             result.7clust.sample9$BCRplot,
             result.7clust.sample10$BCRplot,
             nrow=5)
dev.off()

png("7clustBrandt.jpg",width=20,height=25,units="in", res=72)
grid.arrange(result.7clust.sample1$Brandtplot,
             result.7clust.sample2$Brandtplot,
             result.7clust.sample3$Brandtplot,
             result.7clust.sample4$Brandtplot,
             result.7clust.sample5$Brandtplot,
             result.7clust.sample6$Brandtplot,
             result.7clust.sample7$Brandtplot,
             result.7clust.sample8$Brandtplot,
             result.7clust.sample9$Brandtplot,
             result.7clust.sample10$Brandtplot,
             nrow=5)
dev.off()

png("7BICplots.jpg",width=12,height=15,units="in", res=150)
par(mfrow=c(5,2))
plot(result.7clust.sample1$mclust,"BIC")
plot(result.7clust.sample2$mclust,"BIC")
plot(result.7clust.sample3$mclust,"BIC")
plot(result.7clust.sample4$mclust,"BIC")
plot(result.7clust.sample5$mclust,"BIC")
plot(result.7clust.sample6$mclust,"BIC")
plot(result.7clust.sample7$mclust,"BIC")
plot(result.7clust.sample8$mclust,"BIC")
plot(result.7clust.sample9$mclust,"BIC")
plot(result.7clust.sample10$mclust,"BIC")
dev.off()

png("7clustCoefs.jpg",width=20,height=25,units="in", res=120)
grid.arrange(p.7.1,p.7.2,p.7.3,p.7.4,p.7.5,p.7.6,p.7.7,p.7.8,p.7.9,p.7.10,nrow=5)
dev.off()


png("sample1BCR.jpg", width=20,height=15, units="in",  res=120)
grid.arrange(result.3clust.sample1$BCRplot,p.3.1,result.5clust.sample1$BCRplot,p.5.1,result.7clust.sample1$BCRplot,p.7.1,nrow=3)
dev.off()

png("sample2BCR.jpg", width=20,height=15, units="in",  res=120)
grid.arrange(result.3clust.sample2$BCRplot,p.3.1,result.5clust.sample2$BCRplot,p.5.1,result.7clust.sample2$BCRplot,p.7.1,nrow=3)
dev.off()

png("sample3BCR.jpg", width=20,height=15, units="in",  res=120)
grid.arrange(result.3clust.sample3$BCRplot,p.3.1,result.5clust.sample3$BCRplot,p.5.1,result.7clust.sample3$BCRplot,p.7.1,nrow=3)
dev.off()

png("sample4BCR.jpg", width=20,height=15, units="in",  res=120)
grid.arrange(result.3clust.sample4$BCRplot,p.3.1,result.5clust.sample4$BCRplot,p.5.1,result.7clust.sample4$BCRplot,p.7.1,nrow=3)
dev.off()

png("sample5BCR.jpg", width=20,height=15, units="in",  res=120)
grid.arrange(result.3clust.sample5$BCRplot,p.3.1,result.5clust.sample5$BCRplot,p.5.1,result.7clust.sample5$BCRplot,p.7.1,nrow=3)
dev.off()

png("sample6BCR.jpg", width=20,height=15, units="in",  res=120)
grid.arrange(result.3clust.sample6$BCRplot,p.3.1,result.5clust.sample6$BCRplot,p.5.1,result.7clust.sample6$BCRplot,p.7.1,nrow=3)
dev.off()

png("sample7BCR.jpg", width=20,height=15, units="in",  res=120)
grid.arrange(result.3clust.sample7$BCRplot,p.3.1,result.5clust.sample7$BCRplot,p.5.1,result.7clust.sample7$BCRplot,p.7.1,nrow=3)
dev.off()

png("sample8BCR.jpg", width=20,height=15, units="in",  res=120)
grid.arrange(result.3clust.sample8$BCRplot,p.3.1,result.5clust.sample8$BCRplot,p.5.1,result.7clust.sample8$BCRplot,p.7.1,nrow=3)
dev.off()

png("sample9BCR.jpg", width=20,height=15, units="in",  res=120)
grid.arrange(result.3clust.sample9$BCRplot,p.3.1,result.5clust.sample9$BCRplot,p.5.1,result.7clust.sample9$BCRplot,p.7.1,nrow=3)
dev.off()

png("sample10BCR.jpg", width=20,height=15, units="in",  res=120)
grid.arrange(result.3clust.sample10$BCRplot,p.3.1,result.5clust.sample10$BCRplot,p.5.1,result.7clust.sample10$BCRplot,p.7.1,nrow=3)
dev.off()


## Export shapefiles ####

writeOGR(obj=result.5clust.sample2$pres.ggwr, dsn="output_shapefiles", layer="sample2_5clust", driver="ESRI Shapefile")

pres.ggwr_exponential3west<-pres.ggwr_exponential3[c(which(pres.ggwr_exponential3$classification==2),which(pres.ggwr_exponential3$classification==3)),]



writeOGR(obj=pres.ggwr_exponential3west, dsn="output_shapefiles", layer="west", driver="ESRI Shapefile")
