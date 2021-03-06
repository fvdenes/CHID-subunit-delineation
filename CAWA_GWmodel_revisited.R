# Revisiting sub-unit delineation analysis - 08/07/19 

# Load packages and data
library(raster)
library(sf)
library(maptools)
library(dplyr)
library(data.table)
library(reshape2)
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

# Prepare dataset ####
## updated avian dataset including recent atlas data

load("D:/CHID subunit delineation/BAM_data_package_August2019.RData")

head(PKEYcombo)
head(PCcombo)
head(SScombo)
head(offcombo)

#
SScombo2<-SpatialPointsDataFrame(coords=SScombo[,2:3],data=SScombo,proj4string = LCC)


##MODIS-based landcover (250-m)
nalc2005<-raster("M:/DataStuff/SpatialData/LCC05_NALCMS2010/Land_Cover_MXD/NA_LandCover_2005/data/NA_LandCover_2005/NA_LandCover_2005_LCC.img")
SScombo2<-spTransform(SScombo2,proj4string(nalc2005))
SScombo <- cbind(SScombo,"nalc"=extract(nalc2005,SScombo2@coords))

## Tree cover (density, Global 30m Landsat Tree Canopy Version 4)
treecover<-raster("D:/TreeCover/treecoverlcc1.tif")
SScombo2<-spTransform(SScombo2,proj4string(treecover))
SScombo <- cbind(SScombo,"treecover"=extract(treecover,SScombo2@coords))
SScombo$treecover[which(SScombo$treecover>250)]<-0 # converting values >250 (which represent no forest data) into 0s

## Forest height
height<-raster("D:/ForestHeightLidar/Simard_Pinto_3DGlobalVeg_L3C.tif")
SScombo2<-spTransform(SScombo2,proj4string(height))
SScombo <- cbind(SScombo,"height"=extract(height,SScombo2@coords))

## Road on/off
road<- raster("D:/roadonoff/roadonoff1.tif")
mr <- c(1, 2500000, 1,  NA, NA, 0)
rcroad <- matrix(mr, ncol=3, byrow=TRUE)
rrc <- reclassify(road,rcroad)
SScombo2<-spTransform(SScombo2,proj4string(road))
SScombo <- cbind(SScombo,"road"=extract(rrc,SScombo2@coords))  
  
## Climate data- upload and resample to 250m resolution to match other layers 
climate2010 <- list.files("D:/ClimateAdaptWest/baseline19812010/",pattern="asc$")
setwd("D:/ClimateAdaptWest/baseline19812010/")
clim2010 <- stack(raster(climate2010[1]))
for (i in 2:length(climate2010)) { clim2010 <- addLayer(clim2010, raster(climate2010[i]))}

LCC2<-CRS( "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
proj4string(clim2010)<-LCC2

SScombo2<-spTransform(SScombo2,LCC2)
SScombo <- cbind(SScombo,extract(clim2010,SScombo2@coords))

climate.tifs<-list.files("D:/ClimateAdaptWest/baseline19812010/",pattern="tif")
climtif <- stack(raster(climate.tifs[1]))
climtif <- addLayer(climtif, raster(climate.tifs[2]))

SScombo2<-spTransform(SScombo2,proj4string(climtif))
SScombo <- cbind(SScombo,extract(climtif,SScombo2@coords))


## Landform and topoedaphic 
lf <- raster("D:/topo/lf_lcc1.tif")
SScombo2<-spTransform(SScombo2,LCC)
SScombo <- cbind(SScombo,"landform"=extract(lf,SScombo2@coords))

SScombo$landform[which(SScombo$landform>0&SScombo$landform<1)]<-0
lf_classes<-data.frame(value=0:9,label=factor(c("water","valley","hilltop.in.valley","headwaters", "ridges.and.peaks","plain","local.ridge.in.plain","local.valley.in.plain","gentle.slopes","steep.slopes")))
lfclass<-lf_classes$label[match(SScombo$landform,lf_classes$value)]
SScombo$landform<-lfclass





#save.image("D:/CHID subunit delineation/subunits_revisited.RData")
##load("D:/CHID subunit delineation/subunits_revisited.RData")

#Terrain Ruggedness Index (TRI) is the mean of the absolute differences in elevation between a focal cell and its 8 surrounding cells. It quantifies the total elevation change across the 3×3 cells. Flat areas have a value of zero whereas mountain areas with steep ridges have positive values, which can be greater than 2000m in the Himalaya area. 
#Topographic Position Index (TPI) is the difference between the elevation of a focal cell and the mean of its 8 surrounding cells. Positive and negative values correspond to ridges and valleys, respectively, while zero values correspond generally to flat areas (with the exception of a special case where a focal cell with a value 5 can have surrounding cells with values of 4×1 and 4×9, resulting in a TPI of 0). 
#Roughness is expressed as the largest inter-cell difference of a focal cell and its 8 surrounding cells.
#Slope and aspect are decomposed into 3-dimensional vector components (in the x, y, and z directions) using standard trigonometric operators, and by calculating the resultant vector magnitude within a user-specified moving window size 

TPI <- raster("D:/topo/tpiLCC.tif")
SScombo <- cbind(SScombo,"TPI"=extract(TPI,SScombo2@coords))
TRI <- raster("D:/topo/triLCC.tif")
SScombo <- cbind(SScombo,"TRI"=extract(TRI,SScombo2@coords))
slope <- raster("D:/topo/slopeLCC.tif")
SScombo <- cbind(SScombo,"slope"=extract(slope,SScombo2@coords))
roughness <- raster("D:/topo/roughLCC.tif")
SScombo <- cbind(SScombo,"roughness"=extract(roughness,SScombo2@coords))


## Preparing data table
PC<-inner_join(PCcombo,PKEYcombo,by=c("PKEY"))
PC<-inner_join(PC,SScombo,by="SS")

survey <- aggregate(PC$ABUND, by=list("PKEY"=PC$PKEY,"SS"=PC$SS,"SPECIES"=PC$SPECIES, "ARU"=PC$ARU, " YEAR"=PC$YEAR), FUN=sum) 

speclist<-"CAWA"

for (j in 1:length(speclist)) {
  specoff <- filter(offcombo, SPECIES==as.character(speclist[j]))
  specoff <- distinct(specoff,PKEY, .keep_all=TRUE) 
  
  specdat <- filter(PC, SPECIES == as.character(speclist[j]))
  specdatx <- aggregate(specdat$ABUND,by=list("PKEY"=specdat$PKEY,"SS"=specdat$SS), FUN=sum)
  names(specdatx)[3] <- "ABUND"
  dat1 <- right_join(specdatx,survey[,1:4],by=c("SS","PKEY")) 
  dat1$SPECIES <- as.character(speclist[j])
  dat1$ABUND <- as.integer(ifelse(is.na(dat1$ABUND),0,dat1$ABUND))
  dat11 <- distinct(dat1,SS,.keep_all=TRUE) 
  s1 <- left_join(dat11,specoff, by=c("SPECIES","PKEY"))
  SScombo2 <- distinct(SScombo, SS,.keep_all = TRUE)
  d1 <- left_join(s1, SScombo2, by=c("SS"))
  }

#load("D:/CHID subunit delineation/subunits_revisited.RData")

## Optimize landcover class variable with union of samples 
ltnalc <- read.csv("C:/Users/voeroesd/Documents/Repos/bamanalytics/lookup/nalcms.csv")
d1$LCclass<-ltnalc$Label[match(d1$nalc, ltnalc$Value)]
d1$LCclass<-factor(as.character(d1$LCclass))
d1<-d1[-which(is.na(d1$LCclass)),]


d1<-d1[complete.cases(d1),]

## Map and sample points
# basemap<-readOGR("D:/CHID subunit delineation/province_state_lcc.shp") ## Basemap for point plots
# BCRs<- readOGR("D:/CHID subunit delineation/bcrfinallcc.shp")
# brandt<-readOGR("D:/CHID subunit delineation/BRANDT_diss_type.shp")
# brandt<-spTransform(brandt,CRS("+proj=lcc +lat_1=49 +lat_2=77 +lat_0=0 +lon_0=-95 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"))
# 
# 
# allpoints <- SpatialPointsDataFrame(coords=cbind(d1$X,d1$Y),proj4string =  LCC,data=d1)
# 
# spplot(allpoints,"ABUND",do.log = F,legendEntries=as.character(0:4),
#        key.space=list(x=0.02,y=0.3,corner=c(0,1)),
#        sp.layout=list(basemap))
# 
# 
# spplot(allpoints[which(allpoints$ABUND>0),],"ABUND",do.log = F,cuts=4,legendEntries=as.character(1:4),
#        key.space=list(x=0.5,y=0.9,corner=c(0,1)),
#        sp.layout=list(basemap))


## Obtain 10 quadrant-stratified random samples from full dataset: 

## A function to randomly sample from quadrants defined by latitude and longitude intervals.
sample_quad<-function(data,N=10,Xintervals=100,Yintervals=60,seed=123,onlypres=F){
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
    data2<-data[which(data$ABUND>0),]
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
## A function to plot sampled points
plot_sampled<-function(obj,index,onlypres=T,show="ABUND"){
  thin_cawa_mm<-obj[index,]
  mmsp <- SpatialPointsDataFrame(coords=cbind(thin_cawa_mm$X,thin_cawa_mm$Y),proj4string =  LCC,data=thin_cawa_mm)
  
  if(onlypres==T){
    spplot(mmsp[which(mmsp$ABUND>0),],show,do.log = F,cuts=4,legendEntries=as.character(1:4),
           key.space=list(x=0.5,y=0.9,corner=c(0,1)),
           sp.layout=list(basemap))
  }
  else{
    spplot(mmsp,show,do.log = F,legendEntries=as.character(0:4),
           key.space=list(x=0.02,y=0.3,corner=c(0,1)),
           sp.layout=list(basemap))
  }
}



## Apply sampling function 10 times
samples<-as.list(rep(NA,10))
names(samples)<-paste0("sample",1:10)

samples$sample1<- sample_quad(d1,seed=1)
samples$sample2<- sample_quad(d1,seed=2)
samples$sample3<- sample_quad(d1,seed=3)
samples$sample4<- sample_quad(d1,seed=4)
samples$sample5<- sample_quad(d1,seed=5)
samples$sample6<- sample_quad(d1,seed=6)
samples$sample7<- sample_quad(d1,seed=7)
samples$sample8<- sample_quad(d1,seed=8)
samples$sample9<- sample_quad(d1,seed=9)
samples$sample10<- sample_quad(d1,seed=10)

d2 <-lapply(samples,function(x){d1[x,]})


union_samples<-union(union(union(union(union(union(union(union(union(samples$sample1,samples$sample2),samples$sample3),samples$sample4),samples$sample5),samples$sample6),samples$sample7),samples$sample8),samples$sample9),samples$sample10)
d1_union<-d1[union_samples,]

ol<- optilevels(y=d1_union$ABUND, x=d1_union$LCclass, dist="poisson", offset=d1_union$logoffset)

d1_union$LCclass2<-d1_union$LCclass
levels(d1_union$LCclass2)<-ol$level[[length(ol$level)]]

lc_optim<-function(d2){
  d2$LCclass2<-d2$LCclass
  levels(d2$LCclass2)<-ol$level[[length(ol$level)]]
  return(d2)
}
d2<-lapply(d2,lc_optim)

# # centering  variables for each sample based on their respective means from the union of all samples
centercovs<-function(sample){
  sample$TPI<-scale(sample$TPI)
  sample$TRI<-scale(sample$TRI)
  sample$slope<-scale(sample$slope)
  sample$roughness<-scale(sample$roughness)
  return(sample)
}
d2<-lapply(d2,centercovs)

# Create spatial dataframe from data samples

sp_fun<-function(d2){
  sp<- SpatialPointsDataFrame(coords=cbind(d2$X,d2$Y),proj4string = LCC,data=d2)
}

d2sp <- lapply(d2,sp_fun)

save.image("D:/CHID subunit delineation/subunits_revisited.RData")


rm(list=setdiff(ls(),c("d1","d1_union","d2","LCC","sample_all","plot_sampled","sample_quad","basemap","d2","d2sp"))) # UPDATE OBJECT LIST HERE!
gc()


#save.image("D:/CHID subunit delineation/subunits_revisited.RData")
#load("D:/CHID subunit delineation/subunits_revisited.RData")

# Geographically weighted models ####


## Specify distance matrices and optimize bandwidths ####
DM <- gw.dist(dp.locat = coordinates(d2sp[[1]]), longlat=TRUE)

# define formulas for bandwidth optimization
formula_habclim <- formula(ABUND ~  (LCclass2-1) + treecover + height + TPI + road + EMT + DD_0 + DD5 + TD + MSP + CMI + CMIJJA + ARU + offset(d2sp[[1]]$logoffset))
formula_hab<- formula(ABUND ~  (LCclass2-1) + treecover + height + TPI + road + ARU + offset(d2sp[[1]]$logoffset))

# bandwidths
#bw.ggwr.exponential.habclim.adap <- bw.ggwr(formula_habclim, data = d2sp[[1]], approach = "AICc", kernel = "exponential", adaptive = TRUE, family="poisson", longlat = TRUE, dMat=DM) 
#save.image("D:/CHID subunit delineation/subunits_revisited.RData")
bw.ggwr.gaussian.habclim.adap <- bw.ggwr(formula_habclim, data = d2sp[[1]], approach = "AICc", kernel = "gaussian", adaptive = TRUE, family="poisson", longlat = TRUE, dMat=DM) 
save.image("D:/CHID subunit delineation/subunits_revisited.RData")

#bw.ggwr.exponential.hab <- bw.ggwr(formula_hab, data = d2sp[[1]], approach = "AICc", kernel = "exponential", adaptive = FALSE, family="poisson", longlat = TRUE, dMat=DM) 
#save.image("D:/CHID subunit delineation/subunits_revisited.RData")
bw.ggwr.gaussian.hab <- bw.ggwr(formula_hab, data = d2sp[[1]], approach = "AICc", kernel = "gaussian", adaptive = FALSE, family="poisson", longlat = TRUE, dMat=DM) 
save.image("D:/CHID subunit delineation/subunits_revisited.RData")

#bw.ggwr.exponential.habclim <- bw.ggwr(formula_habclim, data = d2sp[[1]], approach = "AICc", kernel = "exponential", adaptive = FALSE, family="poisson", longlat = TRUE, dMat=DM) 
#save.image("D:/CHID subunit delineation/subunits_revisited.RData")
bw.ggwr.gaussian.habclim <- bw.ggwr(formula_habclim, data = d2sp[[1]], approach = "AICc", kernel = "gaussian", adaptive = FALSE, family="poisson", longlat = TRUE, dMat=DM) 
save.image("D:/CHID subunit delineation/subunits_revisited.RData")

#bw.ggwr.exponential.hab.adap <- bw.ggwr(formula_hab, data = d2sp[[1]], approach = "AICc", kernel = "exponential", adaptive = TRUE, family="poisson", longlat = TRUE, dMat=DM) 
#save.image("D:/CHID subunit delineation/subunits_revisited.RData")
bw.ggwr.gaussian.hab.adap <- bw.ggwr(formula_hab, data = d2sp[[1]], approach = "AICc", kernel = "gaussian", adaptive = TRUE, family="poisson", longlat = TRUE, dMat=DM) 
save.image("D:/CHID subunit delineation/subunits_revisited.RData")

rm(d2,DM,basemap)
gc()

# define functions to fit models
ggwr_sample<-function(formula,index,kernel="gaussian", bw=NULL,adaptive=FALSE,out.model=FALSE,save.path="D:/CHID subunit delineation/output_revisited/"){
  
  d3<-d2sp[[index]]
  
  DM <- gw.dist(dp.locat = coordinates(d3), longlat=TRUE)
  
  samplename<-paste0("ggwr_",kernel,index)
  
  # Fit models
  if(kernel=="gaussian"){
    
    if(is.null(bw)==TRUE){
      bw.gaussian <- bw.ggwr(formula, data = d3, approach = "AICc", kernel = "exponential", adaptive = adaptive, family="poisson", longlat = TRUE, dMat=DM) 
      save(bw.gaussian, file=paste0(save.path,paste0(samplename,"bandwidth.csv")))
    }
    
    else{
      bw.gaussian<-bw
    }
    
    ggwr_gaussian<-ggwr.basic(formula, data = d3, bw = bw.gaussian, kernel = "gaussian", adaptive = adaptive, longlat=TRUE, family="poisson", dMat=DM)
    saveRDS(ggwr_gaussian, file=paste0(save.path,paste0(samplename,".R")))
    if(out.model==TRUE){
      return(ggwr_gaussian)
    }
  }
  
  if(kernel=="exponential"){
    
    if(is.null(bw)==TRUE){
      bw.exponential <- bw.ggwr(formula, data = d3, approach = "AICc", kernel = "exponential", adaptive = adaptive, family="poisson", longlat = TRUE, dMat=DM) 
      save(bw.exponential_grid_east, file=paste0(save.path,paste0(samplename,"bandwidth.csv")))
    }
    
    else{
      bw.exponential<-bw
    }
    
    ggwr_exponential<-ggwr.basic(formula, data = d3, bw = bw.exponential, kernel = "exponential", adaptive = adaptive, longlat=TRUE, family="poisson", dMat=DM)
    saveRDS(ggwr_exponential, file=paste0(save.path,paste0(samplename,".R")))
    if(out.model==TRUE){
      return(ggwr_exponential)
    }
  }
} 

ggwr_sample_grid<-function(formula,index,kernel="gaussian", bw=NULL,adaptive=FALSE,out.model=FALSE,save.path="D:/CHID subunit delineation/output_revisited/"){
  
  d3<-d2sp[[index]]
  
  grd <-SpatialGrid(points2grid(d3,tolerance=0.49),proj4string = LCC)
  
  DM <- gw.dist(dp.locat = coordinates(d3),rp.locat=coordinates(grd), longlat=TRUE)
  
  samplename<-paste0("ggwr_grid_",kernel,index)
  
  # Fit models
  if(kernel=="gaussian"){
    
    if(is.null(bw)==TRUE){
      bw.gaussian <- bw.ggwr(formula, data = d3, approach = "AICc", kernel = "exponential", adaptive = adaptive, family="poisson", longlat = TRUE, dMat=DM) 
      save(bw.gaussian, file=paste0(save.path,paste0(samplename,"bandwidth.csv")))
    }
    
    else{
      bw.gaussian<-bw
    }
    
    ggwr_gaussian<-ggwr.basic(formula, data = d3, bw = bw.gaussian, kernel = "gaussian", adaptive = adaptive, longlat=TRUE, family="poisson", dMat=DM)
    saveRDS(ggwr_gaussian, file=paste0(save.path,paste0(samplename,"_grid",".R")))
    if(out.model==TRUE){
      return(ggwr_gaussian)
    }
  }
  
  if(kernel=="exponential"){
    
    if(is.null(bw)==TRUE){
      bw.exponential <- bw.ggwr(formula, data = d3, approach = "AICc", kernel = "exponential", adaptive = adaptive, family="poisson", longlat = TRUE, dMat=DM) 
      save(bw.exponential_grid_east, file=paste0(save.path,paste0(samplename,"bandwidth.csv")))
    }
    
    else{
      bw.exponential<-bw
    }
    
    ggwr_exponential<-ggwr.basic(formula, data = d3, bw = bw.exponential, kernel = "exponential", adaptive = adaptive, longlat=TRUE, family="poisson", dMat=DM)
    saveRDS(ggwr_exponential, file=paste0(save.path,paste0(samplename,"_grid",".R")))
    if(out.model==TRUE){
      return(ggwr_exponential)
    }
  }
} 

# redefine formulas (update the offset term to work with the above-defined functions)
formula_habclim <- formula(ABUND ~  (LCclass2-1) + treecover + height + TPI + road + EMT + DD_0 + DD5 + TD + MSP + CMI + CMIJJA + ARU + offset(
  logoffset))
formula_hab<- formula(ABUND ~  (LCclass2-1) + treecover + height + TPI + road + ARU + offset(logoffset))

# fit models ####
ggwr_sample(formula= formula_habclim ,index=1,kernel="gaussian", adaptive=T, bw = bw.ggwr.gaussian.habclim.adap)
save.image("D:/CHID subunit delineation/subunits_revisited.RData")
ggwr_sample(formula= formula_hab ,index=3,kernel="gaussian", adaptive=T, bw = bw.ggwr.gaussian.hab.adap)
gc()
ggwr_sample(formula= formula_hab ,index=5,kernel="gaussian", adaptive=T, bw = bw.ggwr.gaussian.hab.adap)
ggwr_sample(formula= formula_hab ,index=7,kernel="gaussian", adaptive=T, bw = bw.ggwr.gaussian.hab.adap)
gc()
ggwr_sample(formula= formula_hab ,index=9,kernel="gaussian", adaptive=T, bw = bw.ggwr.gaussian.hab.adap)
ggwr_sample(formula= formula_hab ,index=2,kernel="gaussian", adaptive=T, bw = bw.ggwr.gaussian.hab.adap)
gc()
ggwr_sample(formula= formula_hab ,index=4,kernel="gaussian", adaptive=T, bw = bw.ggwr.gaussian.hab.adap)
ggwr_sample(formula= formula_hab ,index=6,kernel="gaussian", adaptive=T, bw = bw.ggwr.gaussian.hab.adap)
gc()
ggwr_sample(formula= formula_hab ,index=8,kernel="gaussian", adaptive=T, bw = bw.ggwr.gaussian.hab.adap)
ggwr_sample(formula= formula_hab ,index=10,kernel="gaussian", adaptive=T, bw = bw.ggwr.gaussian.hab.adap)
gc()

ggwr_sample(index=1,samplename="ggwr_expon_1",kernel="exponential", cell.means = T, bw = bw.ggwr.exponential)
ggwr_sample(index=3,samplename="ggwr_expon_3",kernel="exponential", cell.means = T, bw = bw.ggwr.exponential)
#save.image("D:/CHID subunit delineation/subunits_revisited.RData")
gc()
ggwr_sample(index=5,samplename="ggwr_expon_5",kernel="exponential", cell.means = T, bw = bw.ggwr.exponential)
ggwr_sample(index=7,samplename="ggwr_expon_7",kernel="exponential", cell.means = T, bw = bw.ggwr.exponential)
#save.image("D:/CHID subunit delineation/subunits_revisited.RData")
gc()
ggwr_sample(index=9,samplename="ggwr_expon_9",kernel="exponential", cell.means = T, bw = bw.ggwr.exponential)
ggwr_sample(index=2,samplename="ggwr_expon_2",kernel="exponential", cell.means = T, bw = bw.ggwr.exponential)
#save.image("D:/CHID subunit delineation/subunits_revisited.RData")
gc()
ggwr_sample(index=4,samplename="ggwr_expon_4",kernel="exponential", cell.means = T, bw = bw.ggwr.exponential)
ggwr_sample(index=6,samplename="ggwr_expon_6",kernel="exponential", cell.means = T, bw = bw.ggwr.exponential)
#save.image("D:/CHID subunit delineation/subunits_revisited.RData")
gc()
ggwr_sample(index=8,samplename="ggwr_expon_8",kernel="exponential", cell.means = T, bw = bw.ggwr.exponential)
ggwr_sample(index=10,samplename="ggwr_expon_10",kernel="exponential", cell.means = T, bw = bw.ggwr.exponential)
#save.image("D:/CHID subunit delineation/subunits_revisited.RData")
gc()

## load back fitted model objects
ggwr_expon1<-readRDS("D:/CHID subunit delineation/output_revisited/ggwr_expon_1.R")
ggwr_expon2<-readRDS("D:/CHID subunit delineation/output_revisited/ggwr_expon_2.R")
ggwr_expon3<-readRDS("D:/CHID subunit delineation/output_revisited/ggwr_expon_3.R")
ggwr_expon4<-readRDS("D:/CHID subunit delineation/output_revisited/ggwr_expon_4.R")
ggwr_expon5<-readRDS("D:/CHID subunit delineation/output_revisited/ggwr_expon_5.R")
ggwr_expon6<-readRDS("D:/CHID subunit delineation/output_revisited/ggwr_expon_6.R")
ggwr_expon7<-readRDS("D:/CHID subunit delineation/output_revisited/ggwr_expon_7.R")
ggwr_expon8<-readRDS("D:/CHID subunit delineation/output_revisited/ggwr_expon_8.R")
ggwr_expon9<-readRDS("D:/CHID subunit delineation/output_revisited/ggwr_expon_9.R")
ggwr_expon10<-readRDS("D:/CHID subunit delineation/output_revisited/ggwr_expon_10.R")

ggwr_gaussian1<-readRDS("D:/CHID subunit delineation/output_revisited/ggwr_gaussian1.R")
ggwr_gaussian2<-readRDS("D:/CHID subunit delineation/output_revisited/ggwr_gaussian2.R")
ggwr_gaussian3<-readRDS("D:/CHID subunit delineation/output_revisited/ggwr_gaussian3.R")
ggwr_gaussian4<-readRDS("D:/CHID subunit delineation/output_revisited/ggwr_gaussian4.R")
ggwr_gaussian5<-readRDS("D:/CHID subunit delineation/output_revisited/ggwr_gaussian5.R")
ggwr_gaussian6<-readRDS("D:/CHID subunit delineation/output_revisited/ggwr_gaussian6.R")
ggwr_gaussian7<-readRDS("D:/CHID subunit delineation/output_revisited/ggwr_gaussian7.R")
ggwr_gaussian8<-readRDS("D:/CHID subunit delineation/output_revisited/ggwr_gaussian8.R")
ggwr_gaussian9<-readRDS("D:/CHID subunit delineation/output_revisited/ggwr_gaussian9.R")
ggwr_gaussian10<-readRDS("D:/CHID subunit delineation/output_revisited/ggwr_gaussian10.R")


ggwr_gaussian1
image(ggwr_gauss1$SDF,"CMI")
plot(d2sp$sample1,add=T)
plot(basemap,add=T)


## AVeraging GWmodel predictions among samples
### create "blank" model raster
b2011 <- list.files("D:/Beaudoin/2011",pattern="tif$")
setwd("D:/Beaudoin/2011/")
bed1<-raster(b2011[1])
plot(bed1)
bed1
bed2<-raster(ext=extent(bed1),crs=proj4string(bed1),resolution=1000)
rm(bed1)
gc()

### fill raster stacks with estimated coefficient values from GWmodel from each data sample (only from points where species was detected)

files<-list.files("D:/CHID subunit delineation/output_revisited",pattern="gauss") # select GW model save.files

files

for(k in 1:length(files)){ 
  gwmod<- readRDS(paste0("D:/CHID subunit delineation/output_revisited/",files[k]))
  grd<-stack(rasterize(x=gwmod$SDF[which(gwmod$SDF$y>0),],y=bed2,field=gwmod$SDF@data[which(gwmod$SDF$y>0),1],fun=mean))
  for (i in 2:23) {grd <- addLayer(grd, rasterize(x=gwmod$SDF,y=bed2,field=gwmod$SDF@data[,i],fun=mean))} # habitat covariates are columns in the SpatialPointsDataFrame
  names(grd)<-names(gwmod$SDF)[1:23]
  
  assign(x=paste0("grid",k),grd)
  rm(grd,gwmod)
  gc()
}

### Average coefficient estimates within each raster cell

set<-ls()[grep("grid",ls())]
z<-mean(stack(lapply(mget(set),function(x){x[[1]]})),na.rm=TRUE)
for(i in 2:nlayers(grid1)){
  z<-addLayer(z,mean(stack(lapply(mget(set),function(x){x[[i]]})),na.rm=TRUE))
}
names(z)<-names(grid1)
writeRaster(z, filename="D:/CHID subunit delineation/output_revisited/means.grd", format="raster",overwrite=TRUE)

stack1<-stack("D:/CHID subunit delineation/output_revisited/means.grd")
stack1
plot(stack1)


# Clustering ####

## Extract data from grid cells into table
df1<-as.data.frame(stack1,na.rm=TRUE,xy=T)
names(df1)
head(df1)

mclustBIC<- mclustBIC(data=df1[,1:12],G=1:20)
mclustBIC
summary(mclustBIC)
plot(mclustBIC)

# explanation of why highest BIC is best model in Mclust: https://stats.stackexchange.com/questions/237220/mclust-model-selection 

mclustMOD<- Mclust(data=df1[,1:12],G=1:20, x=mclustBIC) 
summary(mclustMOD)

## Create raster with classifications from clustering model:
spdf1<-SpatialPointsDataFrame(coords=df1[,24:25],proj4string = CRS(proj4string(stack1)),data=data.frame(Classification=as.factor(mclustMOD$classification)))
spplot(spdf1)

rast<-stack1$LCclass2Agr.Shrub

class<-as.factor(mclustMOD$classification)

rast2 <- calc(rast, function(x) ifelse(is.na(x)==F, class, NA) ) 


