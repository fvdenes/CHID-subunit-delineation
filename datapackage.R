## Define root folder where data are stored
ROOT <- "e:/peter/bam/Apr2016"
## Load required packages
library(mefa4)
## Load functions kept in separate file
source("~/repos/bamanalytics/R/dataprocessing_functions.R")

load(file.path(ROOT, "out", paste0("data_package_2016-12-01.Rdata")))
load(file.path(ROOT, "out", "offsets-v3_2017-04-19.Rdata"))

TAX <- nonDuplicated(TAX, Species_ID, TRUE)
TAX <- droplevels(TAX[SPP,])
OFF <- OFF[,SPP]

YY <- Xtab(ABUND ~ PKEY + SPECIES, PCTBL)
YY <- YY[,SPP]
YYSS <- Xtab(ABUND ~ SS + SPECIES, PCTBL)
YYSS[YYSS > 1] <- 1
YYSS <- YYSS[,SPP]
rm(PCTBL)

PKEY <- PKEY[,c("PKEY", "SS", "PCODE", "METHOD", "SITE", "ROUND", "YEAR",
    "ROAD")]

with(SS, table(ROADCLASS, PAVSTATUS))
with(SS, table(NBRLANES, PAVSTATUS))
with(SS, table(ROADCLASS, NBRLANES))

SS$ROAD_OK <- SS$NBRLANES <= 2
SS$HAB_OK <- !is.na(SS$HAB_NALC2)

SS$d2road <- NULL
SS$ROADCLASS <- NULL
SS$NBRLANES <- NULL
SS$PAVSTATUS <- NULL

hist(SS$BEADTotalL)
hist(SS$BEADtotPol)
## linear features
SS$LIN <- log(SS$BEADTotalL + 1)
SS$BEADTotalL <- NULL
## Polygon is fine as it is (0-1)
SS$POL <- SS$BEADtotPol
SS$BEADtotPol <- NULL

## slope, cti, elev
SS$SLP <- sqrt(SS$slp90)
SS$SLP2 <- SS$SLP^2
SS$CTI <- (SS$cti90 - 8) / 4
SS$CTI2 <- SS$CTI^2
SS$ELV <- SS$elv90 / 1000
SS$ELV2 <- SS$ELV^2
SS$slp90 <- NULL
SS$cti90 <- NULL
SS$elv90 <- NULL

## height
SS$HGTorig <- SS$HEIGHTSIMARD / 25
SS$HEIGHTSIMARD <- NULL

SS$TR3 <- SS$TREE3
SS$TREE3 <- NULL

## climate
SS$CMIJJA <- (SS$CMIJJA - 0) / 50
SS$CMI <- (SS$CMI - 0) / 50
SS$TD <- (SS$TD - 300) / 100
SS$DD0 <- (SS$DD01 - 1000) / 1000
SS$DD5 <- (SS$DD51 - 1600) / 1000
SS$EMT <- (SS$EMT + 400) / 100
SS$MSP <- (SS$MSP - 400) / 200

SS$DD02 <- SS$DD0^2
SS$DD52 <- SS$DD5^2

SS$CMD <- NULL
#CMI
#CMIJJA
SS$DD01 <- NULL
SS$DD51 <- NULL
#EMT
SS$FFP <- NULL
SS$ID <- NULL
SS$MAP <- NULL
SS$MAT <- NULL
SS$MCMT <- NULL
#MSP
SS$MWMT <- NULL
SS$NFFD <- NULL
SS$PAS <- NULL
SS$PET <- NULL
SS$PPT_sm <- NULL
SS$PPT_wt <- NULL
#TD

## lat/lon
SS$xlon <- (SS$Xcl - 285400) / 1500000
SS$xlat <- (SS$Ycl - 1320000) / 740000
SS$xlon2 <- SS$xlon^2
SS$xlat2 <- SS$xlat^2

SS$FIRE_HA <- NULL
SS$YearLoss[SS$YearLoss == 0] <- NA

SS$YearFire[is.na(SS$YearFire)] <- 2016 - 200
SS$YearLoss[is.na(SS$YearLoss)] <- 2016 - 200

XYSS <- as.matrix(SS[,c("X","Y","Xcl","Ycl")])
rownames(XYSS) <- SS$SS
ii <- intersect(rownames(YYSS), rownames(XYSS))
YYSS <- YYSS[ii,]
XYSS <- XYSS[ii,]
rm(ii)

DAT <- data.frame(PKEY, SS[match(PKEY$SS, SS$SS),])
DAT$SS.1 <- NULL
DAT$PCODE.1 <- NULL
rownames(DAT) <- DAT$PKEY

## years since fire
DAT$YSF <- DAT$YEAR - DAT$YearFire
DAT$YSF[DAT$YSF < 0] <- 100
## years since loss
DAT$YSL <- DAT$YEAR - DAT$YearLoss
DAT$YSL[DAT$YSL < 0] <- 100
## years since most recent burn or loss
DAT$YSD <- pmin(DAT$YSF, DAT$YSL)

DAT$BRN <- ifelse(DAT$YSF <= 10, 1L, 0L)
DAT$LSS <- ifelse(DAT$YSL <= 10, 1L, 0L)
DAT$LSS[DAT$YEAR < 2000] <- NA
DAT$DTB <- ifelse(DAT$YSD <= 10, 1L, 0L)
DAT$DTB[DAT$YEAR < 2000] <- NA


## decid + mixed
DAT$isDM <- ifelse(DAT$HAB_NALC2 %in% c("Decid", "Mixed"), 1L, 0L)
## non-forest (wet etc)
DAT$isNF <- ifelse(DAT$HAB_NALC1 %in%
    c("Agr", "Barren", "Devel", "Grass", "Shrub",
    "WetOpen", "DecidOpen", "ConifOpen", "MixedOpen"), 1L, 0L)
DAT$isDev <- ifelse(DAT$HAB_NALC2 %in% c("Agr", "Devel"), 1L, 0L)
DAT$isOpn <- DAT$isNF
DAT$isOpn[DAT$isDev == 1] <- 0
DAT$isWet <- ifelse(DAT$HAB_NALC2 %in% c("Wet"), 1L, 0L)
DAT$isDec <- ifelse(DAT$HAB_NALC2 %in% c("Decid"), 1L, 0L)
DAT$isMix <- ifelse(DAT$HAB_NALC2 %in% c("Mixed"), 1L, 0L)

## see above for filtering
DAT <- DAT[DAT$YEAR >= 1997,]

## year effect
#plot(table(DAT$YEAR))
#DAT$YR <- DAT$YEAR - 1997

#keep <- rep(TRUE, nrow(DAT))
#keep[DAT$BCR %in% c(1,5,11,17,22,24,26,28,29,30)] <- FALSE
#keep[DAT$BCR %in% c(9, 10)] <- FALSE
#keep[DAT$BCR %in% c(9, 10) & DAT$JURS %in% c("AB","BC")] <- TRUE
#keep[DAT$BOREALLOC != "OUT"] <- TRUE
#keep[is.na(DAT$BCR) | DAT$BCR == "0"] <- FALSE

keep <- !is.na(DAT$BOREALLOC) & DAT$BOREALLOC != "OUT"
keep[!is.null(DAT$BCR) & DAT$BCR %in% c("12", "13","14")] <- TRUE

table(DAT$BCR,keep)

DAT <- droplevels(DAT[keep,])

## fill-in BCR based on closes known value
length(which(DAT$BCR == "0" | is.na(DAT$BCR)))
ssDAT <- nonDuplicated(DAT, SS)
ssDAT$xBCR <- ssDAT$BCR
ii <- which(ssDAT$BCR == "0" | is.na(ssDAT$BCR))
for (i in ii) {
    d <- sqrt((ssDAT$X - ssDAT$X[i])^2 + (ssDAT$Y - ssDAT$Y[i])^2)
    d[ii] <- Inf
    ssDAT$xBCR[i] <- ssDAT$BCR[which.min(d)]
}
DAT$xBCR <- ssDAT$xBCR[match(DAT$SS, ssDAT$SS)]
length(which(DAT$BCR == "0" | is.na(DAT$BCR)))
length(which(DAT$xBCR == "0" | is.na(DAT$xBCR)))

table(DAT$JURS, DAT$xBCR)
DAT$Units <- factor(NA, c("AK","N","W","8E","Mt","HW", "At"))
DAT$Units[DAT$xBCR %in% c("2","4")] <- "AK"
DAT$Units[DAT$xBCR %in% c("3","7")] <- "N" # Arctic/Shield
DAT$Units[DAT$xBCR %in% c("5","9","10")] <- "Mt" # Rockies

DAT$Units[DAT$xBCR %in% c("6","11")] <- "W" # West
DAT$Units[DAT$xBCR %in% c("8") & DAT$JURS %in% c("MB","SK")] <- "W"

DAT$Units[DAT$xBCR %in% c("8") & !(DAT$JURS %in% c("MB","SK"))] <- "8E" # BCR8 east

DAT$Units[DAT$xBCR %in% c("12","13","23")] <- "HW" # Hardwood
DAT$Units[DAT$xBCR %in% c("14")] <- "At" # Atlantic

table(DAT$xBCR, DAT$Units, useNA="a")
table(DAT$Units, useNA="a")

ssDAT <- nonDuplicated(DAT, SS)
table(ssDAT$Units, useNA="a")


#library(RColorBrewer)
#plot(DAT[,c("Xcl","Ycl")], col=DAT$REG, pch=19, cex=0.2)
#plot(DAT[,c("Xcl","Ycl")], col=brewer.pal(12, "Set3")[DAT$BCR_JURS], pch=19, cex=0.2)
#plot(DAT[,c("Xcl","Ycl")], col=ifelse(as.integer(DAT$BCR_JURS)==8,2,1), pch=19, cex=0.2)
#plot(DAT[,c("Xcl","Ycl")], col=DAT$EW, pch=19, cex=0.2)

## resampling blocks
DAT$YR5 <- 0
DAT$YR5[DAT$YEAR > 2000] <- 1
DAT$YR5[DAT$YEAR > 2004] <- 2
DAT$YR5[DAT$YEAR > 2009] <- 3
table(DAT$YEAR,DAT$YR5)
table(DAT$YR5)
DAT$bootg <- interaction(DAT$Units, DAT$YR5, drop=TRUE)

## exclude same year revisits
#DAT$ROUND[is.na(DAT$ROUND)] <- 1
#DAT <- droplevels(DAT[DAT$ROUND == 1,])
DAT$SS_YR <- interaction(DAT$SS, DAT$YEAR, drop=TRUE)
table(table(DAT$SS_YR))
dup <- unique(DAT$SS_YR[duplicated(DAT$SS_YR)])
tmp1 <- DAT[DAT$SS_YR %in% dup, c("PKEY","SS_YR")]
rownames(tmp1) <- tmp1$PKEY
dim(tmp1)
tmp2 <- nonDuplicated(tmp1[sample.int(nrow(tmp1)),], SS_YR)
rownames(tmp2) <- tmp2$PKEY
dim(tmp2)
dim(tmp1) - dim(tmp2)
tmp1 <- tmp1[setdiff(rownames(tmp1), rownames(tmp2)),]
dim(tmp1)
dim(tmp1) + dim(tmp2)
DAT <- DAT[!(DAT$PKEY %in% rownames(tmp1)),]
DAT$SS_YR <- NULL

## unit of resampling is site x year to maintain temporal distribution
DAT$SITE <- as.character(DAT$SITE)
DAT$SITE[is.na(DAT$SITE)] <- paste0("gr.", as.character(DAT$gridcode[is.na(DAT$SITE)]))
DAT$SITE <- as.factor(DAT$SITE)

DAT$SITE_YR <- interaction(DAT$SITE, DAT$YEAR, drop=TRUE, sep="_")

DAT$ROAD_OK[is.na(DAT$ROAD_OK)] <- 1L

DAT$ARU <- ifelse(DAT$PCODE %in% levels(DAT$PCODE)[grep("EMC", levels(DAT$PCODE))],
    1L, 0L)

stopifnot(all(!is.na(DAT$PKEY)))

## placeholder for distance to nearest detection
#DAT$ND2 <- 0

summary(DAT$YSD)
summary(DAT$YSF)
summary(DAT$YSL)
AGEMAX <- 50
DAT$YSD <- pmax(0, 1 - (DAT$YSD / AGEMAX))
DAT$YSF <- pmax(0, 1 - (DAT$YSF / AGEMAX))
DAT$YSL <- pmax(0, 1 - (DAT$YSL / AGEMAX))

DAT$HAB <- DAT$HAB_NALC2
DAT$HABTR <- DAT$HAB_NALC1
DAT$HGT <- DAT$HGTorig
DAT$HGT[DAT$HAB %in% c("Agr","Barren","Devel","Grass", "Shrub")] <- 0
DAT$HGT2 <- DAT$HGT^2
DAT$HGT05 <- sqrt(DAT$HGT)

DAT <- DAT[DAT$ARU == 0,]

data.frame(n=colSums(is.na(DAT)))


keep <- rep(TRUE, nrow(DAT))
keep[DAT$ROAD_OK < 1] <- FALSE
keep[is.na(DAT$HABTR)] <- FALSE
keep[is.na(DAT$CMI)] <- FALSE
keep[is.na(DAT$HGT)] <- FALSE
keep[is.na(DAT$SLP)] <- FALSE
keep[DAT$YEAR < 2001] <- FALSE # GFW years
keep[DAT$YEAR > 2013] <- FALSE # GFW years

data.frame(n=colSums(is.na(DAT[keep,])))

DAT <- droplevels(DAT[keep,])
DAT$YR <- DAT$YEAR - min(DAT$YEAR)

## check NAs, exclude ARU

library(detect)
source("~/repos/bamanalytics/R/analysis_mods.R")
source("~/repos/bragging/R/glm_skeleton.R")

nmin <- 25
B <- 239
Extra <- c("SS", "SITE_YR", "X", "Y", "Xcl", "Ycl", "Units", "xBCR", "JURS")
Save <- c("DAT", "YY", "OFF", "TAX", "mods", "BB")
Date <- "2016-12-01"

pk <- intersect(rownames(DAT), rownames(YY))
DAT <- DAT[pk,]

set.seed(1234)
## make sure that time intervals are considered as blocks
## keep out 10% of the data for validation
id2 <- list()
for (l in levels(DAT$bootg)) {
    sset <- which(DAT$bootg == l)
    id2[[l]] <- sample(sset, floor(length(sset) * 0.9), FALSE)
}
KEEP_ID <- unname(unlist(id2))
HOLDOUT_ID <- setdiff(seq_len(nrow(DAT)), KEEP_ID)

DATk <- DAT[KEEP_ID,]
DAT <- DAT[c(KEEP_ID, HOLDOUT_ID),]
DATk$SITE <- droplevels(DATk$SITE)
BB <- hbootindex(DATk$SITE, DATk$bootg, B=B)

#pk <- intersect(rownames(DAT), rownames(YY))
pk <- rownames(DAT)
#DAT <- DAT[pk,]
YY <- YY[pk,]
YY <- YY[,colSums(YY) >= nmin]
apply(as.matrix(YY), 2, max)
#apply(as.matrix(YY), 2, table)

OFF <- OFF[pk,colnames(YY)]
TAX <- TAX[colnames(YY),]

DAT$CMI2 <- DAT$CMI^2
DAT$CMIJJA2 <- DAT$CMIJJA^2

DAT <- DAT[,c(Extra, getTerms(mods, "list"))]

save(list = Save,
    file=file.path(ROOT, "out", "data", paste0("pack_", Date, ".Rdata")))



ROOT <- "e:/peter/bam/Apr2016"
load(file=file.path(ROOT, "out", "data", paste0("pack_2016-04-18.Rdata")))

library(rworldmap)
X0 <- nonDuplicated(DAT, SS)
plot(getMap(resolution = "low"),
    xlim = c(-193, -48), ylim = c(38, 72), asp = 1)
points(X0[, c("X","Y")], pch=19, cex=0.25, col=as.integer(X0$Units)+1)

library(lattice)
densityplot(~ I(HGT*25) | HABTR, DAT[DAT$HGT>0,])

