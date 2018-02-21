ROOT <- "e:/peter/bam/Apr2016"
library(mefa4)
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
SS$CMI2 <- SS$CMI^2
SS$CMIJJA2 <- SS$CMIJJA^2

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

#keep <- !is.na(DAT$BOREALLOC) & DAT$BOREALLOC != "OUT"
#keep[!is.null(DAT$BCR) & DAT$BCR %in% c("12", "13","14")] <- TRUE

DAT$isBBS <- startsWith(rownames(DAT), "BBS")
table(DAT$isBBS, DAT$ROAD)
DAT$JBCR <- interaction(DAT$JURS, DAT$BCR, sep="::", drop=TRUE)
DAT$RoadBBS <- interaction(DAT$ROAD, DAT$isBBS, sep="::", drop=TRUE)

JB <- c(
    #"AK::5", "BC::5",
    #"BC::9",
    #"AB::10", "BC::10",
    #"MN::23",
    "AK::4", "BC::4", "YK::4",
    "AB::6", "BC::6", "MB::6", "NT::6", "SK::6",
    "MB::8", "NL::8", "ON::8", "QC::8",
    "MB::12", "MN::12", "ON::12", "QC::12", "WI::12",
    "ON::13", "QC::13",
    "NB::14", "NS::14", "PEI::14", "QC::14")

if (FALSE) {
tab_fun <- function(i) {
    j <- DAT$JBCR == i & DAT$RoadBBS != "1::FALSE"
    with(DAT[j,], table(NALC=HAB_NALC2, BBS=isBBS))
}
dets <- lapply(JB, tab_fun)
names(dets) <- JB
dets

ss <- !is.na(DAT$HAB_NALC1) & DAT$ROAD == 0 & DAT$JBCR %in% JB
dat <- DAT[ss,]
rn <- intersect(rownames(dat), rownames(YY))
yy <- ifelse(as.matrix(YY[rn,])>0, 1, 0)
dat <- dat[rn,]
yy <- yy[,colSums(yy)>=20]
yyy <- yy[sample(nrow(yy), 5000),]
yyy <- yyy[,colSums(yyy)>=100]
off <- OFF[rownames(yyy), colnames(yyy)]

library(opticut)
oc <- opticut(yyy ~ 1, strata=dat[rownames(yyy), "HAB_NALC2"], dist="binomial")
plot(hclust(vegan::vegdist(t(summary(oc)$bestpart), "jaccard")), hang=-1, sub="", xlab="")
plot(oc, cex.axis=0.75)
1-round(vegan::vegdist(t(summary(oc)$bestpart), "jaccard"),3)

ftable(DAT$JBCR, DAT$RoadBBS)
ftable(DAT$JBCR, DAT$RoadBBS, DAT$HAB_NALC1)

spp <- "OVEN"
ol <- optilevels(y=yyy[,spp], x=dat[rownames(yyy), "HAB_NALC1"],
    dist="poisson", offset=off[,spp])

ol <- optilevels(y=yyy[,spp], x=dat[rownames(yyy), "HAB_NALC2"],
    dist="poisson", offset=off[,spp])
round(data.frame(D=sort(exp(coef(bestmodel(ol))))),4)
}

library(opticut)
mm <- Mefa(YY, DAT, TAX, "inner")
mm <- mm[!is.na(samp(mm)$RoadBBS) & !is.na(samp(mm)$JBCR) & !is.na(samp(mm)$HAB_NALC1),]
#mm <- mm[,colSums(xtab(mm)>0) >= 100]

tab_data <- function(spp, region) {
    r <- samp(mm)$RoadBBS %in% c("1::TRUE", "0::FALSE")
    mmi <- mm[samp(mm)$JBCR == region & r,]
    y <- factor(ifelse(as.numeric(xtab(mmi)[,spp]) > 0, 1, 0), 0:1)
    x <- samp(mmi)$HAB_NALC1
    list(off=table(NALC=samp(mmi)$HAB_NALC1[samp(mmi)$ROAD==0], Det=y[samp(mmi)$ROAD==0]),
        on=table(NALC=samp(mmi)$HAB_NALC1[samp(mmi)$ROAD==1], Det=y[samp(mmi)$ROAD==1]))
}
get_subset <- function(region, road=FALSE) {
    r <- if (road)
        samp(mm)$RoadBBS == "1::TRUE" else samp(mm)$RoadBBS == "0::FALSE"
    mm[samp(mm)$JBCR == region & r,]
}
find_levels <- function(spp, region, road=FALSE, m=1000) {
    mmi <- get_subset(region, road)
    j <- rep(FALSE, nrow(mmi))
    for (k in levels(droplevels(samp(mm)$HAB_NALC1))) {
        w <- which(samp(mmi)$HAB_NALC1 == k)
        if (length(w) < m) {
            j[w] <- TRUE
        } else {
            j[sample(w, m)] <- TRUE
        }
    }
    mmi <- mmi[j,]
    y <- as.numeric(xtab(mmi)[,spp])
    x <- droplevels(samp(mmi)$HAB_NALC1)
    ol <- optilevels(y=y, x=x, dist="poisson", offset=OFF[rownames(mmi),spp])
    ol
}
'logLik.try-error' <- function (object, ...) {
    structure(-.Machine$double.xmax^(1/3), df = 1,
        nobs = 1, class = "logLik")
}
logLik.glm_skeleton <- function (object, ...) {
    structure(object$logLik, df = object$df,
        nobs = object$nobs, class = "logLik")
}
glm_skeleton <- function(object, ..., CAICalpha=0.5, keep_call=TRUE, vcov=FALSE) {
    if (inherits(object, "try-error"))
        return(structure(as.character(object), class="try-error"))
    out <- structure(list(
        call=object$call,
        formula=formula(object),
        coef=coef(object),
        vcov=NULL,
        converge=object$converge,
        logLik=as.numeric(logLik(object)),
        df=attr(logLik(object), "df"),
        nobs=nobs(object)), class="glm_skeleton")
    if (!out$converge)
        return(structure("glm did not converge", class="try-error"))
    if (!keep_call) {
        out$call <- out$formula <- NULL
    }
    if (vcov)
        out$vcov <- vcov(object)
    out$class0 <- class(object)[1L]
    out$aic <- -2*out$logLik + 2*out$df
    out$bic <- -2*out$logLik + log(out$nobs)*out$df
    out$caic <- CAICalpha * out$aic + (1-CAICalpha) * out$bic
    out
}

## null model
model_null <- function(spp, region, road=FALSE, trend=FALSE) {
    ff0 <- y ~ 1
    ff1 <- y ~ yr
    ff <- if (trend)
        ff1 else ff0
    mmi <- get_subset(region, road)
    y <- as.numeric(xtab(mmi)[,spp])
    yr <- samp(mmi)$YEAR
    glm_skeleton(try(
        glm(ff, family="poisson", offset=OFF[rownames(mmi),spp])),
        keep_call=FALSE, vcov=TRUE)
}
model_lcc <- function(spp, region, road=FALSE, reclass=NULL, trend=FALSE) {
    ff0 <- y ~ x
    ff1 <- y ~ yr + x
    ff <- if (trend)
        ff1 else ff0
    mmi <- get_subset(region, road)
    y <- as.numeric(xtab(mmi)[,spp])
    x <- droplevels(samp(mmi)$HAB_NALC1)
    if (!is.null(reclass))
        levels(x) <- reclass[levels(x)]
    yr <- samp(mmi)$YEAR
    glm_skeleton(try(
        glm(ff, family="poisson", offset=OFF[rownames(mmi),spp])),
        keep_call=FALSE, vcov=TRUE)
}
model_clim <- function(spp, region, road=FALSE, trend=FALSE) {
    ff0 <- y ~ CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 +
        CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT
    ff1 <- y ~ yr + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 +
        CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT
    ff <- if (trend)
        ff1 else ff0
    mmi <- get_subset(region, road)
    y <- as.numeric(xtab(mmi)[,spp])
    yr <- samp(mmi)$YEAR
    glm_skeleton(try(glm(ff, data=samp(mmi),
        family="poisson", offset=OFF[rownames(mmi),spp])),
        keep_call=FALSE, vcov=TRUE)
}
model_lcclim <- function(spp, region, road=FALSE, reclass=NULL, trend=FALSE) {
    ff0 <- y ~ x + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 +
        CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT
    ff1 <- y ~ yr + x + CMI + CMIJJA + DD0 + DD5 + EMT + MSP + TD + DD02 + DD52 + CMI2 + CMIJJA2 +
        CMIJJA:DD0 + CMIJJA:DD5 + EMT:MSP + CMI:DD0 + CMI:DD5 + MSP:TD + MSP:EMT
    ff <- if (trend)
        ff1 else ff0
    mmi <- get_subset(region, road)
    y <- as.numeric(xtab(mmi)[,spp])
    x <- droplevels(samp(mmi)$HAB_NALC1)
    if (!is.null(reclass))
        levels(x) <- reclass[levels(x)]
    yr <- samp(mmi)$YEAR
    glm_skeleton(try(glm(ff, data=samp(mmi),
        family="poisson", offset=OFF[rownames(mmi),spp])),
        keep_call=FALSE, vcov=TRUE)
}

model_all <- function(region, spp) {
    t0 <- proc.time()
    tab <- tab_data(spp, region)
    Err <- structure("0 offroad detections", class="try-error")
    if (sum(tab$off[,"1"]) > 0) {
        set.seed(1)
        ol <- find_levels(spp, region, road=FALSE, m=1000) # use subset of offroad data
        rc <- ol$levels[[length(ol$levels)]]
        D00 <- model_null(spp, region, road=FALSE)
        D01 <- model_lcc(spp, region, road=FALSE, rc)
        D02 <- model_clim(spp, region, road=FALSE)
        D03 <- model_lcclim(spp, region, road=FALSE, rc)
        T00 <- model_null(spp, region, road=FALSE, trend=TRUE)
        T01 <- model_lcc(spp, region, road=FALSE, reclass=rc, trend=TRUE)
        T02 <- model_clim(spp, region, road=FALSE, trend=TRUE)
        T03 <- model_lcclim(spp, region, road=FALSE, reclass=rc, trend=TRUE)

        if (sum(tab$on[,"1"]) > 0) {
            D10 <- model_null(spp, region, road=TRUE)
            D11 <- model_lcc(spp, region, road=TRUE, reclass=rc)
            D12 <- model_clim(spp, region, road=TRUE)
            D13 <- model_lcclim(spp, region, road=TRUE, reclass=rc)
            T10 <- model_null(spp, region, road=TRUE, trend=TRUE)
            T11 <- model_lcc(spp, region, road=TRUE, reclass=rc, trend=TRUE)
            T12 <- model_clim(spp, region, road=TRUE, trend=TRUE)
            T13 <- model_lcclim(spp, region, road=TRUE, reclass=rc, trend=TRUE)
        } else {
            D10 <- D11 <- D12 <- D13 <- Err
            T10 <- T11 <- T12 <- T13 <- Err
        }
    } else {
        D00 <- D01 <- D02 <- D03 <- Err
        T00 <- T01 <- T02 <- T03 <- Err
        if (sum(tab$on[,"1"]) > 0) {
            set.seed(1)
            ol <- find_levels(spp, region, road=TRUE, m=1000) # use subset of offroad data
            rc <- ol$levels[[length(ol$levels)]]
            D10 <- model_null(spp, region, road=TRUE)
            D11 <- model_lcc(spp, region, road=TRUE, reclass=rc)
            D12 <- model_clim(spp, region, road=TRUE)
            D13 <- model_lcclim(spp, region, road=TRUE, reclass=rc)
            T10 <- model_null(spp, region, road=TRUE, trend=TRUE)
            T11 <- model_lcc(spp, region, road=TRUE, reclass=rc, trend=TRUE)
            T12 <- model_clim(spp, region, road=TRUE, trend=TRUE)
            T13 <- model_lcclim(spp, region, road=TRUE, reclass=rc, trend=TRUE)
        } else {
            ol <- NULL
            D10 <- D11 <- D12 <- D13 <- Err
            T10 <- T11 <- T12 <- T13 <- Err
        }
    }
    out <- list(
        species=spp,
        region=region,
        levels=ol,
        table=tab,
        time=as.numeric(proc.time() - t0)[3L],
        density=list(
            null_off=D00,
            lcc_off=D01,
            clim_off=D02,
            lcclim_off=D03,
            null_on=D10,
            lcc_on=D11,
            clim_on=D12,
            lcclim_on=D13),
        trend=list(
            null_off=T00,
            lcc_off=T01,
            clim_off=T02,
            lcclim_off=T03,
            null_on=T10,
            lcc_on=T11,
            clim_on=T12,
            lcclim_on=T13))
    out
}

#SPP2 <- c("ALFL", "OVEN", "CAWA", "OSFL")
for (spp in SPP) {
    cat(spp, "-------------\n")
    res1 <- list()
    for (v in JB) {
        cat(spp, v, "\n");flush.console()
        res1[[v]] <- model_all(v, spp)
    }
    save(res1, file=paste0("e:/peter/bam/2017/foam/foam-results_", spp, ".Rdata"))
}

res1[["AB::6"]]
h <- function(x) {
    if (inherits(x, "try-error"))
        return(0)
    data.frame(D=round(sort(exp(c(x$coef[1], x$coef[1]+x$coef[-1]))), 4))
}

xx <- lapply(res1, function(z) h(z$density$lcc_on))

## todo:
## OK - add year effect
## - calculate geographic discrepancies
## - compare D_null and T_null with geographic sampling
## - calculate fit in other regions
## - process prediction grid
## - check spatial patterns and change climate if needed
## - pl/cl???

fl <- list.files("e:/peter/bam/2017/foam")
SPP <- substr(sapply(strsplit(fl, "_"), "[[", 2), 1, 4)

LEV <- c("ConifDense", "Agr", "ConifOpen", "ConifSparse", "DecidDense",
    "DecidOpen", "DecidSparse", "Devel", "Grass", "MixedDense", "MixedOpen",
    "MixedSparse", "Shrub", "WetDense", "WetOpen", "WetSparse")
h2 <- function(x) {
    if (inherits(x, "try-error"))
        return(structure(rep(NA, length(LEV)), names=LEV))
    x <- data.frame(D=round(sort(exp(c(x$coef[1], x$coef[1]+x$coef[-1]))), 4))
    out <- list()
    for (i in 1:length(x$D)) {
        a <- rownames(x)[i]
        a <- substr(a, 2, nchar(a))
        a <- strsplit(a, "\\+")[[1]]
        out[[i]] <- structure(rep(x$D[i], length(a)), names=a)
    }
    out <- unlist(out)
    names(out)[names(out) == "Intercept)"] <- LEV[1]
    out <- out[match(LEV, names(out))]
    names(out) <- LEV
    out
}

Den <- list()
#spp <- "CAWA"
for (spp in SPP) {
    fn <- paste0("e:/peter/bam/2017/foam/foam-results_", spp, ".Rdata")
    load(fn)
    res4 <- res1[c("ON::8", "QC::8", "ON::12", "QC::12")]
    Den[[spp]] <- sapply(res4, function(z) h2(z$density$lcc_on))
}
save(Den, file="e:/peter/bam/2017/Density-ON-QC.RData")

f <- function(x) {
    x <- t(x)
    x / apply(x, 1, max, na.rm=TRUE)
}
i <- "CAWA"
barplot(f(Den[[i]]), main=i, beside=TRUE, legend.text=TRUE,
    col=c("tomato", "gold", "grey", "turquoise"), ylab="density / max density")

barplot(t(Den[[i]]), main=names(Den)[i], beside=TRUE)
