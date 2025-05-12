library(MASS)
library(stats)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(lmtest)
library(stringr)
library(lubridate)
library(ggplot2)
library(forecast)
###read nfdrs data

setwd("G:\\Workspace\\ercm")
stationcat <- read.csv("virtualStnCatalog_update.csv")
statnum <- stationcat[,2]
statnum.0 <- str_pad(statnum, 6, pad = "0")
stationcat[,2] <- statnum.0
write.table(stationcat, "virtualStnCatalog_update_zero.csv", row.names = FALSE, sep = ",", quote= F)
setwd("G:\\Workspace\\ercm\\NFDRS_output")

flist <- Sys.glob("*.csv")
flist <- flist[grep("NFDRS", flist)]

nstations <- length(flist)

yearseq <- 2006:2020

ndays.all <- 365*length(yearseq)

ercmat <- matrix(NA, nrow = ndays.all, ncol = nstations)
namevec <- rep(NA, nstations)
ercmat.base <- ercmat

dayseq.all <- rep(1:365, length(yearseq))
yearseq.all <- rep(yearseq, each = 365)

for(curstation in 1:nstations)
{
  curdat <- read.csv(flist[curstation])
  ##subset to 2004:2018
  
  datevals <- curdat$Date
  datev <- as.Date(datevals)
  yearv <- format(datev, "%Y")
  yearn <- as.numeric(yearv)
  subyears <- yearn %in% yearseq
  curdat <- curdat[subyears,]
  tempdate <- curdat$Date
  doy <- strftime(tempdate, format = "%j")
  
  #ditch leap days
  curdat.sub <- curdat[doy != "366",]
  
  sumnas <- sum(is.na(doy))
  # testval <- dim(curdat.sub)[1] * (sumnas + 1)
  # if(testval == 10220){    
  #   ercmat[,curstation] <- curdat.sub$ERC
  #   tempname <- flist[curstation]
  #   tempname <- gsub(".csv", "", tempname)
  #   namevec[curstation] <- tempname  }
  ercmat[,curstation] <- curdat.sub$ERC
  
  tempname <- gsub(".csv", "", flist[curstation])
  namevec[curstation] <- tempname
  print(curstation/nstations)
  
  
  
}

ercmat.base <- ercmat



dailymeanmat.base <- matrix(NA, nrow = 365, ncol = nstations)
dailymeanmat.rcp8 <- matrix(NA, nrow = 365, ncol = nstations)

dailysdmat <- matrix(NA, nrow = 365, ncol = nstations)
residsmat.base <- matrix(NA, nrow = nrow(ercmat.base), ncol = ncol(ercmat.base))
residsmat.rcp8 <- matrix(NA, nrow = nrow(ercmat.base), ncol = ncol(ercmat.base))
residsmat.filt.base <- matrix(NA, nrow = nrow(ercmat.base), ncol = ncol(ercmat.base))
residsmat.filt.rcp8 <- matrix(NA, nrow = nrow(ercmat.base), ncol = ncol(ercmat.base))

armod.list.base <- vector("list", nstations)
mamod.list.base <- vector("list", nstations)

armod.list.rcp8 <- vector("list", nstations)

sdvec.resids.base <- rep(NA, nstations)
sdvec.filt.base <- rep(NA, nstations)

sdvec.resids.rcp8 <- rep(NA, nstations)
sdvec.filt.rcp8 <- rep(NA, nstations)


for(curstation in 1:nstations)
{
  tempvec <- ercmat.base[,curstation]
  tempvec <- na.exclude(tempvec)
  tempsd <- sd(tempvec)
  ercdf <- data.frame(erc = ercmat.base[,curstation], year = as.character(yearseq.all), day = (dayseq.all))
  ercdf <- na.exclude(ercdf)
  
  if(tempsd > 1e-16)
  {
    ercsummary <- ercdf %>% 
      group_by(day) %>% 
      summarize(mean = mean(erc),
                sd  = sd(erc),
                perc80 = quantile(erc, .8),
                max = max(erc)
                
      )
    
    # 
    # plot(ercsummary$day, ercsummary$mean)
    # plot(ercsummary$day, ercsummary$sd)
    # 
    meanvec <- rep(ercsummary$mean, length(yearseq))
    dailymeanmat.base[,curstation]  <- ercsummary$mean
    ##get smoothed spline
    
    #erc.spline <- smooth.spline(ercdf$erc ~ercdf$day, df = 20)
    # plot(ercsummary$day, ercsummary$mean)
    #lines(erc.spline, col = "red")
    residvec.base <- ercmat.base[,curstation] - meanvec
    armod.base <- ar(residvec.base, order.max = 20, method = "yule-walker")
    residsmat.filt.base[,curstation] <- armod.base$resid
    armod.list.base[[curstation]] <- armod.base
    
    residsmat.base[,curstation] <- residvec.base
    
    
    sdvec.resids.base[curstation] <- sd(residvec.base)
    sdvec.filt.base[curstation] <- sd(na.exclude(armod.base$resid))
    # 
    # residvec.rcp8 <- ercmat.rcp8[,curstation] - meanvec
    # armod.rcp8 <- ar(residvec.rcp8, order.max = 20, method = "yule-walker")
    # armod.list.rcp8[[curstation]] <- armod.rcp8
    # 
    # sdvec.resids.rcp8[curstation] <- sd(residvec.rcp8)
    # sdvec.filt.rcp8[curstation] <- sd(na.exclude(armod.rcp8$resid))
    # 
    # ercdf <- data.frame(erc = ercmat.rcp8[,curstation], year = as.character(yearseq.all), day = as.character(dayseq.all))
    # ercdf <- na.exclude(ercdf)
    # 
    #   ercsummary <- ercdf %>% 
    #     group_by(day) %>% 
    #     summarize(mean = mean(erc),
    #               sd  = sd(erc),
    #               perc80 = quantile(erc, .8),
    #               max = max(erc)
    #               
    #     )
    # 
    # # 
    # # plot(ercsummary$day, ercsummary$mean)
    # # plot(ercsummary$day, ercsummary$sd)
    # # 
    # meanvec <- rep(ercsummary$mean, 15)
    # dailymeanmat.rcp8[,curstation]  <- ercsummary$mean
    # 
    # residvec.rcp8 <- ercmat.rcp8[,curstation] - meanvec
    # residsmat.rcp8[,curstation] <- residvec.rcp8
    # armod.rcp8 <- ar(residvec.rcp8, 30)
    # armod.list.rcp8[[curstation]] <- armod.rcp8
    # residsmat.filt.rcp8[,curstation] <- armod.rcp8$resid
    # 
  }
  print(curstation / nstations)
}




validmat.filt.base <- data.frame(residsmat.filt.base)

validmat.filt.base[is.na(validmat.filt.base)] <- 0

# validmat.filt.rcp8 <- data.frame(residsmat.filt.rcp8)

 numnas.base <- rep(NA, dim(validmat.filt.base)[1])
# numnas.rcp8 <- rep(NA, dim(validmat.filt.rcp8)[1])
for(i in 1:length(numnas.base))
{
  numnas.base[i] <- sum(is.na(validmat.filt.base[i,]))
  # numnas.rcp8[i] <- sum(is.na(validmat.filt.rcp8[i,]))
  
  print(i / length(numnas.base))
}

# validmat.filt.base <- validmat.filt.base[numnas.base < 57,]

# validmat.filt.rcp8 <- validmat.filt.rcp8[numnas.rcp8 < 57,]



validmat.filt.base %>% drop_na()
#              gene hsap mmul mmus rnor cfam
# 2 ENSG00000199674    0    2    2    2    2
# 6 ENSG00000221312    0    1    2    3    2
colsd.base <- apply(validmat.filt.base, 2, "sd")
#colsd.rcp8 <- apply(validmat.filt.rcp8, 2, "sd")

colsd.base[is.na(colsd.base)] <- 0
# colsd.rcp8[is.na(colsd.rcp8)] <- 0

nonzero.base <- validmat.filt.base[,colsd.base > 1e-16]
# nonzero.rcp8 <- validmat.filt.rcp8[,colsd.base > 1e-16]

nonzero.ercmean.base <- dailymeanmat.base[,colsd.base > 1e-16]
# nonzero.ercmean.rcp8<- dailymeanmat.rcp8[,colsd.rcp8 > 1e-16]
namevec.sub <- namevec[colsd.base > 1e-16]

armod.valid.base <- armod.list.base[colsd.base > 1e-16]

# armod.valid.rcp8 <- armod.list.rcp8[colsd.base > 1e-16]

covmat.base <- cov(nonzero.base)
# covmat.rcp8 <- cov(nonzero.rcp8)

cormat.base <- cor(nonzero.base)
# cormat.rcp8 <- cor(nonzero.rcp8)

melted_cormat.base <- melt(cormat.base)
head(melted_cormat.base)


nyear.sim <- 1001
n.sim <- nyear.sim * 365
nstations.sim <- dim(covmat.base)[1]
t1 <- Sys.time()

residsim.base <- mvrnorm(n = n.sim, rep(0, nstations.sim), covmat.base)
#residsim.rcp8 <- mvrnorm(n = n.sim, rep(0, nstations.sim), covmat.rcp8)
t2 <- Sys.time()

t2 - t1


###Filter correlated residuals

validstations <- dim(residsim.base)[2]

simresids.base <- matrix(NA, nrow = nrow(residsim.base), ncol = ncol(residsim.base))
# simresids.rcp8 <- matrix(NA, nrow = nrow(residsim.base), ncol = ncol(residsim.base))
output.ercmean.base <-  matrix(NA, nrow = nrow(residsim.base), ncol = ncol(residsim.base))
# output.ercmean.rcp8 <- matrix(NA, nrow = nrow(residsim.base), ncol = ncol(residsim.base))

for(curstation in 1:validstations)
{
  armod.base <- armod.valid.base[[curstation]] 
  
  sr.b <- sdvec.resids.base[curstation]
  sf.b <- sdvec.filt.base[curstation]

   sr.r <- sdvec.resids.rcp8[curstation]
   sf.r <- sdvec.filt.rcp8[curstation]
  
  r.b <- sr.b / sf.b
  r.r <- sr.r / sf.r
  
  curfilt.base <- stats::filter(residsim.base[,curstation], armod.base$ar, sides = 1, method = "recursive")

  var.est <- armod.base$var.pred
  sd.est <- sqrt(var.est)
  curfilt.sd <- sd(na.exclude(curfilt.base))
  rtemp <- sf.b / curfilt.sd
  curfilt.base <- curfilt.base / curfilt.sd
  curfilt.base <- curfilt.base * sd.est ###To ensure proper variance of ouptut residuals
  
  simresids.base[,curstation] <- curfilt.base
  outvec.mean <- dailymeanmat.base[,curstation]
  outvec.mean <- rep(outvec.mean,nyear.sim)
  outvec.erc <- outvec.mean + curfilt.base
  outvec.erc[outvec.erc < 0] <- 0
  output.ercmean.base[,curstation] <- outvec.erc
  
  
  # 
  # armod.rcp8 <- armod.valid.rcp8[[curstation]] 
  # 
  # curfilt.rcp8 <- stats::filter(residsim.rcp8[,curstation], armod.rcp8$ar, sides = 1, method = "recursive")
  # curfilt.sd <- sd(na.exclude(curfilt.rcp8))
  # rtemp <- sf.r / curfilt.sd
  # 
  # curfilt.rcp8 <- curfilt.rcp8 * rtemp
  # curfilt.rcp8 <- curfilt.rcp8 * r.r
  # 
  # simresids.rcp8[,curstation] <- curfilt.rcp8
  # outvec.mean <- dailymeanmat.rcp8[,curstation]
  # outvec.mean <- rep(outvec.mean,nyear.sim)
  # outvec.erc <- outvec.mean + curfilt.rcp8
  # outvec.erc[outvec.erc < 0] <- 0
  # output.ercmean.rcp8[,curstation] <- outvec.erc
  # 
  print(curstation / validstations)
}

colnames(output.ercmean.base) <- namevec




