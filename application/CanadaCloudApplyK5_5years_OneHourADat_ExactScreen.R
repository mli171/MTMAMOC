library(VGAM)
library(lubridate)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(data.table)
library(readr)
library(ggplot2)

library(foreach)
library(doMC)

sourceCpp("~/Desktop/MTM/MTMAMOC/src/QuasiNewton.cpp")
source("~/Desktop/MTM/MTMAMOC/tools/rfunc.R")

#-------------------------- Detection single job
one.MTM.LogL <- function(job=NULL, nTasks=NULL, nCore=NULL,
                         y_log=NULL, mm=NULL,
                         tauClc=NULL,DesignXH0=NULL,
                         stepsize=NULL, mytol=NULL){
  Ti <- dim(y_log)[1]
  K  <- dim(y_log)[2]
  len_par <- (K-1)*(dim(DesignXH0)[2]+1)+(K-1)^2
  
  nSubtasks <- nTasks/nCore
  RES <- data.frame(job=job, task=1:nSubtasks, tim=rep(NA, nSubtasks),
                    logLHa=rep(NA, nSubtasks), tau=rep(NA, nSubtasks))
  ParamMatrix <- matrix(0, nrow=nSubtasks, ncol=len_par)
  RES <- cbind(RES, ParamMatrix)
  
  #################################################################
  #               Parameter Estimation under Ha                   #
  #################################################################
  for(subid in 1:nSubtasks){
    
    tim1 <- Sys.time()
    tauEst     <- tauClc[(job-1)*nSubtasks+subid]
    CpValue    <- c(rep(0, tauEst-1), rep(1, Ti-tauEst+1))
    DesignXHa  <- cbind(DesignXH0, CpValue)
    parinitialHa  <- ParInitCalcu(y_log=y_log, m=mm, DesignX=DesignXHa)
    resEstHa   <- SeqQuasiNWCpp(parinitialHa, y_log, mm, DesignXHa, 
                                stepsize, 1000, mytol)
    parHa      <- resEstHa$parEst
    logLHa     <- LoglikCalcCpp(parHa, DesignXHa, y_log, mm, 
                                conv = mytol, stepNR=stepsize[2])
    tim2 <- Sys.time()
    timecost <- difftime(tim2, tim1, units="mins")
    cat("\n No.job", job, "task", subid, 
        "within", nSubtasks, "(subtasks) completed within", 
        timecost, "mins!")
    #################################################################
    #---------------------- Result Return --------------------------#
    #################################################################
    RES$tim[subid]            <- timecost
    RES$logLHa[subid]         <- logLHa
    RES$tau[subid]            <- tauEst
    RES[subid, 6:(5+len_par)] <- parHa
  }
  
  return(RES)
}

#-------------------------- Paraellel computing implementation
useMultiCore <- function(nTasks=NULL, nCore=NULL,
                         y_log=NULL, mm=NULL,
                         tauClc=NULL, DesignXH0=NULL,
                         stepsize=NULL, mytol=NULL){
  cat("Multicores working, please wait ... \n")
  registerDoMC(cores = nCore)
  tim.start = Sys.time()
  FinalRes <- foreach(i=1:nCore, .combine = "rbind") %dopar%
    one.MTM.LogL(job=i, nTasks=nTasks, nCore=nCore,
                 y_log=y_log, mm=mm, tauClc=tauClc, 
                 DesignXH0=DesignXH0, stepsize=stepsize, 
                 mytol=mytol)
  tim.end = Sys.time()
  cat("Done.\n")
  cat("\n\n nTasks =", nTasks, "\t nCore =", nCore,
      "\t Aveg. time =", mean(FinalRes$tim),
      "\t Total Time =", difftime(tim.end, tim.start, units="hours"),
      "hours \n\n")
  return(FinalRes)
}

###################################
##### Main Script Starts Here #####
###################################

datdir = "~/OneDrive - University of Louisiana Lafayette/data/"
dat = read.table(file=paste0(datdir, "CloudCover/FortStJohnAirport/hourly_day.txt"))
dat$Date = ISOdate(year = dat$V1, month = dat$V2, day = dat$V3, hour = dat$V4)

## filter1: 30 years data from 1965 to 1994
dat1 = dat[dat$V1>=1976 & dat$V1 <=1980,]
dat1$V5 = as.numeric(as.character(dat1$V5))

## filter2: to exclude Feb 29 for exact period T=365
premain = which(!(dat1$V2 == 2 & dat1$V3 == 29))
dat1 = dat1[premain, ]

# Computation setting
ss = 365
nCoreuse = 10
mytol    = 1e-05
stepsize = c(1,0.1)

TargetHourClc = c(9, 15)
TargetHourClc = c(15)

for(iii in 1:length(TargetHourClc)){
  
  #---- DownSampling
  TargetHour = TargetHourClc[iii]
  dat2 = dat1[dat1$V4==TargetHour, ]
  dat2$V6 = catg_agg(dat2$V5)
  mm = dat2$V6
  
  Ti    = length(mm)
  K     = max(mm) + 1
  y_log = matrix(0, nrow=Ti, ncol=K)
  for (tt in 1:Ti) {y_log[tt,mm[tt]+1] = 1}
  
  #---- Candidate changepoint set
  lc     = floor(0.05*Ti)
  uc     = floor(0.95*Ti)
  tmp    = lc:uc
  tauClc = lc:(lc+(floor(length(tmp)/nCoreuse) + 1)*nCoreuse-1)
  nTasks = length(tauClc)
  
  #---- Design matrix
  AValue     = rep(1,Ti)
  TrendValue = 1:Ti/Ti
  BValue     = cos(2*pi*(1:Ti)/ss)
  DValue     = sin(2*pi*(1:Ti)/ss)
  
  # #---- Parallel computing (Model 1)
  # DesignXH0  = cbind(AValue, BValue, DValue)
  
  #---- Parallel computing (Model 2)
  DesignXH0  = cbind(AValue, TrendValue, BValue, DValue)
  
  CPCloudHaClc = useMultiCore(nTasks=nTasks, nCore=nCoreuse,
                              y_log=y_log, mm=mm,
                              tauClc=tauClc, DesignXH0=DesignXH0,
                              stepsize=stepsize, mytol=mytol)
  
  # save(CPCloudHaClc, file=paste0("~/Desktop/MTM/MTMAMOC/application/res/MTMApplication_CandaCloudK5_SEA_5Years76to80_Hour", TargetHour,".RData")) # (Model 1)
  save(CPCloudHaClc, file=paste0("~/Desktop/MTM/MTMAMOC/application/res/MTMApplication_CandaCloudK5_FULL_5Years76to80_Hour", TargetHour,".RData")) # (Model 2)
}

