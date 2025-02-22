rm(list = ls())

library(VGAM)
library(lubridate)
library(TSA)
library(mvtnorm)
library(Rcpp)
library(RcppArmadillo)
library(data.table)
library(readr)
library(ggplot2)
library(car)
# change the opacity of lines in plot
library(scales)

sourceCpp("src/QuasiNewton.cpp")
source("tools/rfunc.R")

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

datdir = "~/OneDrive - University of Louisiana Lafayette/data/"
dat = read.table(file=paste0(datdir, "CloudCover/FortStJohnAirport/hourly_day.txt"))
dat$Date = ISOdate(year = dat$V1, month = dat$V2, day = dat$V3, hour = dat$V4)

## filter1: 30 years data from 1965 to 1994
dat1 = dat[dat$V1>=1965 & dat$V1 <=1994,]
dat1$V5 = as.numeric(as.character(dat1$V5))

## filter2: to exclude Feb 29 for exact period T=365
premain = which(!(dat1$V2 == 2 & dat1$V3 == 29))
dat1 = dat1[premain, ]

## filter3: to ONLY include the target hour

TargetHour = 9

dat2 = dat1[dat1$V4==TargetHour, ]
dat2$V6 = catg_agg(dat2$V5)


################################################################################
################################################################################
##########                       Model Fitting                        ##########
################################################################################
################################################################################

mm = dat2$V6
m = dat2$V5

Ti     <- length(mm)
K      <- max(mm) + 1
y_log <- matrix(0, nrow=Ti, ncol=K)
for (tt in 1:Ti) {y_log[tt,mm[tt]+1] <- 1}

ss = 365

nCoreuse <- 40
lc  <- floor(0.05*Ti)
uc  <- floor(0.95*Ti)
tmp <- lc:uc
tauClc <- lc:(lc+(floor(length(tmp)/nCoreuse) + 1)*nCoreuse-1)

# Design Matrix
AValue     <- rep(1,Ti)
TrendValue <- 1:Ti/Ti
BValue     <- cos(2*pi*(1:Ti)/ss)
DValue     <- sin(2*pi*(1:Ti)/ss)

# Computation setting
mytol    = 1e-05
stepsize = c(1, 1)
mymaxit  = 1000

#################################################################
#               Parameter Estimation under H0                   #   Model 2 with no changepoint
#################################################################
DesignXEst_M2_H0  <- cbind(AValue, TrendValue, BValue, DValue)
parinitial_M2_H0  <- ParInitCalcu(y_log=y_log, m=mm, DesignX=DesignXEst_M2_H0)
resEst_M2_H0      <- SeqQuasiNWCpp(parinitial_M2_H0, y_log, mm, DesignXEst_M2_H0, stepsize, mymaxit, mytol)
par_est_M2_H0     <- resEst_M2_H0$parEst
mySEtheta_M2_H0   <- sqrt(diag(solve(resEst_M2_H0$Hesstheta)))
mySExi_M2_H0      <- sqrt(diag(solve(resEst_M2_H0$Hessxi)))
mySE_M2_H0        <- c(mySEtheta_M2_H0, mySExi_M2_H0)

parMatrix_M2_H0 <- cbind(par_est_M2_H0[1:length(par_est_M2_H0)],
                         par_est_M2_H0[1:length(par_est_M2_H0)] - 2*mySE_M2_H0,
                         par_est_M2_H0[1:length(par_est_M2_H0)] + 2*mySE_M2_H0)
colnames(parMatrix_M2_H0) <- c("Param", "LowerB", "UpperB")
rownames(parMatrix_M2_H0) <- c(paste("A", 1:(K-1), sep=""),
                               paste("beta", 1:(K-1), sep=""),
                               paste("B", 1:(K-1), sep=""),
                               paste("D", 1:(K-1), sep=""),
                               rep("xi", (K-1)^2))
# round(parMatrix_M2_H0, digits=4)

#------- Fitted Value with VS. without Changepoint for each categories
len_mean_M2_H0 <- dim(DesignXEst_M2_H0)[2]
theta_M2_H0    <- matrix(par_est_M2_H0[1:((K-1)*len_mean_M2_H0)], 
                         nrow=K-1, ncol=len_mean_M2_H0)
eta_M2_H0      <- DesignXEst_M2_H0%*%t(theta_M2_H0)
gamm_M2_H0     <- exp(eta_M2_H0)/(1+exp(eta_M2_H0))
pm_M2_H0       <- matrix(0, Ti, K)
for (i in 1:(K-1)){
  if(i==1){
    pm_M2_H0[,i] <- gamm_M2_H0[,i]
  }else{
    pm_M2_H0[,i] <- gamm_M2_H0[,i] - gamm_M2_H0[,i-1]
  }
}
pm_M2_H0[,K] <- 1 - gamm_M2_H0[,K-1]


logL_M2_H0 <- LoglikCalcCpp(resEst_M2_H0$parEst, DesignXEst_M2_H0, y_log, mm, conv = mytol, stepNR=stepsize[2])
# data from server
load(paste0("application/MTMApplication_CandaCloudK5_FULL_30Years65to94_Hour", TargetHour, ".RData"))
LambdaStat_M2 <- -2*(logL_M2_H0 - CPCloudHaClc$logLHa)
tauest_M2_H0  <- which.max(LambdaStat_M2)


#################################################################
#               Parameter Estimation under Ha                   #   Model 2 with changepoint
#################################################################

tauEst_M2          <- tauClc[tauest_M2_H0]
CpValue            <- c(rep(0, tauEst_M2-1), rep(1, Ti-tauEst_M2+1))
DesignXEst_M2_Ha   <- cbind(DesignXEst_M2_H0, CpValue)
parinitial_M2_Ha   <- ParInitCalcu(y_log=y_log, m=mm, DesignX=DesignXEst_M2_Ha)
resEst_M2_Ha       <- SeqQuasiNWCpp(parinitial_M2_Ha, y_log, mm, DesignXEst_M2_Ha, stepsize, 1000, mytol)
par_est_M2_Ha      <- resEst_M2_Ha$parEst
mySEtheta_M2_Ha    <- sqrt(diag(solve(resEst_M2_Ha$Hesstheta)))
mySExi_M2_Ha       <- sqrt(diag(solve(resEst_M2_Ha$Hessxi)))
mySE_M2_Ha         <- c(mySEtheta_M2_Ha, mySExi_M2_Ha)
parMatrix_M2_Ha   <- cbind(par_est_M2_Ha[1:length(par_est_M2_Ha)],
                           par_est_M2_Ha[1:length(par_est_M2_Ha)] - 2*mySE_M2_Ha,
                           par_est_M2_Ha[1:length(par_est_M2_Ha)] + 2*mySE_M2_Ha)
colnames(parMatrix_M2_Ha) <- c("Param", "LowerB", "UpperB")
rownames(parMatrix_M2_Ha) <- c(paste("A", 1:(K-1), sep=""),
                               paste("beta", 1:(K-1), sep=""),
                               paste("B", 1:(K-1), sep=""),
                               paste("D", 1:(K-1), sep=""),
                               paste("Delta", 1:(K-1), sep=""),
                               rep("xi", (K-1)^2))
# round(parMatrix_M2_Ha, digits=4)


logL_M2_Ha <- LoglikCalcCpp(par_est_M2_Ha, DesignXEst_M2_Ha, y_log, mm, conv = mytol, stepNR=stepsize[2])
LambdaStat.Est_M2_Ha <- -2*(logL_M2_H0 - logL_M2_Ha)
myresid_M2_Ha = MTMResidCpp1(parEst=par_est_M2_Ha, y_log=y_log, mm=mm, DesignX=DesignXEst_M2_Ha, K=K, Ti=Ti, mytol=1e-5, 1)

#------- Fitted Value with VS. without Changepoint for each categories
len_mean_M2_Ha <- dim(DesignXEst_M2_Ha)[2]
theta_M2_Ha    <- matrix(par_est_M2_Ha[1:((K-1)*len_mean_M2_Ha)], 
                         nrow=K-1, ncol=len_mean_M2_Ha)
eta_M2_Ha      <- DesignXEst_M2_Ha%*%t(theta_M2_Ha)
gamm_M2_Ha     <- exp(eta_M2_Ha)/(1+exp(eta_M2_Ha))
pm_M2_Ha       <- matrix(0, Ti, K)
for (i in 1:(K-1)){
  if(i==1){
    pm_M2_Ha[,i] <- gamm_M2_Ha[,i]
  }else{
    pm_M2_Ha[,i] <- gamm_M2_Ha[,i] - gamm_M2_Ha[,i-1]
  }
}
pm_M2_Ha[,K] <- 1 - gamm_M2_Ha[,K-1]


save(dat2, 
     DesignXEst_M2_H0, par_est_M2_H0, parMatrix_M2_H0, pm_M2_H0, gamm_M2_H0, 
     tauClc, LambdaStat_M2, 
     DesignXEst_M2_Ha, par_est_M2_Ha, parMatrix_M2_Ha, pm_M2_Ha, gamm_M2_Ha,
     tauEst_M2, LambdaStat.Est_M2_Ha, myresid_M2_Ha,
     file = paste0("Application/ApplicationNeed_Hour", TargetHour, "_Model2.RData"))


















#################################################################
#               Parameter Estimation under H0                   #   Model 1 with no changepoint
#################################################################

DesignXEst_M1_H0  <- cbind(AValue, BValue, DValue)
parinitial_M1_H0  <- ParInitCalcu(y_log=y_log, m=mm, DesignX=DesignXEst_M1_H0)
resEst_M1_H0      <- SeqQuasiNWCpp(parinitial_M1_H0, y_log, mm, DesignXEst_M1_H0, stepsize, mymaxit, mytol)
par_est_M1_H0     <- resEst_M1_H0$parEst
mySEtheta_M1_H0   <- sqrt(diag(solve(resEst_M1_H0$Hesstheta)))
mySExi_M1_H0      <- sqrt(diag(solve(resEst_M1_H0$Hessxi)))
mySE_M1_H0        <- c(mySEtheta_M1_H0, mySExi_M1_H0)

parMatrix_M1_H0   <- cbind(par_est_M1_H0[1:length(par_est_M1_H0)],
                           par_est_M1_H0[1:length(par_est_M1_H0)] - 2*mySE_M1_H0,
                           par_est_M1_H0[1:length(par_est_M1_H0)] + 2*mySE_M1_H0)
colnames(parMatrix_M1_H0) <- c("Param", "LowerB", "UpperB")
rownames(parMatrix_M1_H0) <- c(paste("A", 1:(K-1), sep=""),
                               paste("B", 1:(K-1), sep=""),
                               paste("D", 1:(K-1), sep=""),
                               rep("xi", (K-1)^2))
# round(parMatrix_M1_H0, digits=4)


#------- Fitted Value with VS. without Changepoint for each categories
len_mean_M1_H0 <- dim(DesignXEst_M1_H0)[2]
theta_M1_H0    <- matrix(par_est_M1_H0[1:((K-1)*len_mean_M1_H0)], nrow=K-1, ncol=len_mean_M1_H0)
eta_M1_H0      <- DesignXEst_M1_H0%*%t(theta_M1_H0)
gamm_M1_H0     <- exp(eta_M1_H0)/(1+exp(eta_M1_H0))
pm_M1_H0       <- matrix(0, Ti, K)
for (i in 1:(K-1)){
  if(i==1){
    pm_M1_H0[,i] <- gamm_M1_H0[,i]
  }else{
    pm_M1_H0[,i] <- gamm_M1_H0[,i] - gamm_M1_H0[,i-1]
  }
}
pm_M1_H0[,K] <- 1 - gamm_M1_H0[,K-1]


logL_M1_H0 <- LoglikCalcCpp(resEst_M1_H0$parEst, DesignXEst_M1_H0, y_log, mm, conv = mytol, stepNR=stepsize[2])
# data from server
load(paste0("application/MTMApplication_CandaCloudK5_SEA_30Years65to94_Hour", TargetHour, ".RData"))
LambdaStat_M1 <- -2*(logL_M1_H0 - CPCloudHaClc$logLHa)
tauest_M1_H0  <- which.max(LambdaStat_M1)


#################################################################
#               Parameter Estimation under Ha                   #   Model 1 with changepoint
#################################################################

tauEst_M1          <- tauClc[tauest_M1_H0]
CpValue            <- c(rep(0, tauEst_M1-1), rep(1, Ti-tauEst_M1+1))
DesignXEst_M1_Ha   <- cbind(DesignXEst_M1_H0, CpValue)
parinitial_M1_Ha   <- ParInitCalcu(y_log=y_log, m=mm, DesignX=DesignXEst_M1_Ha)
resEst_M1_Ha       <- SeqQuasiNWCpp(parinitial_M1_Ha, y_log, mm, DesignXEst_M1_Ha, stepsize, 1000, mytol)
par_est_M1_Ha      <- resEst_M1_Ha$parEst
mySEtheta_M1_Ha    <- sqrt(diag(solve(resEst_M1_Ha$Hesstheta)))
mySExi_M1_Ha       <- sqrt(diag(solve(resEst_M1_Ha$Hessxi)))
mySE_M1_Ha         <- c(mySEtheta_M1_Ha, mySExi_M1_Ha)
parMatrix_M1_Ha   <- cbind(par_est_M1_Ha[1:length(par_est_M1_Ha)],
                           par_est_M1_Ha[1:length(par_est_M1_Ha)] - 2*mySE_M1_Ha,
                           par_est_M1_Ha[1:length(par_est_M1_Ha)] + 2*mySE_M1_Ha)
colnames(parMatrix_M1_Ha) <- c("Param", "LowerB", "UpperB")
rownames(parMatrix_M1_Ha) <- c(paste("A", 1:(K-1), sep=""),
                               paste("B", 1:(K-1), sep=""),
                               paste("D", 1:(K-1), sep=""),
                               paste("Delta", 1:(K-1), sep=""),
                               rep("xi", (K-1)^2))
# round(parMatrix_M1_Ha, digits=4)

logL_M1_Ha <- LoglikCalcCpp(par_est_M1_Ha, DesignXEst_M1_Ha, y_log, mm, conv = mytol, stepNR=stepsize[2])
LambdaStat.Est_M1_Ha <- -2*(logL_M1_H0 - logL_M1_Ha)
myresid_M1_Ha = MTMResidCpp1(parEst=par_est_M1_Ha, y_log=y_log, mm=mm, DesignX=DesignXEst_M1_Ha, K=K, Ti=Ti, mytol=1e-5, 1)

#------- Fitted Value with VS. without Changepoint for each categories
len_mean_M1_Ha <- dim(DesignXEst_M1_Ha)[2]
theta_M1_Ha    <- matrix(par_est_M1_Ha[1:((K-1)*len_mean_M1_Ha)], nrow=K-1, ncol=len_mean_M1_Ha)
eta_M1_Ha      <- DesignXEst_M1_Ha%*%t(theta_M1_Ha)
gamm_M1_Ha     <- exp(eta_M1_Ha)/(1+exp(eta_M1_Ha))
pm_M1_Ha       <- matrix(0, Ti, K)
for (i in 1:(K-1)){
  if(i==1){
    pm_M1_Ha[,i] <- gamm_M1_Ha[,i]
  }else{
    pm_M1_Ha[,i] <- gamm_M1_Ha[,i] - gamm_M1_Ha[,i-1]
  }
}
pm_M1_Ha[,K] <- 1 - gamm_M1_Ha[,K-1]

save(dat2, 
     DesignXEst_M1_H0, par_est_M1_H0, parMatrix_M1_H0, pm_M1_H0, gamm_M1_H0, 
     tauClc, LambdaStat_M1, 
     DesignXEst_M1_Ha, par_est_M1_Ha, parMatrix_M1_Ha, pm_M1_Ha, gamm_M1_Ha,
     tauEst_M1, LambdaStat.Est_M1_Ha, myresid_M1_Ha,
     file = paste0("Application/ApplicationNeed_Hour", TargetHour, "_Model1.RData"))





par.tab = data.frame(Model2H0=rep("-", 36), 
                     Model2Ha=rep("-", 36),
                     Model1H0=rep("-", 36),
                     Model1Ha=rep("-", 36))

##### Model2H0
for(i in 1:16){
  par.tab[i,1] = paste0(round(parMatrix_M2_H0[i,1],4), 
                        " (", round(parMatrix_M2_H0[i,2],4), 
                        ", ", round(parMatrix_M2_H0[i,3],4),
                        ")")
}
for(i in 17:NROW(parMatrix_M2_H0)){
  par.tab[i+4,1] = paste0(round(parMatrix_M2_H0[i,1],4), 
                          " (", round(parMatrix_M2_H0[i,2],4), 
                          ", ", round(parMatrix_M2_H0[i,3],4),
                          ")")
}


##### Model2Ha
for(i in 1:NROW(parMatrix_M2_Ha)){
  par.tab[i,2] = paste0(round(parMatrix_M2_Ha[i,1],4), 
                        " (", round(parMatrix_M2_Ha[i,2],4), 
                        ", ", round(parMatrix_M2_Ha[i,3],4),
                        ")")
}


##### Model1H0

for(i in 1:4){
  par.tab[i,3] = paste0(round(parMatrix_M1_H0[i,1],4), 
                        " (", round(parMatrix_M1_H0[i,2],4), 
                        ", ", round(parMatrix_M1_H0[i,3],4),
                        ")")
}
for(i in 5:12){
  par.tab[i+4,3] = paste0(round(parMatrix_M1_H0[i,1],4), 
                          " (", round(parMatrix_M1_H0[i,2],4), 
                          ", ", round(parMatrix_M1_H0[i,3],4),
                          ")")
}
for(i in 13:NROW(parMatrix_M1_H0)){
  par.tab[i+8,3] = paste0(round(parMatrix_M1_H0[i,1],4), 
                          " (", round(parMatrix_M1_H0[i,2],4), 
                          ", ", round(parMatrix_M1_H0[i,3],4),
                          ")")
}


##### Model1Ha

for(i in 1:4){
  par.tab[i,4] = paste0(round(parMatrix_M1_Ha[i,1],4), 
                        " (", round(parMatrix_M1_Ha[i,2],4), 
                        ", ", round(parMatrix_M1_Ha[i,3],4),
                        ")")
}
for(i in 5:NROW(parMatrix_M1_Ha)){
  par.tab[i+4,4] = paste0(round(parMatrix_M1_Ha[i,1],4), 
                          " (", round(parMatrix_M1_Ha[i,2],4), 
                          ", ", round(parMatrix_M1_Ha[i,3],4),
                          ")")
}


library(xtable)

print(xtable(par.tab[,c(1,2,4)]), include.rownames=FALSE)
print(xtable(par.tab), include.rownames=FALSE)
