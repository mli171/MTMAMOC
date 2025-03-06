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
sourceCpp("src/myfuncClipPaper.cpp")

UniCpDetect <- function(res_sim = NULL, 
                        DesignXEst = NULL,
                        stepsize=NULL, 
                        BoundAdjust=NULL){
  
  #--------------- Input objects ---------------#
  X_hour <- res_sim$X_hour
  X_hour_wide <- res_sim$X_hour_wide
  Ts <- dim(X_hour_wide)[1]
  K <- dim(X_hour_wide)[2]
  
  #----------------- Estimation ----------------#
  cout <- summary(factor(X_hour))
  ci_initial <- as.vector(qnorm(cumsum(cout/sum(cout))))[1:(K-1)]
  p <- which(!is.finite(ci_initial))
  if(length(p)!=0){
    if(p==1){
      ci_initial[p] <- -6
    }else{
      ci_initial[p] <-  6
    }
  }
  if(sum(DesignXEst)==0){
    # stationary
    par_initial <- ci_initial
  }else{
    # non-stationary
    fri_initial <- rep(0, dim(DesignXEst)[2])
    par_initial <- c(ci_initial, fri_initial)
  }
  
  #------- Marginal Parameter Estimation
  res_NW <- NW_cpp(par=par_initial, X=X_hour_wide, DesignX=DesignXEst,
                   stepsize=stepsize, conv=1e-05)
  Marg_est <- res_NW$par
  MyHS <- res_NW$hs
  
  #------- Dependence Paramater Estimation
  phi_initial <- acf(X_hour, plot=F)$acf[2,1,1]
  MySampleCov <- sum(diag(t(X_hour_wide[1:(Ts-1),])%*%X_hour_wide[2:Ts,]))/Ts
  # Bound Adjustment
  myBd <- c(-0.99, 0.99)
  if(phi_initial < BoundAdjust){
    myBd <- c(-0.99, phi_initial)
  }
  OptimSecond2 <- optim(par=acf(X_hour, plot=F)$acf[2,1,1], fn=diffMOM,
                        method = "Brent", lower=-0.99, upper=0.99,
                        control=list(reltol=1e-05), Marg_est=Marg_est,
                        K=K, DesignX=DesignXEst, MySampleCov=MySampleCov)
  phi_est <- OptimSecond2$par
  par_est <- c(Marg_est, phi_est)
  
  #---------------- CUSUM on Y ----------------#
  seg_est <- c(-Inf, par_est[1:(K-1)], Inf)
  # Univariate Expectation
  if(length(par_est)==K){
    mst <- rep(0, Ts)
  }else{
    mst <- DesignXEst%*%par_est[K:(length(par_est)-1)]
  }
  EXL <- rep(0, Ts)
  for(t in 1:Ts){EXL[t] <- K - sum(pnorm(seg_est-mst[t])[2:K])}
  InnovRes <- UniInnovRcpp(EXL=EXL, mst=mst, Marg_est = Marg_est,
                           phi_est = phi_est, K=K, numCoef=1500)
  HQ <- InnovRes$HQ
  V <- InnovRes$V
  mylag <- InnovRes$lag
  rm(InnovRes)
  # One Step Ahead Prediction
  MyPred <- ClipPred(EXL=EXL, X_hour=X_hour, mylag=mylag, HQ=HQ)
  ItY    <- MyPred$Innovation/sqrt(V)
  # CUSUM Test Statistics
  MyCp   <- MyCusum(ErrErr=MyPred$Innovation, V=V)
  CpCsmY <- MyCp$Csm
  CpLocY <- MyCp$Location
  CpValY <- MyCp$Maximum
  
  #---------------- CUSUM on Z ----------------#
  # Estimation situation with tau1hat modified
  # Latent variable reconstruction
  ResLatent <- Re_latent(par = par_est, DesignX = DesignXEst, X = X_hour_wide)
  Z_expct <- ResLatent$Z_expct
  mst <- ResLatent$mst
  epsilonHat1 <- c(0,(Z_expct[2:Ts] - mst[2:Ts]) - phi_est*(Z_expct[1:(Ts-1)] - mst[1:(Ts-1)]))
  ResCusum1 <- rep(0, Ts)
  for(tau in 1:Ts){
    ResCusum1[tau] <- (sum(epsilonHat1[1:tau])-tau/Ts*sum(epsilonHat1[1:Ts]))/sqrt(Ts)
  }
  # tau1 variance approximation
  qn <- floor((Ts)^(1/3))
  epsilonbar <- mean(epsilonHat1)
  tmpsum <- 0
  for(s in 1:qn){
    tmpsum1 <- 0
    for(t in 1:(Ts-s)){
      tmpsum1 <- tmpsum1 + (epsilonHat1[t]-epsilonbar)*(epsilonHat1[t+s]-epsilonbar)
    }
    # tmpsum <- tmpsum + (1-s/(qn+1))*tmpsum1/(Ts-s)
    tmpsum <- tmpsum + (1-s/(qn+1))*tmpsum1/Ts
  }
  tau1hat <- sum((epsilonHat1-epsilonbar)^2)/Ts + 2*tmpsum
  sigmahat1 <- sqrt(tau1hat)
  ItZ <- abs(ResCusum1/sigmahat1)
  CpLocZ <- which.max(ItZ)
  CpValZ <- max(ItZ)
  
  res <- list(CpLocY, CpValY, ItY, CpCsmY, CpLocZ, CpValZ, ItZ, par_est, res_NW$hs)
  names(res) <- c("CpLocY", "CpValY", "ItY", "CpCsmY", "CpLocZ", "CpValZ", "ItZ", "ParEst", "Hessian")
  return(res)
}

# Expectation of Truncated Normal
expct_TN <- function(mu=NULL, sigma=NULL, cl=NULL, cu=NULL){
  return(mu+sqrt(sigma)*(dnorm(x=cl,mean=mu,sd=sqrt(sigma))-dnorm(x=cu,mean=mu,sd=sqrt(sigma)))/
           (pnorm(q=cu,mean=mu,sd=sqrt(sigma))-pnorm(q=cl,mean=mu,sd=sqrt(sigma))))
}

# Latent Variable Reconstruction
Re_latent <- function(par=NULL, DesignX=NULL, X=NULL){
  
  Ts <- dim(X)[1]
  K <- dim(X)[2]
  
  # Parameters Extraction
  ci <-  par[1:(K-1)]
  seg <- c(-Inf, ci, Inf)
  len_par <- length(par)
  phi <- par[len_par]
  
  if(length(par)==K){
    mst <- rep(0, Ts)
  }else{
    mst <- DesignX%*%par[K:(len_par-1)]
  }
  
  sig2 <- 1 # latent variable variance fix to be 1
  # allocate
  Z_expct <- rep(0, Ts)
  for(h in 1:Ts){
    pp <- which(X[h,]==1)
    Z_expct[h] <- expct_TN(mu=mst[h], sigma=sig2, cl=seg[pp], cu=seg[pp+1])
  }
  
  res <- list(mst, Z_expct)
  names(res) <- c("mst", "Z_expct")
  return(res)
}


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

datdir = "~/OneDrive - University of Louisiana Lafayette/data/"
dat = read.table(file=paste0(datdir, "CloudCover/FortStJohnAirport/hourly_day.txt"))
dat$Date = ISOdate(year = dat$V1, month = dat$V2, day = dat$V3, hour = dat$V4)

## filter1: 30 years data from 1975 to 1980
dat1 = dat[dat$V1 >= 1976 & dat$V1 <=1980,]
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

nCoreuse <- 10
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

##### Design Matrix
res_sim = list(X_hour=mm, X_hour_wide=y_log)
DesignXH0  = cbind(AValue, BValue, DValue)
FitRes11 = try(UniCpDetect(res_sim = res_sim, DesignXEst=DesignXH0, stepsize = 1, BoundAdjust = -0.5))


# setEPS()
# postscript("ApplyFigure3.eps", width = 14, height = 7)
par(cex.lab=1, cex.axis=1.2, mfrow=c(1,2))
plot(x=1:Ti, y=FitRes11$CpCsm, type="l", col="black",
     xlab="Year", ylab=expression('CUSUM'['I'](tau)),
     xaxt="n",
     xlim=c(-100, Ti+100), ylim=c(0, 2), xaxs="i")
abline(h=1.355, lty="dashed", col="black")
text(x=350, y=0.98, labels="95th Percentile")
tmp <- seq(from=min(dat1$V1), to=max(dat1$V1), by=5)
px <- rep(0, length(tmp))
for(i in 1:length(tmp)){
  px[i] <- which(dat2$V1 == tmp[i])[1]
}
axis(1, at = px, labels = tmp, las = 1)
# dev.off()
abline(v=FitRes11$CpLocY, col="red")



# setEPS()
# postscript("ApplyFigure4.eps", width = 14, height = 7)
# par(cex.lab=1, cex.axis=1.2, mfrow=c(1,1))
plot(x=1:Ti, y=FitRes11$ItZ, type="l", col="black",
     xlab="Day", 
     # ylab=expression('CUSUM'['e'](tau)),
     ylab="CUSUM Statistics",
     xaxt="n",
     xlim=c(-30, Ti+30), ylim=c(0, 2), xaxs="i")
abline(h=1.355, lty="dashed", col="black")
text(x=350, y=0.96, labels="95th Percentile")
tmp <- seq(from=min(dat1$V1), to=max(dat1$V1), by=5)
px <- rep(0, length(tmp))
for(i in 1:length(tmp)){
  px[i] <- which(dat2$V1 == tmp[i])[1]
}
axis(1, at = px, labels = tmp, las = 1)
# dev.off()
abline(v=FitRes11$CpLocZ, col="red")



dat2[FitRes11$CpLocY,]
dat2[FitRes11$CpLocZ,]

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
load(paste0("application/res/MTMApplication_CandaCloudK5_FULL_5Years76to80_Hour", TargetHour, ".RData"))
LambdaStat_M2 <- -2*(logL_M2_H0 - CPCloudHaClc$logLHa) # 19.599
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
     file = paste0("application/ApplicationNeed_Hour", TargetHour, "_Model2_5years1976to1980.RData"))
















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
load(paste0("application/res/MTMApplication_CandaCloudK5_SEA_5Years76to80_Hour", TargetHour, ".RData"))
LambdaStat_M1 <- -2*(logL_M1_H0 - CPCloudHaClc$logLHa) # 17.757
tauest_M1_H0  <- which.max(LambdaStat_M1) # [1] 939

max(LambdaStat_M1) # [1] 35.49049
dat2[tauest_M1_H0,]
#         V1 V2 V3 V4 V5                Date V6
# 42577 1978  7 28  9  9 1978-07-28 09:00:00  3

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
     file = paste0("application/ApplicationNeed_Hour", TargetHour, "_Model1_5years1976to1980.RData"))





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
