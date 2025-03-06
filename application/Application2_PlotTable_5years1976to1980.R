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

myTSK11_K5_mean = function(dat2, TargetHour, pm1, pm2, tauEst, mytitles){
  
  p = which(dat2$V1>=1976 & dat2$V1 <=1979)
  n = length(p)
  K = max(dat2$V6) + 1
  
  y = dat2[p,]
  pm1 = pm1[p,]
  pm2 = pm2[p,]
  
  ## plot parameters
  mylwd1 = 3
  mylwd2 = 3.5
  
  ## K = 11 Data
  # plot(c(1,n), y=c(0, 10), type='n',
  #      xlab="", ylab="Categories",
  #      xaxt="n", axes=FALSE, cex.lab=2.5)
  # box()
  # lines(y[,5], col="grey50", lwd=mylwd1)
  px1 <- rep(0, 12)
  for(i in 1:12){
    px1[i] <- which(as.numeric(y$V1) == 1978 & as.numeric(y$V2) == i)[1]
  }
  px2 <- rep(0, 12)
  for(i in 1:12){
    px2[i] <- which(as.numeric(y$V1) == 1979 & as.numeric(y$V2) == i)[1]
  }
  px <- c(px1, px2)
  # axis(1, at=px, labels=(format(as.Date(y$Date), "%Y-%m"))[px], padj=0.5, las=1.5, cex.lab=2.5)
  # axis(2, at=0:10, labels=1:11, las=1.5, cex.lab=2.5)
  # # title(paste0(mytitles, " (K=11)"), line = 1, cex.main = 2.5)
  # title(xlab="Time", line = 4, cex.main = 2.5)
  
  ## K = 5 Fitted Mean
  EXL1 = EXL2 = rep(0, n)
  for(t in 1:n){
    EXL1[t] = sum(pm1[t,]*(0:(K-1)))
    EXL2[t] = sum(pm2[t,]*(0:(K-1)))
  }
  
  ## K = 5 Data
  tmpp <- which(y$V1==dat2$V1[tauEst] & y$V2==dat2$V2[tauEst] & y$V3==dat2$V3[tauEst])
  plot(c(1,n), y=c(-0.5, 4), type='n', xlab="", ylab="Categories",
       xaxt="n", axes=FALSE, cex.lab=2.5)
  lines(y$V6, col="grey50", lwd=mylwd1)
  lines(1:(tmpp-1), EXL2[1:(tmpp-1)], col="darkorange", lty="solid", lwd=mylwd2)
  lines(tmpp:n, EXL2[tmpp:n], col="blue", lty="solid", lwd=mylwd2)
  lines(1:n, EXL1, col="red", lty="longdash", lwd=mylwd2)
  box()
  axis(1, at=px, labels=(format(as.Date(y$Date), "%Y-%m"))[px], las=1.5, padj=0.5, cex.lab=2.5)
  axis(2, at=0:4, labels=1:5, las=1.5)
  legend("bottom", horiz = T, legend=c("Without Changepoint", "Before Changepoint", "After Changepoint"),
         col=c("red", "darkorange", "blue"), lty=c("longdash", "solid", "solid"),
         lwd=rep(mylwd2,3), cex=2.5, bty="n", seg.len=1.2, x.intersp=0.4, text.width=210)
  # title(paste0(mytitles, " (K=5)"), line = 1, cex.main = 2.5)
  title(xlab="Time", line = 4, cex.main = 2.5)
}


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

TargetHour = 9                       ### change to 15 to obtain appendix results
mytitles = paste0("Hour ", TargetHour)

################################################################################
#######################         Figure 4               #########################
################################################################################
## 1978 and 1979 (Two Years Plot with fitted mean structure)

setEPS()
postscript(paste0("application/CloudDat1979mean_Hour", TargetHour, ".eps"), width =  16, height = 8)


# m <- matrix(c(1,2,3,4,5,5), nrow = 3, ncol = 2, byrow = TRUE)
# layout(mat = m,heights = c(0.48,0.48,0.04))

par(mfrow=c(1,1), cex.axis=2.5, cex.lab=2.5, mar=c(5,5,4,4))

load(paste0("application/ApplicationNeed_Hour", TargetHour, "_Model1_5years1976to1980.RData"))
myTSK11_K5_mean(dat2, TargetHour, pm_M1_H0, pm_M1_Ha, tauEst_M1, mytitles=mytitles)

# load(paste0("application/ApplicationNeed_Hour", TargetHour, "_Model2.RData"))
# myTSK11_K5_mean(dat2, TargetHour, pm_M2_H0, pm_M2_Ha, tauEst_M2, mytitles=mytitles)

# par(cex.axis=2.5, cex.lab=2.5, mar=c(1,1,1,1))
# plot(1, type = "n", axes=FALSE, xlab="", ylab="")
# legend(x = "center",inset = 0, horiz = T, legend=c("Without Changepoint", "Before Changepoint", "After Changepoint"),
#        col=c("red", "darkorange", "blue"), lty=c("longdash", "solid", "solid"), 
#        lwd=rep(3.5,3), cex=2.5, bty="n", seg.len=3)

dev.off()


################################################################################
#######################         Figure 6               #########################
################################################################################

load(paste0("application/ApplicationNeed_Hour", TargetHour, "_Model1_5years1976to1980.RData"))
Ti = dim(dat2)[1]
K = 5

setEPS()
postscript(paste0("application/fittedprob_Hour", TargetHour, ".eps"), width = 18, height = 8)

par(mfrow=c(2,3), mar=c(4,5,4,4))

tmpyear <- 1976:1980
y <- dat2[dat2$V1>=1976 & dat2$V1 <=1980,]
px <- rep(0, length(tmpyear))
for(i in 1:length(tmpyear)){
  px[i] <- which(y$V1 == tmpyear[i])[1]
}

for(k in 1:K){
  plot(c(0, Ti), c(0, 0.6), type="n", xlab="Time", ylab="Probability",
       axes=FALSE, main = paste("Category", k), cex.main=2.5, cex.lab=2.5, las=1)
  axis(1, at=px[seq(0, 30, 5)], labels=tmpyear[seq(0, 30, 5)], cex.axis=2.5)
  axis(2, at=seq(0, 0.6, 0.1), labels=seq(0, 0.6, 0.1), cex.axis=2.5)
  axis(2, at=seq(0, 0.6, 0.1), labels=F, cex.axis=2.5)
  box()
  # lines(1:Ti, pm_M2_H0[,k], col="black", lwd = 2)
  lines(1:Ti, pm_M1_Ha[,k], col="green4", lwd = 2)
  
  if(k==3){
    legend("topright", horiz = F, legend=c("Model 2 without changepoint", 
                                           "Model 1 with changepoint"),
           lty=c(1,1), col=c("black", "green4"), bty="n", cex=2.5, lwd=c(2,2))
  }
}

# plot(c(0, Ti), c(0, 0.6), type="n", xlab="", ylab="", axes=F)


Ti = NROW(dat2)
plot(x=c(0, Ti), y=c(0, 90), type="n", col="black", axes=FALSE,
     xlab="Time", ylab=expression(Lambda(tau)), xaxt="n", xaxs="i", 
     cex.main=2.5, cex.lab=2.5, las=1)
px <- rep(0, length(tmp))
for(i in 1:length(tmp)){
  px[i] <- which(dat2$V1 == tmp[i])[1]
}
axis(1, at=px, labels=1976:1980, cex.axis=2.5)
axis(2, at=seq(0, 90, 20), labels=seq(0, 90, 20), cex.axis=2.5)
box()
lines(x=tauClc, y=LambdaStat_M1)
abline(h=17.911, col="red", lty="dashed", lwd=2)
text(x=1800, y=17.911+9, labels="95th Percentile", cex = 2)
title("Log-likelihood ratio statistics", line = 1, cex.main = 2.5)
# title(xlab = "Time", line = 4, cex.lab = 2.5)


dev.off()


################################################################################
#######################         Figure S6              #########################
################################################################################

setEPS()
postscript(paste0("application/Lambda_max_Hour", TargetHour, "_sup.eps"), width =  16, height = 8)

par(mfrow=c(1,1), cex.axis=2.5, cex.lab=2.5, mar=c(6,6,4,4))

Ti = NROW(dat2)
plot(x=c(1, Ti), y=c(0, 230), type="n", col="black",
     xlab="", ylab=expression(Lambda(tau)), xaxt="n", xaxs="i", cex.lab=2.5)
lines(x=tauClc, y=LambdaStat_M1)
tmp <- 1976:1980
px <- rep(0, length(tmp))
for(i in 1:length(tmp)){
  px[i] <- which(dat2$V1 == tmp[i])[1]
}
axis(1, at = px, labels = tmp, las = 1, padj=0.5)
abline(h=17.911, col="red", lty="dashed", lwd=2)
text(x=1600, y=17.911+9, labels="95th Percentile", cex = 2)
title(xlab = "Time", line = 4, cex.main = 2.5)
title("Log-likelihood ratio statistics", line = 1, cex.main = 2.5)

dev.off()


################################################################################
#######################        Inference               #########################
################################################################################
## logit cumulative odds

K = 5
# average absolute change of DELTA in the logit of cumulative probabilities
mean(abs(par_est_M1_Ha[13:16]))
# [1] 0.462971

p.b <- p.a <- numeric(0)
for (i in 1:(K-1)){
  p.b[i] <- mean(gamm_M1_Ha[1:(tauEst_M1-1), i])
  p.a[i] <- mean(gamm_M1_Ha[tauEst_M1:Ti, i])
}
u <- numeric(0)
for (i in 1:(K-1)){
  u[i] <- log(p.a[i]/(1-p.a[i]))-log(p.b[i]/(1-p.b[i]))
}
mean(abs(u))
# [1] 0.4560516


################################################################################
#######################         Figure 7               #########################
################################################################################
## PIT Histograms


myPIT_hist = function(TargetHour, H){
  
  load(paste0("application/ApplicationNeed_Hour", TargetHour, "_Model1_5years1976to1980.RData"))
  
  Ti = dim(dat2)[1]
  K = 5
  mm = dat2$V6
  m = dat2$V5
  
  Pdist = t(apply(cbind(rep(0, Ti), myresid_M1_Ha$OneStepPred), 1, cumsum))
  
  Ti = NROW(Pdist)
  # H = 100
  
  Ftubar = rep(0, H)
  for(h in 1:H){
    u = h/H
    for(t in 1:Ti){
      pytminus1 = Pdist[t,mm[t]+1]
      pyt = Pdist[t,mm[t]+2]
      if(u < pytminus1) {
        Ftubar[h] = Ftubar[h] + 0
      }else if(u >= pyt){
        Ftubar[h] = Ftubar[h] + 1
      }else{
        Ftubar[h] = Ftubar[h] + (u-pytminus1)/(pyt-pytminus1)
      }
    }
    Ftubar[h] = Ftubar[h]/Ti
  }
  
  tmpvalue = c(0, Ftubar, 1)
  tmpvalue = tmpvalue[2:(H+1)] - tmpvalue[1:H]
  
  tmpx = barplot(tmpvalue, plot=T, col="lightblue", ylim=c(0,0.025), 
                 xlab="PIT histograms", ylab="Relative Frequency", yaxt = "n")
  axis(1, at=tmpx, label=seq(from=1/H, to=1, by=1/H), tick=F)
  axis(2, at=c(0, 0.005, 0.010, 0.015, 0.020, 0.025), label=c(0, 0.005, 0.010, 0.015, 0.020, 0.025), tick=T)
  abline(h=1/H, col="red", lty="dashed", lwd=2)
  # title(paste0("Hour ", TargetHour), line = 1, cex.main = 2.5)
  box()
  
}



setEPS()
postscript(paste0("application/CanadaCloudPIT_", TargetHour, ".eps"), width = 16, height = 8)

par(mfrow=c(1,1), cex.axis=2.5, cex.lab=2.5, mar=c(5,5,4,4))
myPIT_hist(TargetHour=TargetHour, H=50)

dev.off()


################################################################################
#######################         Figure 8               #########################
################################################################################
# one-step-ahead prediction residuals plots

setEPS()
postscript(paste0("application/Pc_Acf_", TargetHour, ".eps"), width = 16, height = 10)

par(mfcol=c(2,2), cex.axis=2.5, cex.lab=2.5, mar=c(6,7,5,5))
load(paste0("application/ApplicationNeed_Hour", TargetHour, "_Model1_5years1976to1980.RData"))
Ti = dim(dat2)[1]
K = 5
## Conditional Residuals ACF plots
for(iii in 1:(K-1)){
  numlag = 30
  tmpacf = stats::acf(myresid_M1_Ha$OneStepResid[,iii], lag=numlag, plot = F)
  plot(tmpacf$acf[2:numlag], 
       type="h", main=paste0("Category ", iii),
       xlab="Lag", ylim=c(-0.2, 0.2), ylab="",
       las=1, xaxt="n", lwd=3, cex.main=2.5)
  abline(h=0)
  title(ylab="ACF", line=4.5, cex.lab=2.5)
  # Add labels to the x-axis
  axis(1, at=c(1:numlag), labels=c(1:numlag))
  # Add 5% critical levels
  abline(h=c(2/sqrt(Ti),-2/sqrt(Ti)), lty=c(2,2), col="blue", lwd=3)
}

dev.off()