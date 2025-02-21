library(xtable)

load("2ParameterEstimation/MTM_TRDSEA_K5_T3650_SS365_ParamEst_PositiveCorr_H0_1th.RData")


# computation time (units in "minutes")
round(mean(ParamFinal$tim), digits = 2)
# [1] 0.04
ParamMatrix <- ParamFinal[, 7:38]
# True
TruePos <- c(-1.3863, -0.4055, 0.4055, 1.3863, 0.1, 0.1, 0.1, 0.1, -0.1, -0.2, -0.15, -0.3, 0.2, 0.1, 0.15, 0.3, 
             2.8,  2.2,  1.9,  1.0,
             1.3,  1.2,  0.9,  0.6,
             0.8,  0.8,  0.5,  0.7,
             0.5,  0.4,  0.3,  0.3)

# Mean
MeanPos <- colMeans(ParamMatrix)
# Bias
BiasPos <- (TruePos - MeanPos)/TruePos
# Standard Deviation
SdPos   <- apply(ParamMatrix, 2, sd)
# Latex
xtable(data.frame(TruePos, MeanPos, BiasPos, SdPos), digits=4)








tmp <- apply(ParamMatrix, 2, range)
tmp



setEPS()
postscript("Paper/fig/PosParam_Mean.eps", width = 14, height = 10)
mycexfont <- 2.5
myupperY  <- 5000
# mylabelX  <- 0.15
mylabelY  <- 4500
Meantmp <- MeanPos
par(cex.lab=1.5, cex.axis=1.5, mfrow=c(4,4))

# (A1)
mylabelX <- 0.35
myxlim <- c(Meantmp[1] - 0.7, Meantmp[1] + 0.7)
mybreaks <- seq(from=floor(tmp[1,1]*100)/100,
                to=ceiling(tmp[2,1]*100)/100,
                by=0.01)
hist(ParamMatrix$'1', 
     xlim=myxlim,
     ylim=c(0,myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(-1.3863+mylabelX-0.2, mylabelY, labels = expression(alpha[1]==-1.3863), cex=mycexfont)
abline(v=Meantmp[1], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (A2)
mylabelX <- 0.35
myxlim <- c(Meantmp[2] - 0.7, Meantmp[2] + 0.7)
mybreaks <- seq(from=floor(tmp[1,2]*100)/100,
                to=ceiling(tmp[2,2]*100)/100,
                by=0.01)
hist(ParamMatrix$'2', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(-0.4055+mylabelX-0.2, mylabelY, labels = expression(alpha[2]==-0.4055), cex=mycexfont)
abline(v=Meantmp[2], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (A3)
mylabelX <- 0.35
myxlim <- c(Meantmp[3] - 0.7, Meantmp[3] + 0.7)
mybreaks <- seq(from=floor(tmp[1,3]*100)/100,
                to=ceiling(tmp[2,3]*100)/100,
                by=0.01)
hist(ParamMatrix$'3', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(0.4055+mylabelX-0.2, mylabelY, labels = expression(alpha[3]==0.4055), cex=mycexfont)
abline(v=Meantmp[3], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (A4)
mylabelX <- 0.35
myxlim <- c(Meantmp[4] - 0.7, Meantmp[4] + 0.7)
mybreaks <- seq(from=floor(tmp[1,4]*100)/100,
                to=ceiling(tmp[2,4]*100)/100,
                by=0.01)
hist(ParamMatrix$'4', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(1.3863+mylabelX-0.2, mylabelY, labels = expression(alpha[4]==1.3863), cex=mycexfont)
abline(v=Meantmp[4], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (B1)
mylabelX <- 0.35
myxlim <- c(Meantmp[9] - 0.7, Meantmp[9] + 0.7)
mybreaks <- seq(from=floor(tmp[1,9]*100)/100,
                to=ceiling(tmp[2,9]*100)/100,
                by=0.01)
hist(ParamMatrix$'9', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(-0.1+mylabelX, mylabelY, labels = expression('B'[1]==-0.1), cex=mycexfont)
abline(v=Meantmp[9], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (B2)
mylabelX <- 0.35
myxlim <- c(Meantmp[10] - 0.7, Meantmp[10] + 0.7)
mybreaks <- seq(from=floor(tmp[1,10]*100)/100,
                to=ceiling(tmp[2,10]*100)/100,
                by=0.01)
hist(ParamMatrix$'10', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(-0.2+mylabelX, mylabelY, labels = expression('B'[2]==-0.2), cex=mycexfont)
abline(v=Meantmp[10], col="red", lwd=1.5)
box(which = "plot", lty = "solid")

# (B3)
mylabelX <- 0.35
myxlim <- c(Meantmp[11] - 0.7, Meantmp[11] + 0.7)
mybreaks <- seq(from=floor(tmp[1,11]*100)/100,
                to=ceiling(tmp[2,11]*100)/100,
                by=0.01)
hist(ParamMatrix$'11', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(-0.15+mylabelX-0.1, mylabelY, labels = expression('B'[3]==-0.15), cex=mycexfont)
abline(v=Meantmp[11], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (B4)
mylabelX <- 0.35
myxlim <- c(Meantmp[12] - 0.7, Meantmp[12] + 0.7)
mybreaks <- seq(from=floor(tmp[1,12]*100)/100,
                to=ceiling(tmp[2,12]*100)/100,
                by=0.01)
hist(ParamMatrix$'12', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(-0.3+mylabelX, mylabelY, labels = expression('B'[4]==-0.3), cex=mycexfont)
abline(v=Meantmp[12], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (D1)
mylabelX <- 0.35
myxlim <- c(Meantmp[13] - 0.7, Meantmp[13] + 0.7)
mybreaks <- seq(from=floor(tmp[1,13]*100)/100,
                to=ceiling(tmp[2,13]*100)/100,
                by=0.01)
hist(ParamMatrix$'13', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(0.2+mylabelX, mylabelY, labels = expression('D'[1]==0.2), cex=mycexfont)
abline(v=Meantmp[13], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (D2)
mylabelX <- 0.35
myxlim <- c(Meantmp[14] - 0.7, Meantmp[14] + 0.7)
mybreaks <- seq(from=floor(tmp[1,14]*100)/100,
                to=ceiling(tmp[2,14]*100)/100,
                by=0.01)
hist(ParamMatrix$'14', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(0.1+mylabelX, mylabelY, labels = expression('D'[2]==0.1), cex=mycexfont)
abline(v=Meantmp[14], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (D3)
mylabelX <- 0.35
myxlim <- c(Meantmp[15] - 0.7, Meantmp[15] + 0.7)
mybreaks <- seq(from=floor(tmp[1,15]*100)/100,
                to=ceiling(tmp[2,15]*100)/100,
                by=0.01)
hist(ParamMatrix$'15', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(0.15+mylabelX, mylabelY, labels = expression('D'[3]==0.15), cex=mycexfont)
abline(v=Meantmp[15], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (D4)
mylabelX <- 0.35
myxlim <- c(Meantmp[16] - 0.7, Meantmp[16] + 0.7)
mybreaks <- seq(from=floor(tmp[1,16]*100)/100,
                to=ceiling(tmp[2,16]*100)/100,
                by=0.01)
hist(ParamMatrix$'16', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(0.3+mylabelX, mylabelY, labels = expression('D'[4]==0.3), cex=mycexfont)
abline(v=Meantmp[16], col="red", lwd=1.5)
box(which = "plot", lty = "solid")

# (beta1)
mylabelX <- 0.25
myxlim <- c(Meantmp[5] - 0.7, Meantmp[5] + 0.7)
mybreaks <- seq(from=floor(tmp[1,5]*100)/100,
                to=ceiling(tmp[2,5]*100)/100,
                by=0.01)
hist(ParamMatrix$'5', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(0.1+mylabelX, mylabelY, labels = expression(beta[1]==0.1), cex=mycexfont)
abline(v=Meantmp[5], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (beta2)
mylabelX <- 0.25
myxlim <- c(Meantmp[6] - 0.7, Meantmp[6] + 0.7)
mybreaks <- seq(from=floor(tmp[1,6]*100)/100,
                to=ceiling(tmp[2,6]*100)/100,
                by=0.01)
hist(ParamMatrix$'6', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(0.1+mylabelX, mylabelY, labels = expression(beta[2]==0.1), cex=mycexfont)
abline(v=Meantmp[6], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (beta3)
mylabelX <- 0.25
myxlim <- c(Meantmp[7] - 0.7, Meantmp[7] + 0.7)
mybreaks <- seq(from=floor(tmp[1,7]*100)/100,
                to=ceiling(tmp[2,7]*100)/100,
                by=0.01)
hist(ParamMatrix$'7', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(0.1+mylabelX, mylabelY, labels = expression(beta[3]==0.1), cex=mycexfont)
abline(v=Meantmp[7], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (beta4)
mylabelX <- 0.25
myxlim <- c(Meantmp[8] - 0.7, Meantmp[8] + 0.7)
mybreaks <- seq(from=floor(tmp[1,8]*100)/100,
                to=ceiling(tmp[2,8]*100)/100,
                by=0.01)
hist(ParamMatrix$'8', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(0.1+mylabelX, mylabelY, labels = expression(beta[4]==0.1), cex=mycexfont)
abline(v=Meantmp[8], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


dev.off()







setEPS()
postscript("Paper/fig/PosParam_Dependence.eps", width = 14, height = 10)

par(cex.lab=1.5, cex.axis=1.5, mfrow=c(4,4))
# (xi11)
mylabelX <- 0.5
myxlim <- c(Meantmp[17] - 1, Meantmp[17] + 1)
mybreaks <- seq(from=floor(tmp[1,17]*100)/100,
                to=ceiling(tmp[2,17]*100)/100,
                by=0.01)
hist(ParamMatrix$'17', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(2.8+mylabelX, mylabelY, labels = expression(xi[11]==2.8), cex=mycexfont)
abline(v=Meantmp[17], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (xi12)
mylabelX <- 0.5
myxlim <- c(Meantmp[18] - 1, Meantmp[18] + 1)
mybreaks <- seq(from=floor(tmp[1,18]*100)/100,
                to=ceiling(tmp[2,18]*100)/100,
                by=0.01)
hist(ParamMatrix$'18', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(2.2+mylabelX, mylabelY, labels = expression(xi[12]==2.2), cex=mycexfont)
abline(v=Meantmp[18], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (xi13)
mylabelX <- 0.5
myxlim <- c(Meantmp[19] - 1, Meantmp[19] + 1)
mybreaks <- seq(from=floor(tmp[1,19]*100)/100,
                to=ceiling(tmp[2,19]*100)/100,
                by=0.01)
hist(ParamMatrix$'19', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(1.9+mylabelX, mylabelY, labels = expression(xi[13]==1.9), cex=mycexfont)
abline(v=Meantmp[19], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (xi14)
mylabelX <- 0.5
myxlim <- c(Meantmp[20] - 1, Meantmp[20] + 1)
mybreaks <- seq(from=floor(tmp[1,20]*100)/100,
                to=ceiling(tmp[2,20]*100)/100,
                by=0.01)
hist(ParamMatrix$'20', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(1.0+mylabelX, mylabelY, labels = expression(xi[14]==1.0), cex=mycexfont)
abline(v=Meantmp[20], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (xi21)
mylabelX <- 0.5
myxlim <- c(Meantmp[21] - 1, Meantmp[21] + 1)
mybreaks <- seq(from=floor(tmp[1,21]*100)/100,
                to=ceiling(tmp[2,21]*100)/100,
                by=0.01)
hist(ParamMatrix$'21', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(1.3+mylabelX, mylabelY, labels = expression(xi[21]==1.3), cex=mycexfont)
abline(v=Meantmp[21], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (xi22)
mylabelX <- 0.5
myxlim <- c(Meantmp[22] - 1, Meantmp[22] + 1)
mybreaks <- seq(from=floor(tmp[1,22]*100)/100,
                to=ceiling(tmp[2,22]*100)/100,
                by=0.01)
hist(ParamMatrix$'22', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(1.2+mylabelX, mylabelY, labels = expression(xi[22]==1.2), cex=mycexfont)
abline(v=Meantmp[22], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (xi23)
mylabelX <- 0.5
myxlim <- c(Meantmp[23] - 1, Meantmp[23] + 1)
mybreaks <- seq(from=floor(tmp[1,23]*100)/100,
                to=ceiling(tmp[2,23]*100)/100,
                by=0.01)
hist(ParamMatrix$'23', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(0.9+mylabelX, mylabelY, labels = expression(xi[23]==0.9), cex=mycexfont)
abline(v=Meantmp[23], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (xi24)
mylabelX <- 0.5
myxlim <- c(Meantmp[24] - 1, Meantmp[24] + 1)
mybreaks <- seq(from=floor(tmp[1,24]*100)/100,
                to=ceiling(tmp[2,24]*100)/100,
                by=0.01)
hist(ParamMatrix$'24', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(0.6+mylabelX, mylabelY, labels = expression(xi[24]==0.6), cex=mycexfont)
abline(v=Meantmp[24], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (xi31)
mylabelX <- 0.5
myxlim <- c(Meantmp[25] - 1, Meantmp[25] + 1)
mybreaks <- seq(from=floor(tmp[1,25]*100)/100,
                to=ceiling(tmp[2,25]*100)/100,
                by=0.01)
hist(ParamMatrix$'25', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(0.8+mylabelX, mylabelY, labels = expression(xi[31]==0.8), cex=mycexfont)
abline(v=Meantmp[25], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (xi32)
mylabelX <- 0.5
myxlim <- c(Meantmp[26] - 1, Meantmp[26] + 1)
mybreaks <- seq(from=floor(tmp[1,26]*100)/100,
                to=ceiling(tmp[2,26]*100)/100,
                by=0.01)
hist(ParamMatrix$'26', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(0.8+mylabelX, mylabelY, labels = expression(xi[32]==0.8), cex=mycexfont)
abline(v=Meantmp[26], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (xi33)
mylabelX <- 0.5
myxlim <- c(Meantmp[27] - 1, Meantmp[27] + 1)
mybreaks <- seq(from=floor(tmp[1,27]*100)/100,
                to=ceiling(tmp[2,27]*100)/100,
                by=0.01)
hist(ParamMatrix$'27', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(0.5+mylabelX, mylabelY, labels = expression(xi[33]==0.5), cex=mycexfont)
abline(v=Meantmp[27], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (xi34)
mylabelX <- 0.5
myxlim <- c(Meantmp[28] - 1, Meantmp[28] + 1)
mybreaks <- seq(from=floor(tmp[1,28]*100)/100,
                to=ceiling(tmp[2,28]*100)/100,
                by=0.01)
hist(ParamMatrix$'28', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(0.7+mylabelX, mylabelY, labels = expression(xi[34]==0.7), cex=mycexfont)
abline(v=Meantmp[28], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (xi41)
mylabelX <- 0.5
myxlim <- c(Meantmp[29] - 1, Meantmp[29] + 1)
mybreaks <- seq(from=floor(tmp[1,29]*100)/100,
                to=ceiling(tmp[2,29]*100)/100,
                by=0.01)
hist(ParamMatrix$'29', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(0.5+mylabelX, mylabelY, labels = expression(xi[41]==0.5), cex=mycexfont)
abline(v=Meantmp[29], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (xi42)
mylabelX <- 0.5
myxlim <- c(Meantmp[30] - 1, Meantmp[30] + 1)
mybreaks <- seq(from=floor(tmp[1,30]*100)/100,
                to=ceiling(tmp[2,30]*100)/100,
                by=0.01)
hist(ParamMatrix$'30', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(0.4+mylabelX, mylabelY, labels = expression(xi[42]==0.4), cex=mycexfont)
abline(v=Meantmp[30], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (xi43)
mylabelX <- 0.5
myxlim <- c(Meantmp[31] - 1, Meantmp[31] + 1)
mybreaks <- seq(from=floor(tmp[1,31]*100)/100,
                to=ceiling(tmp[2,31]*100)/100,
                by=0.01)
hist(ParamMatrix$'31', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(0.3+mylabelX, mylabelY, labels = expression(xi[43]==0.3), cex=mycexfont)
abline(v=Meantmp[31], col="red", lwd=1.5)
box(which = "plot", lty = "solid")


# (xi44)
mylabelX <- 0.5
myxlim <- c(Meantmp[32] - 1, Meantmp[32] + 1)
mybreaks <- seq(from=floor(tmp[1,32]*100)/100,
                to=ceiling(tmp[2,32]*100)/100,
                by=0.01)
hist(ParamMatrix$'32', 
     xlim=myxlim,
     ylim=c(0, myupperY),
     main=NULL, xlab=NULL, freq=T,
     breaks=mybreaks)
text(0.3+mylabelX, mylabelY, labels = expression(xi[44]==0.3), cex=mycexfont)
abline(v=Meantmp[32], col="red", lwd=1.5)
box(which = "plot", lty = "solid")

dev.off()










##### QQ plots



# setEPS()
# postscript("Paper/fig/PosParamQQplots1.eps", width = 16, height = 16)

# setEPS()
jpeg("Paper/fig/PosParamQQplots1.jpeg", width = 4800, height = 4800, res=300)

par(cex.lab=2, cex.axis=2, cex.main = 2.3, mfrow=c(4,4), mar = c(6, 6, 6, 1))

# (A1)
qqnorm(ParamMatrix$'1', pch = 1, main=expression(alpha[1]),  frame = TRUE)
qqline(ParamMatrix$'1', col = "red", lwd = 1.5)
# (A2)
qqnorm(ParamMatrix$'2', pch = 1, main=expression(alpha[2]),  frame = TRUE)
qqline(ParamMatrix$'2', col = "red", lwd = 1.5)
# (A3)
qqnorm(ParamMatrix$'3', pch = 1, main=expression(alpha[3]),  frame = TRUE)
qqline(ParamMatrix$'3', col = "red", lwd = 1.5)
# (A4)
qqnorm(ParamMatrix$'4', pch = 1, main=expression(alpha[4]),  frame = TRUE)
qqline(ParamMatrix$'4', col = "red", lwd = 1.5)
# (B1)
qqnorm(ParamMatrix$'9', pch = 1, main=expression('B'[1]),  frame = TRUE)
qqline(ParamMatrix$'9', col = "red", lwd = 1.5)
# (B2)
qqnorm(ParamMatrix$'10', pch = 1, main=expression('B'[2]),  frame = TRUE)
qqline(ParamMatrix$'10', col = "red", lwd = 1.5)
# (B3)
qqnorm(ParamMatrix$'11', pch = 1, main=expression('B'[3]),  frame = TRUE)
qqline(ParamMatrix$'11', col = "red", lwd = 1.5)
# (B4)
qqnorm(ParamMatrix$'12', pch = 1, main=expression('B'[4]),  frame = TRUE)
qqline(ParamMatrix$'12', col = "red", lwd = 1.5)
# (D1)
qqnorm(ParamMatrix$'13', pch = 1, main=expression('D'[1]),  frame = TRUE)
qqline(ParamMatrix$'13', col = "red", lwd = 1.5)
# (D2)
qqnorm(ParamMatrix$'14', pch = 1, main=expression('D'[2]),  frame = TRUE)
qqline(ParamMatrix$'14', col = "red", lwd = 1.5)
# (D3)
qqnorm(ParamMatrix$'15', pch = 1, main=expression('D'[3]),  frame = TRUE)
qqline(ParamMatrix$'15', col = "red", lwd = 1.5)
# (D4)
qqnorm(ParamMatrix$'16', pch = 1, main=expression('D'[4]),  frame = TRUE)
qqline(ParamMatrix$'16', col = "red", lwd = 1.5)
# (beta1)
qqnorm(ParamMatrix$'5', pch = 1, main=expression(beta[1]),  frame = TRUE)
qqline(ParamMatrix$'5', col = "red", lwd = 1.5)
# (beta2)
qqnorm(ParamMatrix$'6', pch = 1, main=expression(beta[2]),  frame = TRUE)
qqline(ParamMatrix$'6', col = "red", lwd = 1.5)
# (beta3)
qqnorm(ParamMatrix$'7', pch = 1, main=expression(beta[3]),  frame = TRUE)
qqline(ParamMatrix$'7', col = "red", lwd = 1.5)
# (beta4)
qqnorm(ParamMatrix$'8', pch = 1, main=expression(beta[4]),  frame = TRUE)
qqline(ParamMatrix$'8', col = "red", lwd = 1.5)


dev.off()

# setEPS()
# postscript("Paper/fig/PosParamQQplots2.eps",
#            width = 16, height = 16)

jpeg("Paper/fig/PosParamQQplots2.jpeg", width = 4800, height = 4800, res=300)
par(cex.lab=2, cex.axis=2, cex.main = 2.3, mfrow=c(4,4), mar = c(6, 6, 6, 1))


# (xi11)
qqnorm(ParamMatrix$'17', pch = 1, main=expression(xi[11]),  frame = TRUE)
qqline(ParamMatrix$'17', col = "red", lwd = 1.5)
# (xi12)
qqnorm(ParamMatrix$'18', pch = 1, main=expression(xi[12]),  frame = TRUE)
qqline(ParamMatrix$'18', col = "red", lwd = 1.5)
# (xi13)
qqnorm(ParamMatrix$'19', pch = 1, main=expression(xi[13]),  frame = TRUE)
qqline(ParamMatrix$'19', col = "red", lwd = 1.5)
# (xi14)
qqnorm(ParamMatrix$'20', pch = 1, main=expression(xi[14]),  frame = TRUE)
qqline(ParamMatrix$'20', col = "red", lwd = 1.5)
# (xi21)
qqnorm(ParamMatrix$'21', pch = 1, main=expression(xi[21]),  frame = TRUE)
qqline(ParamMatrix$'21', col = "red", lwd = 1.5)
# (xi22)
qqnorm(ParamMatrix$'22', pch = 1, main=expression(xi[22]),  frame = TRUE)
qqline(ParamMatrix$'22', col = "red", lwd = 1.5)
# (xi23)
qqnorm(ParamMatrix$'23', pch = 1, main=expression(xi[23]),  frame = TRUE)
qqline(ParamMatrix$'23', col = "red", lwd = 1.5)
# (xi24)
qqnorm(ParamMatrix$'24', pch = 1, main=expression(xi[24]),  frame = TRUE)
qqline(ParamMatrix$'24', col = "red", lwd = 1.5)
# (xi31)
qqnorm(ParamMatrix$'25', pch = 1, main=expression(xi[31]),  frame = TRUE)
qqline(ParamMatrix$'25', col = "red", lwd = 1.5)
# (xi32)
qqnorm(ParamMatrix$'26', pch = 1, main=expression(xi[32]),  frame = TRUE)
qqline(ParamMatrix$'26', col = "red", lwd = 1.5)
# (xi33)
qqnorm(ParamMatrix$'27', pch = 1, main=expression(xi[33]),  frame = TRUE)
qqline(ParamMatrix$'27', col = "red", lwd = 1.5)
# (xi34)
qqnorm(ParamMatrix$'28', pch = 1, main=expression(xi[34]),  frame = TRUE)
qqline(ParamMatrix$'28', col = "red", lwd = 1.5)
# (xi41)
qqnorm(ParamMatrix$'29', pch = 1, main=expression(xi[41]),  frame = TRUE)
qqline(ParamMatrix$'29', col = "red", lwd = 1.5)
# (xi42)
qqnorm(ParamMatrix$'30', pch = 1, main=expression(xi[42]),  frame = TRUE)
qqline(ParamMatrix$'30', col = "red", lwd = 1.5)
# (xi43)
qqnorm(ParamMatrix$'31', pch = 1, main=expression(xi[43]),  frame = TRUE)
qqline(ParamMatrix$'31', col = "red", lwd = 1.5)
# (xi44)
qqnorm(ParamMatrix$'32', pch = 1, main=expression(xi[44]),  frame = TRUE)
qqline(ParamMatrix$'32', col = "red", lwd = 1.5)


dev.off()