library(VGAM)
library(mvtnorm)
library(foreach)
library(doMC)
library(Rcpp)
library(RcppArmadillo)

sourceCpp("Tools/QuasiNewton.cpp")

#-------------------------- Raw Estimation of Dependence Model Parameters
initial_xi <- function(m=NULL, avg_delta=NULL){
  
  # Naive Transition Matrix Estimation
  ct <- matrix(0, K, K)  # previous*current
  for(i in 1:K){
    p <- which(m==(i-1)) # where is i-1 category
    for (j in 1:length(p)){
      ct[i,m[p[j]+1]+1] <- ct[i,m[p[j]+1]+1] + 1
    }
  }
  pc_avg <- ct/rowSums(ct)   # pc_avg: j*k
  
  if(any(pc_avg == 0)){
    cat("Zero Exist in Transition Matrix!!!\n")
    pc_avg[which(pc_avg==0)] <- 0.00001
    # break
  }
  xi <- matrix(0, K-1, K-1) # previous*current
  for(j in 1:(K-1)){
    for(k in 1:(K-1)){
      # calculate xi
      xi[j,k] <- log(pc_avg[j,k]/pc_avg[j,K])-avg_delta[k]
    }
  }
  return(t(xi)) # res: current*previous
}

#-------------------------- Initial Values for All Model Parameters
ParInitCalcu <- function(y_log=NULL, m=NULL, DesignX=NULL){
  
  require(VGAM)
  
  K <- dim(y_log)[2]
  len_mean <- dim(DesignX)[2]
  if(len_mean == 1){
    fit.d <- vglm(y_log ~ 1, cumulative(parallel=FALSE))
  }else{
    fit.d <- vglm(y_log ~ DesignX[,2:dim(DesignX)[2]], cumulative(parallel=FALSE))
  }
  avg_delta <- colMeans(fit.d@fitted.values)
  XiInit <- initial_xi(m=m, avg_delta=avg_delta)
  
  parinit <- c(as.vector(fit.d@coefficients), as.vector(XiInit))
  return(parinit)
}

#-------------------------- Genetic Algorithm fit function
GAfit <- function(tau=NULL, DesignXH0=NULL, y_log=NULL, mm=NULL, 
                  stepsize=NULL, mytol=NULL){
  
  Ti = dim(DesignXH0)[1]
  
  CpValue    <- c(rep(0, tau-1), rep(1, Ti-tau+1))
  DesignXHa  <- cbind(DesignXH0, CpValue)
  parinitialHa  <- ParInitCalcu(y_log=y_log, m=mm, DesignX=DesignXHa)
  resEstHa <- SeqNewtonCpp(parinitialHa, y_log, mm, DesignXHa, 
                           stepsize, 1000, mytol)
  parEstHa <- resEstHa$parEst
  logLHa <- LoglikCalcCpp(parEstHa, DesignXHa, y_log, 
                          mm, conv = mytol, 
                          stepNR=stepsize[2])
  return(logLHa) # logLHa is a negative value
}

#-------------------------- Form new population
NewPop <- function(pop=NULL, fit=NULL, Pc=NULL, Pm=NULL, Mjmp=NULL, maxgen=NULL){ 
  
  IslandSize = length(pop)
  
  count <- 1
  repeat{
    flag <- c(0,0)
    ## step 3: select a pair of parents with linear ranking
    # the least fit (smallest value) has the lowest probability (0)
    # the best fit (largest value) has the largest probability
    probs0 <- 2*(IslandSize - rank(fit))/(IslandSize*(IslandSize-1))
    parent.i <- sample(1:IslandSize, 2, prob=probs0)
    
    # dad has the better fit (larger LogLik return value)
    tmp <- pop[parent.i]
    tmpp <- order(fit[parent.i], decreasing = T)
    dad <- tmp[tmpp[1]]
    mom <- tmp[tmpp[2]]
    
    ## step 4-1: crossover: round the average of dad and mom
    r41 <- runif(1, min=0, max=1)
    if (r41 <= Pc) {
      child <- round((dad+mom)/2)
    } else {
      child <- dad
      flag[1] <- 1
    }
    
    ## step 4-2: mutation: Jump to another day (period amount) randomly "+" or "-"
    r42 <- runif(1, min=0, max=1)
    if(r42 <= Pm){
      tmpsign <- sample(x=c(-1,1), size = 1)
      child  <- child + tmpsign*Mjmp
    }else{
      flag[2] <- 1
    }
    
    ## step 5: form new generation
    # steady state method: replace the least fit in the current pop if child is better
    # with child if child is better.
    
    flagsum <- flag[1] + flag[2] 
    # flag[1]=1 indicating no cross-over
    # flag[2]=1 indicating no mutation
    # flagsum < 2 indicating new individual produced
    
    if (flagsum<2){
      fit.child <- GAfit(tau=child, DesignXH0, y_log, mm, stepsize, mytol)
      if (min(fit)<fit.child) { 
        # indicating child is better than the worst one
        pp <- order(fit)[1]
        pop[pp] <- child
        fit[pp] <- fit.child
      }
      count <- count + 1
    }
    
    # check stopping
    if (count >= maxgen){break}
  }
  
  bestfit <- max(fit)
  bestchrom <- pop[which.max(fit)]
  
  res <- c(bestfit, bestchrom)
  return(res)
}

NewPop.parallel = function(pop, fit, Pc, Pm, Mjmp, maxgen){
  
  # each core calculates for each single island
  nTasks = dim(pop)[1]
  nCore = nTasks
  
  cat("\n          1. Multicores working, please wait ... \n")
  registerDoMC(cores = nCore)
  
  tim.start = Sys.time()
  FinalRes <- foreach(i=1:nCore, .combine = "rbind") %dopar%
    NewPop(pop[i,], fit[i,], Pc, Pm, Mjmp, maxgen)
  
  tim.end = Sys.time()
  
  cat("          2. nTasks =", nTasks, "\t nCore =", nCore,
      "\t Total Time =", difftime(tim.end, tim.start, units="hours"),
      "hours")
  
  return(FinalRes)
}

#-------------------------- Check Convergence

GAfit.one.job = function(job=i, nTasks=nTasks, nCore=nCore,
                         tau, DesignXH0, y_log, mm, stepsize, mytol){
  nSubtasks <- round(nTasks/nCore)
  RES <- rep(NA, nSubtasks)
  for(subid in 1:nSubtasks){
    tauEst = tau[(job-1)*nSubtasks+subid]
    RES[subid] = GAfit(tau=tauEst, DesignXH0, y_log, mm, stepsize, mytol)
  }
  return(RES)
}

GAfit.parallel = function(tau, DesignXH0, y_log, mm, stepsize, mytol, nCore){
  
  nTasks = length(tau)
  cat("\n          1. Multicores working, please wait ... \n")
  registerDoMC(cores = nCore)
  
  tim.start = Sys.time()
  FinalRes <- foreach(i=1:nCore, .combine = "c") %dopar%
    GAfit.one.job(job=i, nTasks=nTasks, nCore=nCore, 
                  tau, DesignXH0, y_log, mm, stepsize, mytol)
  tim.end = Sys.time()
  
  cat("          2. nTasks =", nTasks, "\t nCore =", nCore,
      "\t Total Time =", difftime(tim.end, tim.start, units="hours"),
      "hours")
  
  return(FinalRes)
}

#-------------------------- Genetic Algorithm Main Function
GA.MTM <- function(nIsland=NULL, IslandSize=NULL, Pc=NULL, Pm=NULL, maxMig=NULL, maxgen=NULL,
                   maxconv=NULL, Mjmp=NULL, tauClc=NULL, DesignXH0=NULL, 
                   y_log=NULL, mm=NULL, stepsize=NULL, mytol=NULL, 
                   parallel.initialize=FALSE, nCore=5, parallel.newpop=FALSE){
  
  cat("\n==== Step 1 and 2: Start to initialize population and calculate fitness ...")
  
  tim1 = Sys.time()
  # step 1: initial population
  pp <- sample(x=tauClc, size=nIsland*IslandSize, replace = F)
  # pp = floor(seq(from=min(tauClc), to=max(tauClc), length.out=nIsland*IslandSize))
  island <- matrix(pp, nrow=nIsland, ncol=IslandSize)
  islandFit <- matrix(0, nrow=nIsland, ncol=IslandSize)
  
  # step 2: Evaluate the fitness of the island individuals
  if(parallel.initialize==TRUE){
    islandFit = matrix(GAfit.parallel(tau=pp, DesignXH0, y_log, mm, stepsize, mytol, nCore),
                       nrow=nIsland, ncol=IslandSize)
  }else{
    for(i in 1:nIsland){
      for(j in 1:IslandSize){
        islandFit[i,j] <- GAfit(tau=island[i,j], DesignXH0, y_log, mm, stepsize, mytol)
      }
    }
  }
  tim2 = Sys.time()
  
  cat("\n      Completed within", tim2-tim1, "secs.")
  
  countMig <- 0
  overbest <- 0
  overbestChrom <- rep(0, times=maxMig)
  
  repeat{
    Bfit <- rep(0, times=nIsland)
    Bchrom <- rep(0, times=nIsland)
    
    cat("\n==== Step 3,4 and 5: No.", countMig+1, "migration ...")
    
    tim31 = Sys.time()
    # step 3,4,5: New Pop: select parents, crossover, mutation:
    
    if(parallel.newpop==TRUE){
      output = NewPop.parallel(pop=island, fit=islandFit, Pc, Pm, Mjmp, maxgen)
      Bfit = output[,1]
      Bchrom = output[,2]
    }else{
      for (k in 1:nIsland){
        output <- NewPop(pop=island[k,], fit=islandFit[k,], Pc, Pm, Mjmp, maxgen)
        Bfit[k] <- output[[1]]
        Bchrom[k] <- output[[2]]
      }
    }
    
    tim32 = Sys.time()
    
    cat("\n      New Pop completed within", tim32-tim31, "secs.")
    
    tim41 = Sys.time()
    ## Migration: After every maxgen, apply migration;
    #             replace least fit (smaller logL return value) in each island
    #             with the one randomly selected from the pool 
    #             of best island fit (largest logL return value)
    for (k in 1:nIsland){
      # select least fit individual for island k
      leastfitID = which.min(islandFit[k,])
      
      repeat{
        # random select one island as the source of best fit
        best <- sample(1:nIsland, 1)
        if (k!=best){
          island[k,leastfitID] <- Bchrom[best]
          islandFit[k,leastfitID] <- Bfit[best]
          break
        }
      }
    }
    tim42 = Sys.time()
    
    cat("\n      Migration completed within", tim42-tim41, "secs. \n")
    
    ## reform the Bestfit and associated chromsome
    for (k in 1:nIsland){
      p <- which.max(islandFit[k,])
      Bchrom[k] <- island[k,p]
      Bfit[k] <- islandFit[k,p]
    }
    
    countMig <- countMig + 1
    overbest[countMig] <- max(Bfit)
    overbestChrom[countMig] <- Bchrom[order(Bfit)[nIsland]]
    
    # Check convergence if all individuals has best fit
    if(countMig >= maxconv){
      tmp2 <- overbest[(countMig-maxconv+1):countMig]
      decision <- checkConv(tmp2, maxconv, mytol)
      if (decision == 1){
        bestfit <- overbest[countMig]
        bestchrom <- overbestChrom[countMig]
        break
      }
    }
    
    # Check stopping if reach max number of Migration
    if(countMig >= maxMig){
      bestfit <- overbest[countMig]
      bestchrom <- overbestChrom[countMig]
      break
    }
    # cat("\n      island:", island)
    cat("\n      Bchrom:", Bchrom)
    cat("\n      Chrosome:", overbestChrom)
  }
  
  output <- list(bestfit=bestfit, bestchrom=bestchrom, countMig=countMig)
  return(output)
}


#-------------------------- Check Convergence
checkConv <- function(a=NULL, maxconv=NULL, mytol=NULL){
  
  i<-1
  repeat{
    diff <- abs(a[i+1]-a[i])
    if (diff < mytol){
      i <- i+1
      if (i >= maxconv) {
        decision <- 1
        break
      }
    }
    else{
      decision <- 0
      break
    }
  }
  return(decision)
}


###################################
##### Main Script Starts Here #####
###################################

dat <- read.table(file = "5CanadaCloudData/hourly_day.txt")
dat1 <- dat[dat$V1>=1965 & dat$V1 <=1994,]

m <- as.numeric(as.character(dat1$V5))
dat1$V5 <-m

mm <- m
mm[m==1]  <- 0
mm[m==2]  <- 1
mm[m==3]  <- 1
mm[m==4]  <- 2
mm[m==5]  <- 2
mm[m==6]  <- 2
mm[m==7]  <- 3
mm[m==8]  <- 3
mm[m==9]  <- 4
mm[m==10] <- 4

dat1$V6 <- mm

#-------------------- Form Data Set -------------------#
Ti     <- length(mm)
K      <- max(mm) + 1
ss     <- 8*365

y_log <- matrix(0, nrow=Ti, ncol=K)
for (tt in 1:Ti) {y_log[tt,mm[tt]+1] <- 1}

# Design Matrix
AValue     <- rep(1,Ti)
TrendValue <- 1:Ti/Ti
BValue     <- cos(2*pi*(1:Ti)/ss)
DValue     <- sin(2*pi*(1:Ti)/ss)
DesignXH0  <- cbind(AValue, TrendValue, BValue, DValue)

# Computation setting
mytol    <- 0.00001
stepsize <- c(1,1)


#################################################################
#               Parameter Estimation under H0                   #
#################################################################
tim1 <- Sys.time()
parinitialH0  <- ParInitCalcu(y_log=y_log, m=mm, DesignX=DesignXH0)
resEstH0   <- SeqQuasiNWCpp(parinitialH0, y_log, mm, DesignXH0, 
                            stepsize, 1000, mytol)
logLH0     <- LoglikCalcCpp(resEstH0$parEst, DesignXH0, y_log, 
                            mm, conv = mytol, stepNR=stepsize[2])
par_est1 <- resEstH0$parEst
tim2 <- Sys.time()

### Parameter SE
mySEtheta <- sqrt(diag(solve(resEstH0$Hesstheta)))
mySExi <- sqrt(diag(solve(resEstH0$Hessxi)))
mySE1 <- c(mySEtheta, mySExi)

parMatrix1 <- cbind(par_est1[1:length(par_est1)],
                    par_est1[1:length(par_est1)] - 2*mySE1,
                    par_est1[1:length(par_est1)] + 2*mySE1)
colnames(parMatrix1) <- c("Param", "LowerB", "UpperB")
rownames(parMatrix1) <- c(paste("A", 1:(K-1), sep=""),
                          paste("beta", 1:(K-1), sep=""),
                          paste("B", 1:(K-1), sep=""),
                          paste("D", 1:(K-1), sep=""),
                          rep("xi", (K-1)^2))
round(parMatrix1, digits=4)

logLH0     <- LoglikCalcCpp(resEstH0$parEst, DesignXH0, y_log, 
                            mm, conv = mytol, stepNR=stepsize[2])
logLH0
# [1] -82159.04



#---- Tasks formulate
nCoreuse = 5

lc  <- floor(0.05*Ti) # 4382
uc  <- floor(0.95*Ti) # 83273
tmp <- lc:uc
tauClc <- lc:(lc+(floor(length(tmp)/nCoreuse) + 1)*nCoreuse-1) # From 730 to 83281
range(tauClc)
# [1]  4382 83281

# nIsland    <- nYears   # number of island
nIsland    <- 5
IslandSize <- 20   # number of individuals in each island
Pc         <- 0.9  # prob of corss-over 
Pm         <- 0.1  # prob of mutation
maxgen     <- 4    # for each subpopulation, after maxgen then apply migration
maxMig     <- 15   # maximum migration times
maxconv    <- 4    # number of consecutive migrations
Mjmp       <- 15    # Mutation Jump from original candidate time point
mytol      <- 0.00001

finaltim1 <- Sys.time()
result <- GA.MTM(nIsland, IslandSize, Pc, Pm, maxMig, maxgen, maxconv, Mjmp, tauClc,
                 DesignXH0, y_log, mm, stepsize, mytol, parallel.initialize=TRUE, 
                 nCore=nCoreuse, parallel.newpop=TRUE)
finaltim2 <- Sys.time()


# ==== Step 1 and 2: Start to initialize population and calculate fitness ...
# 1. Multicores working, please wait ... 
# 2. nTasks = 100 	 nCore = 5 	 Total Time = 0.6066463 hours
# Completed within 36.39879 secs.
# ==== Step 3,4 and 5: No. 1 migration ...
# 1. Multicores working, please wait ... 
# 2. nTasks = 5 	 nCore = 5 	 Total Time = 0.0892129 hours
# New Pop completed within 5.352857 secs.
# Migration completed within 0.0001299381 secs. 
# 
# Bchrom: 40613 39462 39462 39462 38066
# Chrosome: 39462 0 0 0 0 0 0 0 0 0 0 0 0 0 0
# ==== Step 3,4 and 5: No. 2 migration ...
# 1. Multicores working, please wait ... 
# 2. nTasks = 5 	 nCore = 5 	 Total Time = 0.0926318 hours
# New Pop completed within 5.557912 secs.
# Migration completed within 9.298325e-05 secs. 
# 
# Bchrom: 39462 39462 39462 39462 39462
# Chrosome: 39462 39462 0 0 0 0 0 0 0 0 0 0 0 0 0
# ==== Step 3,4 and 5: No. 3 migration ...
# 1. Multicores working, please wait ... 
# 2. nTasks = 5 	 nCore = 5 	 Total Time = 0.09259829 hours
# New Pop completed within 5.555903 secs.
# Migration completed within 0.0001058578 secs. 
# 
# Bchrom: 39462 39462 39462 39462 39462
# Chrosome: 39462 39462 39462 0 0 0 0 0 0 0 0 0 0 0 0
# ==== Step 3,4 and 5: No. 4 migration ...
# 1. Multicores working, please wait ... 
# 2. nTasks = 5 	 nCore = 5 	 Total Time = 0.0948703 hours
# New Pop completed within 5.692226 secs.
# Migration completed within 0.0001370907 secs. 
# 
# Bchrom: 40150 39462 39462 39462 39462
# Chrosome: 39462 39462 39462 40150 0 0 0 0 0 0 0 0 0 0 0
# ==== Step 3,4 and 5: No. 5 migration ...
# 1. Multicores working, please wait ... 
# 2. nTasks = 5 	 nCore = 5 	 Total Time = 0.09686014 hours
# New Pop completed within 5.811613 secs.
# Migration completed within 0.0001440048 secs. 
# 
# Bchrom: 40150 39462 39462 39462 39462
# Chrosome: 39462 39462 39462 40150 40150 0 0 0 0 0 0 0 0 0 0
# ==== Step 3,4 and 5: No. 6 migration ...
# 1. Multicores working, please wait ... 
# 2. nTasks = 5 	 nCore = 5 	 Total Time = 0.09187446 hours
# New Pop completed within 5.512473 secs.
# Migration completed within 0.0001070499 secs. 
# 
# Bchrom: 40150 40150 39462 39462 39462
# Chrosome: 39462 39462 39462 40150 40150 40150 0 0 0 0 0 0 0 0 0
# ==== Step 3,4 and 5: No. 7 migration ...
# 1. Multicores working, please wait ... 
# 2. nTasks = 5 	 nCore = 5 	 Total Time = 0.09247447 hours
# New Pop completed within 5.548474 secs.
# Migration completed within 0.0001499653 secs. 


finaltim2 - finaltim1
# Time difference of 1.25718 hours

result
# $bestfit
# [1] -82110.68
# 
# $bestchrom
# [1] 40150
# 
# $countMig
# [1] 7




tauEst = result$bestchrom
CpValue    <- c(rep(0, tauEst-1), rep(1, Ti-tauEst+1))
DesignXHa  <- cbind(DesignXH0, CpValue)
parinitialHa  <- ParInitCalcu(y_log=y_log, m=mm, DesignX=DesignXHa)
resEstHa   <- SeqQuasiNWCpp(parinitialHa, y_log, mm, DesignXHa, 
                            stepsize, 1000, mytol)
parHa      <- resEstHa$parEst
logLHa     <- LoglikCalcCpp(parHa, DesignXHa, y_log, mm, 
                            conv = mytol, stepNR=stepsize[2])
logLHa
# [1] -82110.68
LambdaStat.Est <- -2*(logLH0 - logLHa)
LambdaStat.Est
# [1] 96.71995

# estimated changepoint related date
dat1[tauEst,]
# 43078 1978  9 28 14  7  3