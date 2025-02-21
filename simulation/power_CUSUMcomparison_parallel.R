rm(list=ls())

# set library path
.libPaths(c("userpackages", .libPaths()))

library(VGAM)
library(foreach)
library(doMC)
library(Rcpp)
library(RcppArmadillo)

sourceCpp("src/QuasiNewton.cpp")
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

#-------------------------- Raw Estimation of Dependence Model Parameters
initial_xi <- function(m=NULL, K=NULL, avg_delta=NULL){
  
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
ParInitCalcu <- function(y_log=NULL, m=NULL, DesignX=NULL, DesignXH0=NULL){
  
  require(VGAM)
  
  K <- dim(y_log)[2]
  len_mean <- dim(DesignX)[2]
  if(len_mean == 1){
    fit.d <- try(vglm(y_log ~ 1, cumulative(parallel=FALSE)))
    MarginalParam = as.vector(fit.d@coefficients)
  }else{
    fit.d <- try(vglm(y_log ~ DesignX[,2:dim(DesignX)[2]], cumulative(parallel=FALSE)))
    if(class(fit.d) == "try-error"){
      fit.d <- vglm(y_log ~ DesignXH0[,2:dim(DesignXH0)[2]], cumulative(parallel=FALSE))
      MarginalParam = c(as.vector(fit.d@coefficients), rep(0, K-1))
    }else{
      MarginalParam = as.vector(fit.d@coefficients)
    }
  }
  avg_delta <- colMeans(fit.d@fitted.values)
  XiInit <- initial_xi(m=m, K=K, avg_delta=avg_delta)
  
  parinit <- c(MarginalParam, as.vector(XiInit))
  return(parinit)
}

#-------------------------- MTM LRT fitness function
MTMGAfit = function(chromosome, plen=0, y_log, mm, DesignXH0, logL0, stepsize, maxit, mytol){
  
  # tau here only contain the changepoint locations
  m = chromosome[1]
  tau = chromosome[2]
  
  N = NROW(y_log) #length of the series
  
  if(tau < floor(0.05*N) | tau > ceiling(0.95*N)){
    Lambdatau = -2*(logL0 + 1e6)
  }else{
    CpValue = c(rep(0, tau-1), rep(1, N-tau+1))
    DesignXHa  = cbind(DesignXH0, CpValue)
    parinitialHa  = ParInitCalcu(y_log=y_log, m=mm, DesignX=DesignXHa, DesignXH0=DesignXH0)
    resEstHa = SeqNewtonCpp(parinitialHa, y_log, mm, DesignXHa, stepsize, maxit, mytol)
    parEstHa = resEstHa$parEst
    logLHa = LoglikCalcCpp(parEstHa, DesignXHa, y_log, mm, conv = mytol, stepNR=stepsize[2])
    Lambdatau = -2*(logL0 - logLHa)
  }
  
  return(-Lambdatau)
}

one.job = function(job=NULL, nTasks=NULL, nCore=NULL, condition=NULL){
  
  #----------- Simulation parameters ------------------#
  DeltaT = condition$Delta
  propT  = condition$prop
  nY     = condition$nY
  # m      = condition$m
  
  ##### simulation parameters
  ss     = 365
  N      = nY*ss
  K      = 5
  tauT   = round(propT*N)
  
  # switch (as.character(nY),
  #   "3"  = {CV=19.836},
  #   "5"  = {CV=19.599},
  #   "10" = {CV=19.798},
  #   "20" = {CV=19.794},
  #   "30" = {CV=19.453}
  # )
  
  A     <- c(-1.3863, -0.4055, 0.4055, 1.3863) # equal space
  Alpha <- c( 0.1,  0.1,  0.1,   0.1)
  B     <- c(-0.1, -0.2, -0.15, -0.3)
  D     <- c( 0.2,  0.1,  0.15,  0.3)
  # Delta <- c(rep(0, m), rep(DeltaT, K-1-m))
  # Delta <- c( 0, DeltaT, DeltaT, DeltaT)
  Delta = c(0, DeltaT, 0, DeltaT)
  if(DeltaT == 0){
    parT  <- c(A,
               Alpha,
               B, D, 
               5.9,  4.0,  2.6,  1.7,
               4.7,  3.6,  2.4,  1.7,
               3.2,  2.2,  1.5,  1.2,
               1.8,  1.6,  0.9,  0.6)
               # 3.9,  2.7,  1.7,  1.3,
               # 2.3,  1.5,  1.0,  0.5,
               # 2.1,  1.1,  0.9,  0.7,
               # 1.2,  0.9,  0.6,  0.5)
    # 0, 0, 0, 0, 
    # 0, 0, 0, 0,
    # 0, 0, 0, 0,
    # 0, 0, 0, 0)
  }else{
    parT  <- c(A,
               Alpha,
               B, D, 
               Delta,
               5.9,  4.0,  2.6,  1.7,
               4.7,  3.6,  2.4,  1.7,
               3.2,  2.2,  1.5,  1.2,
               1.8,  1.6,  0.9,  0.6)
               # 3.9,  2.7,  1.7,  1.3,
               # 2.3,  1.5,  1.0,  0.5,
               # 2.1,  1.1,  0.9,  0.7,
               # 1.2,  0.9,  0.6,  0.5)
    # 0, 0, 0, 0, 
    # 0, 0, 0, 0,
    # 0, 0, 0, 0,
    # 0, 0, 0, 0)
  }
  
  
  ##### Design Matrix
  AValue     = rep(1, N)
  TrendValue = (1:N)/N
  BValue     = cos(2*pi*(1:N)/ss)
  DValue     = sin(2*pi*(1:N)/ss)
  CpValue    = c(rep(0, tauT-1), rep(1, N-tauT+1))
  # DesignXT   = cbind(AValue, TrendValue, BValue, DValue, CpValue)
  # DesignXH0  = cbind(AValue, TrendValue, BValue, DValue)
  
  if(DeltaT == 0){
    DesignXT   = cbind(AValue, TrendValue, BValue, DValue)
  }else{
    DesignXT   = cbind(AValue, TrendValue, BValue, DValue, CpValue)
  }
  
  DesignXH0  = cbind(AValue, TrendValue, BValue, DValue)
  
  ##### computation setting
  mytol = 1e-6
  stepsize = c(1, 1)
  maxit = 1000
  
  nSubtasks = round(nTasks/nCore)
  RES = data.frame(job=job, task=1:nSubtasks,
                   tim=rep(NA, nSubtasks),
                   seed=rep(NA, nSubtasks),
                   fit.mtm=rep(NA, nSubtasks),
                   tauhat.mtm=rep(NA, nSubtasks),
                   CpValY.cusum=rep(NA, nSubtasks),
                   CpLocY.cusum=rep(NA, nSubtasks),
                   CpValZ.cusum=rep(NA, nSubtasks),
                   CpLocZ.cusum=rep(NA, nSubtasks))
  # param.mat.mtm = matrix(NA, nrow=nSubtasks, ncol=NCOL(DesignXT)*(K-1)+(K-1)^2)
  # ub.mat.mtm = lb.mat.mtm = param.mat.mtm
  # colnames(param.mat.mtm) = paste0("param", 1:NCOL(param.mat.mtm))
  # colnames(lb.mat.mtm) = paste0("lb", 1:NCOL(lb.mat.mtm))
  # colnames(ub.mat.mtm) = paste0("ub", 1:NCOL(ub.mat.mtm))
  
  for(subid in 1:nSubtasks){
    
    tim.start = Sys.time()
    
    myseed = 0 + job*nSubtasks + subid
    cat("\n No.job", job, "task", subid, "start at seed", myseed, "!")
    set.seed(myseed)
    mm = MTMSimCpp(par=parT, K=K, DesignX=DesignXT, m0=1, conv=mytol, stepNR=stepsize[2])
    
    y_log <- matrix(0, nrow=N, ncol=K)
    for (tt in 1:N) {y_log[tt, mm[tt]+1] <- 1}
    
    ##### CUSUM detection #####
    res_sim = list(X_hour=mm, X_hour_wide=y_log)
    FitRes11 = try(UniCpDetect(res_sim = res_sim, DesignXEst=DesignXH0, stepsize = 1, BoundAdjust = -0.5))
    if(class(FitRes11) == "try-error"){
      CpValY = NA
      CpLocY = NA
      CpValZ = NA
      CpLocZ = NA
    }else{
      CpValY = FitRes11$CpValY
      CpLocY = FitRes11$CpLocY
      CpValZ = FitRes11$CpValZ
      CpLocZ = FitRes11$CpLocZ
    }
    ##### MTM detection #####
    
    ##  Model fit under H0
    parinitialH0  = ParInitCalcu(y_log=y_log, m=mm, DesignX=DesignXH0)
    resEstH0      = try(SeqQuasiNWCpp(parinitialH0, y_log, mm, DesignXH0, stepsize, maxit, mytol))
    if(class(resEstH0) == "try-error"){
      Lambda_max = NA
      tauhat = NA
    }else{
      logLH0        = LoglikCalcCpp(resEstH0$parEst, DesignXH0, y_log, mm, conv = mytol, stepNR=stepsize[2])
      
      ##  Maximum LRT searching under Ha
      tauClc = floor(0.05*N):ceiling(0.95*N)
      logLHaClc = rep(NA, length(tauClc))
      for(qqq in 1:length(tauClc)){
        tmptau = tauClc[qqq]
        CpValue = c(rep(0, tmptau-1), rep(1, N-tmptau+1))
        DesignXHa  = cbind(DesignXH0, CpValue)
        parinitialHa  = ParInitCalcu(y_log=y_log, m=mm, DesignX=DesignXHa, DesignXH0=DesignXH0)
        resEstHa = try(SeqNewtonCpp(parinitialHa, y_log, mm, DesignXHa, stepsize, maxit, mytol))
        if(class(resEstHa) == "try-error"){
          logLHaClc[qqq] = NA
        }else{
          parEstHa = as.vector(resEstHa$parEst)
          logLHaClc[qqq] = LoglikCalcCpp(parEstHa, DesignXHa, y_log, mm, conv = mytol, stepNR=stepsize[2])
        }
      }
      LambdaClc <- -2*(logLH0 - logLHaClc)
      if(all(is.na(LambdaClc))){
        Lambda_max = NA
        tauhat = NA
      }else{
        Lambda_max = max(LambdaClc, na.rm = TRUE)
        tauhat = tauClc[which.max(LambdaClc)]    
      }
    }
    
    # ## Model refit with detected changepoint under Ha
    # CpValue       = c(rep(0, tauhat-1), rep(1, N-tauhat+1))
    # DesignXEst    = cbind(DesignXH0, CpValue)
    # parinitialEst = ParInitCalcu(y_log=y_log, m=mm, DesignX=DesignXEst)
    # resEst        = SeqQuasiNWCpp(parinitialEst, y_log, mm, DesignXEst, stepsize, maxit, mytol)
    # paramEst      = as.vector(resEst$parEst)
    # ## CI Bound
    # prop_sigma = c(sqrt(diag(solve(resEst$Hesstheta))), sqrt(diag(solve(resEst$Hessxi))))
    # ub = paramEst + 1.96*prop_sigma
    # lb = paramEst - 1.96*prop_sigma
    
    tim.end <- Sys.time()
    timeused <- difftime(tim.end, tim.start, units="secs")
    cat("\n No.job", job, "task", subid,
        "completed within", timeused, "at seed", myseed, "!")
    
    ##### Result return #####
    RES$tim[subid] = as.numeric(timeused)
    RES$seed[subid] = myseed
    RES$fit.mtm[subid] = Lambda_max
    RES$tauhat.mtm[subid] = tauhat
    RES$CpValY.cusum[subid] = CpValY
    RES$CpLocY.cusum[subid] = CpLocY
    RES$CpValZ.cusum[subid] = CpValZ
    RES$CpLocZ.cusum[subid] = CpLocZ
    
    # param.mat[subid, ] = paramEst
    # lb.mat[subid, ] = lb
    # ub.mat[subid, ] = ub
  }
  
  # RES = cbind(RES, lb.mat, ub.mat)
  
  return(RES)
  
}

useMultiCore <- function(nTasks=NULL, nCore=NULL, condition=NULL){
  cat("Multicores working, please wait ... \n")
  registerDoMC(cores = nCore)
  tim.start = Sys.time()
  FinalRes = foreach(i=1:nCore, .combine = "rbind") %dopar%
    one.job(job=i, nTasks=nTasks, nCore=nCore, condition=condition)
  tim.end = Sys.time()
  cat("Done.\n")
  cat("\n\n nTasks =", nTasks, "\t nCore =", nCore,
      "\t Aveg. time =", mean(FinalRes$tim),
      "\t Total Time =", difftime(tim.end, tim.start, units="hours"),
      "hours \n\n")
  return(FinalRes)
}

######################################################
#-------------- Simulation Setting ------------------#
######################################################

##### parameters
allDelta = c(-0.3, 0.3, 0.5)
allprop = c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
allnY    = 3
# number of DeltaT = 0
# allm = c(0, 1, 2, 3)
eg = expand.grid(Delta=allDelta, prop=allprop, nY=allnY, stringsAsFactors=FALSE)
eg$ID = 84+1:NROW(eg)


### Warning: !!!!!!!!!!!
# check your computer available cores before run code below!!!!!!!!!!!!!
detectCores()
### Warning: !!!!!!!!!!!
# You have to set "nCores < number of Cores available in your computer" !!!!!
### Warning: !!!!!!!!!!!
# be careful of your computer memory !!!!!!!!!!!


# number of cores you for computation
nCore = 48
# number of simulations
nTasks = 1008

datafile.dir = "power/strpos/"

for(iscen in 1:dim(eg)){
  
  condition = eg[iscen,]
  FinalRes = useMultiCore(nTasks=nTasks, nCore=nCore, condition=condition)
  
  # save(FinalRes, condition, file=paste0(datafile.dir, "multicp_configID", condition$ID, "_strpos_model2.RData"))
}

