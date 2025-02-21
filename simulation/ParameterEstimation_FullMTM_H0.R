# Parallelized perform nTasks simulations to check the parameter estimation 
#    performance by (bias, sd, mse)

library(VGAM)
library(foreach)
library(doMC)
library(Rcpp)
library(RcppArmadillo)

sourceCpp("src/QuasiNewton.cpp")

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

#-------------------------- Server Single job
one.job.param <- function(job=NULL, nTasks=NULL, nCore=NULL, iii=NULL,
                          parT=NULL, Ti=NULL, K=NULL,
                          DesignXT=NULL, DesignXEst=NULL,
                          stepsize=NULL, mytol=NULL){
  
  len_par <- length(parT)
  nSubtasks <- ceiling(nTasks/nCore)
  RES <- data.frame(job=job, task=1:nSubtasks, tim=rep(NA, nSubtasks),
                    seed=rep(NA, nSubtasks), 
                    LogLT=rep(NA, nSubtasks), LogLEst=rep(NA, nSubtasks))
  ParamMatrix <- matrix(0, nrow=nSubtasks, ncol=len_par)
  RES <- cbind(RES, ParamMatrix)
  
  for(subid in 1:nSubtasks){
    
    tim1 <- Sys.time()
    #################################################################
    #                          Simulation                           #
    #################################################################
    # myseed <- 20191111 + job*nSubtasks + subid + iii*100000
    # set.seed(myseed)
    mm <- MTMSimCpp(par=parT, K=K, DesignX=DesignXT, m0=1, 
                    conv=mytol, stepNR=stepsize[2])
    y_log <- matrix(0, nrow=Ti, ncol=K)
    for (tt in 1:Ti) {y_log[tt, mm[tt]+1] <- 1}
    LogLT <- LoglikCalcCpp(parT, DesignXT, y_log, mm = mm, 
                           conv = mytol, stepNR=stepsize[2])
    #################################################################
    #                     Parameter Estimation                      #
    #################################################################
    parinitialEst  <- ParInitCalcu(y_log=y_log, m=mm, DesignX=DesignXEst)
    resEst <- SeqNewtonCpp(parinitialEst, y_log, mm, DesignXEst, 
                           stepsize, 1000, mytol)
    parEst <- as.vector(resEst$parEst)
    LogLEst <- LoglikCalcCpp(parEst, DesignXEst, y_log, mm, 
                             conv = mytol, stepNR=stepsize[2])
    
    tim2 <- Sys.time()
    timecost <- difftime(tim2, tim1, units="mins")
    cat("\n No.job", job, "task", subid, 
        "within", nSubtasks, "(subtasks) completed within", 
        timecost, "mins!")  
    #################################################################
    #                         Result Return                         #
    #################################################################
    RES$tim[subid]     <- timecost
    RES$seed[subid]    <- 0            # indicate complete random
    RES$LogLT[subid]   <- LogLT
    RES$LogLEst[subid] <- LogLEst
    RES[subid, (dim(RES)[2]-len_par + 1):dim(RES)[2]] <- parEst
  }
  
  return(RES)
}

#-------------------------- Parallel computing implementation
useMultiCore <- function(nTasks=NULL, nCore=NULL, iii=NULL,
                         parT=NULL, Ti=NULL, K=NULL,
                         DesignXT=NULL, DesignXEst=NULL,
                         stepsize=NULL, mytol=NULL){
  
  cat("Multicores working, please wait ... \n")
  registerDoMC(cores = nCore)
  tim.start = Sys.time()
  FinalRes <- foreach(i=1:nCore, .combine = "rbind") %dopar%
    one.job.param(job=i, nTasks=nTasks, nCore=nCore, iii=iii,
                  parT=parT, Ti=Ti, K=K,
                  DesignXT=DesignXT, DesignXEst=DesignXEst,
                  stepsize=stepsize, mytol=mytol)
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

#----- Simulation Parameters (K=5)
Ti     <- 600
K      <- 3
ss     <- 3

A     <- c(-0.5,  0.8)
Alpha <- c(-0.8, -0.5)
B     <- c(-0.4, -0.5)
D     <- c( 0.9,  0.9)
parT  <- c(A, Alpha, B, D, 1.7,  0.9,  0.7,  0.4)

#----- computation setting
maxit    <- 1000
stepsize <- c(1, 0.5)
mytol    <- 1e-05

nCore  <- 5
nTasks <- 1000

#----- Design Matrix
AValue     <- rep(1, Ti)
TrendValue <- (1:Ti)/Ti
BValue     <- cos(2*pi*(1:Ti)/ss)
DValue     <- sin(2*pi*(1:Ti)/ss)

DesignXT   <- cbind(AValue, TrendValue, BValue, DValue)
DesignXEst <- cbind(AValue, TrendValue, BValue, DValue)

for(iii in 1:nSim){
  ParamFinal <- useMultiCore(nTasks=nTasks, nCore=nCore, iii=iii,
                             parT=parT, Ti=Ti, K=K,
                             DesignXT=DesignXT, DesignXEst=DesignXEst,
                             stepsize=stepsize, mytol=mytol)
}

# number of NA's 
length(which(is.na(ParamFinal[,7])))/nTasks

# mean
colMeans(ParamFinal[,7:dim(ParamFinal)[2]])
# bias
colMeans(ParamFinal[,7:dim(ParamFinal)[2]]) - parT
# sd
apply(ParamFinal[,7:dim(ParamFinal)[2]], 2, sd)
# mse
sum((ParamFinal[,7:dim(ParamFinal)[2]] - colMeans(ParamFinal[,7:dim(ParamFinal)[2]]))^2)/nTasks



