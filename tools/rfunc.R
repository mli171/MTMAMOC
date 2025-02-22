
#-------------------------- Raw Estimation of Dependence Model Parameters
initial_xi = function(m=NULL, avg_delta=NULL){
  
  # Naive Transition Matrix Estimation
  ct = matrix(0, K, K)  # previous*current
  for(i in 1:K){
    p = which(m==(i-1)) # where is i-1 category
    for (j in 1:length(p)){
      ct[i,m[p[j]+1]+1] = ct[i,m[p[j]+1]+1] + 1
    }
  }
  pc_avg = ct/rowSums(ct)   # pc_avg: j*k
  
  if(any(pc_avg == 0)){
    cat("Zero Exist in Transition Matrix!!!\n")
    pc_avg[which(pc_avg==0)] = 0.00001
    # break
  }
  xi = matrix(0, K-1, K-1) # previous*current
  for(j in 1:(K-1)){
    for(k in 1:(K-1)){
      # calculate xi
      xi[j,k] = log(pc_avg[j,k]/pc_avg[j,K])-avg_delta[k]
    }
  }
  return(t(xi)) # res: current*previous
}

#-------------------------- Initial Values for All Model Parameters
ParInitCalcu = function(y_log=NULL, m=NULL, DesignX=NULL){
  
  require(VGAM)
  
  K = dim(y_log)[2]
  len_mean = dim(DesignX)[2]
  if(len_mean == 1){
    fit.d = vglm(y_log ~ 1, cumulative(parallel=FALSE))
  }else{
    fit.d = vglm(y_log ~ DesignX[,2:dim(DesignX)[2]], cumulative(parallel=FALSE))
  }
  avg_delta = colMeans(fit.d@fitted.values)
  XiInit = initial_xi(m=m, avg_delta=avg_delta)
  
  parinit = c(as.vector(fit.d@coefficients), as.vector(XiInit))
  
  return(parinit)
}

#-------------------------- periodogram frequency plot
myperiodogram = function(mydat, candPeriod=NULL, TableReturn=FALSE){
  
  require(astsa)
  
  u = factor(mydat) # first, input the data as factors and then 
  x = model.matrix(~u-1)[,1:3] # make an indicator matrix
  Var = var(x) # var-cov matrix
  xspec = mvspec(x, spans=c(7,7), plot=FALSE)
  fxxr = Re(xspec$fxx)  # fxxr is real(fxx)
  # compute Q = Var^-1/2
  ev = eigen(Var)
  Q = ev$vectors%*%diag(1/sqrt(ev$values))%*%t(ev$vectors)
  # compute spec envelope and scale vectors
  num = xspec$n.used # sample size used for FFT
  nfreq = length(xspec$freq) # number of freqs used
  specenv = matrix(0,nfreq,1) # initialize the spec envelope 
  beta = matrix(0,nfreq,3) # initialize the scale vectors 
  for (k in 1:nfreq){
    ev = eigen(2*Q%*%fxxr[,,k]%*%Q/num, symmetric=TRUE)
    specenv[k] = ev$values[1] # spec env at freq k/n is max evalue 
    b = Q%*%ev$vectors[,1] # beta at freq k/n
    beta[k,] = b/sqrt(sum(b^2)) # helps to normalize beta
  }
  # output and graphics
  frequency = xspec$freq
  plot(frequency, 100*specenv, type="l", xlim=c(0,0.05),
       xlab="Frequency", ylab="Spectral Envelope (%)") 
  # add significance threshold to plot
  m = xspec$kernel$m
  etainv = sqrt(sum(xspec$kernel[-m:m]^2)) 
  thresh=100*(2/num)*exp(qnorm(.9999)*etainv)
  abline(h=thresh, lty=6, col=4)
  # text(0.1-0.02, thresh+0.02, round(thresh,3), xpd=TRUE, adj=1, col=4)
  
  if(!is.null(candPeriod)){
    for(i in 1:length(candPeriod)){
      abline(v=1/candPeriod[i], col="red", lty="dashed")
      text(1/candPeriod[i]+0.007, 100*max(specenv)-0.01, paste0("Frequence = 1/", candPeriod[i]),xpd=TRUE,adj=1, col="red")
    }
  }
  
  # details
  if(TableReturn){
    output = cbind(frequency, specenv, beta)
    colnames(output) = c("freq","specenv", "A", "C", "G") 
    return(output)
  }
  
}

#-------------------------- category aggregation for application in cloud cover data
catg_agg = function(m){
  
  mm = rep(NA, length(m))
  
  mm[m==0]  = 0
  mm[m==1]  = 1
  mm[m==2]  = 1
  mm[m==3]  = 1
  mm[m==4]  = 2
  mm[m==5]  = 2
  mm[m==6]  = 2
  mm[m==7]  = 3
  mm[m==8]  = 3
  mm[m==9]  = 3
  mm[m==10] = 4
  
  if (any(is.na(mm))) {
    cat("\n Warning: NA exist in the series! \n")
  }
  
  return(mm)
}


#-------------------------- Genetic Algorithm Fitness function
GAfit <- function(chromosome=NULL, plen = 0, DesignXH0=NULL, y_log=NULL, 
                  mm=NULL, stepsize=NULL, mytol=NULL, logLH0=NULL){
  
  m = chromosome[1]
  tau = chromosome[2]
  
  Ti = dim(DesignXH0)[1]
  if(tau < floor(0.05*Ti) | tau > ceiling(0.95*Ti)){
    Lambdatau = -2*(logL0 + 1e6)
  }else{
    CpValue    <- c(rep(0, tau-1), rep(1, Ti-tau+1))
    DesignXHa  <- cbind(DesignXH0, CpValue)
    parinitialHa  <- ParInitCalcu(y_log=y_log, m=mm, DesignX=DesignXHa)
    resEstHa <- SeqNewtonCpp(parinitialHa, y_log, mm, DesignXHa, 
                             stepsize, 1000, mytol)
    parEstHa <- resEstHa$parEst
    logLHa <- LoglikCalcCpp(parEstHa, DesignXHa, y_log, 
                            mm, conv = mytol, 
                            stepNR=stepsize[2])
    Lambdatau = -2*(logLH0 - logLHa)
  }
  
  return(-Lambdatau)
}
