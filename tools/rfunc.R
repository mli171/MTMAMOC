
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

#-------------------------- category aggregation for application in cloud cover data
catg_agg = function(m){
  
  mm = rep(NA, length(m))
  
  mm[m==0]  <- 0
  mm[m==1]  <- 1
  mm[m==2]  <- 1
  mm[m==3]  <- 1
  mm[m==4]  <- 2
  mm[m==5]  <- 2
  mm[m==6]  <- 3
  mm[m==7]  <- 3
  mm[m==8]  <- 3
  mm[m==9]  <- 3
  mm[m==10] <- 4
  
  if (any(is.na(mm))) {
    cat("\n Warning: NA exist in the series! \n")
  }
  
  return(mm)
}

