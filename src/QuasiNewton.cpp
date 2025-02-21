// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>

using namespace Rcpp;
using namespace arma;

int len_mean, len_par, len_Marg, len_Dep;
int Ti, K, i, k, kp, kk, j, a, b, ap, bp, p, q, g, gg, t, lag;


// [[Rcpp::export]]
arma::vec MTMSimCpp(arma::vec par, int K, arma::mat DesignX, int m0, double conv, double stepNR){
  
  // extract parameters
  Ti        = DesignX.n_rows;
  len_mean  = DesignX.n_cols; 
  len_par   = par.size();
  
  //-----------  Allocation Memory  -----------//
  int iconv_nr;
  double diff_nr;
  
  arma::vec m(Ti, fill::zeros);
  arma::vec delta0(K-1, fill::zeros);
  arma::vec delta1(K-1, fill::zeros);
  arma::vec yvec(K, fill::zeros);
  arma::vec f(K, fill::zeros);;
  
  arma::mat theta(K-1, len_mean, fill::zeros);
  arma::mat xi(K-1, K-1, fill::zeros);
  arma::mat eta(Ti, K-1, fill::zeros);
  arma::mat gamm(Ti, K-1, fill::zeros);
  arma::mat pm(Ti, K, fill::zeros);
  
  arma::mat pc_nr(K,K,fill::zeros);
  arma::mat pc_tmp(K-1,K,fill::zeros);
  arma::mat dfsum(K-1,K-1,fill::zeros);
  arma::mat pt(1, 1, fill::zeros);
  
  //-----------  Allocation Memory  -----------//
  
  
  for(i=0;i<len_mean;i++){
    theta.col(i) = par.subvec(i*(K-1),(i+1)*(K-1)-1);
  }
  for(i=0;i<K-1;i++){
    xi.col(i) = par.subvec(len_mean*(K-1)+i*(K-1), \
           len_mean*(K-1)+i*(K-1)+(K-1)-1);
  }
  
  //----------- Marginal Parameter -----------//
  eta = DesignX*trans(theta);
  gamm = exp(eta)/(1 + exp(eta));
  
  // pm.zeros();
  for(i=0;i<K-1;i++){
    if(i==0){
      pm.col(i) = gamm.col(i);
    }else{
      pm.col(i) = gamm.col(i) - gamm.col(i-1);
    }
  }
  pm.col(K-1) = 1- gamm.col(K-2);
  
  // Check marginal probability positive
  while(bool positive_prob=any(vectorise(pm)<0)) {
    Rcout<<"Cpp Notice: Marginal Probability less than 0 !!!"<<std::endl;
    break;
  }
  
  m(1) = m0;
  for(t=1;t<Ti;t++){
    
    iconv_nr = 0;
    diff_nr  = 0.01;
    
    if(t>1){delta0 = delta1;}
    
    while((iconv_nr<1000) & (diff_nr>conv)){
      
      // Conditional Probability: current(row)*previous(col)
      pc_nr.zeros();
      pc_tmp.zeros();
      for(k=0;k<K-1;k++){
        for(j=0;j<K;j++){
          // indicator vector in conditional model
          yvec.zeros();
          yvec(j) = 1;
          pc_tmp(k,j) = delta0(k)+ \
            sum(trans(xi.row(k))%yvec.subvec(0,K-2));
        }
      }
      
      // log link
      for(j=0;j<K;j++){
        for(k=0;k<K-1;k++){
          pc_nr(k,j) = exp(pc_tmp(k,j))/(1+sum(exp(pc_tmp.col(j))));
        }
        pc_nr(K-1,j) = 1/(1+sum(exp(pc_tmp.col(j))));
      }
      
      dfsum.zeros();
      for(k=0;k<K-1;k++){
        for(j=0;j<K-1;j++){
          if(j==k){
            dfsum(k,j) = sum(pc_nr.row(k)%(1-pc_nr.row(k))%pm.row(t-1));
          }else{
            dfsum(k,j) = -sum(pc_nr.row(k)%pc_nr.row(j)%pm.row(t-1));
          }
        }
      }
      
      f = pc_nr*trans(pm.row(t-1)) - trans(pm.row(t));
      delta1 = delta0 - solve(dfsum,f.subvec(0,K-2),solve_opts::fast)*stepNR;
      diff_nr = sqrt(sum((delta1-delta0)%(delta1-delta0)));
      delta0 = delta1;
      iconv_nr = iconv_nr + 1;
    }
    
    if(iconv_nr>1000){
      Rcout << "Delta is not converged at time=" << t << std::endl;
    }
    
    // Generate Count Data Randomly //
    pt = randu(1);
    uvec tmp = find(pt(0,0)>cumsum(pc_nr.col(m(t-1))));
    m(t) = tmp.size();
  }
  
  return(m);
}

// [[Rcpp::export]]
void NRDeltaCpp(arma::mat xi, arma::mat& pm, arma::mat& DesignX, arma::mat& delta, arma::cube& pc, int K, double conv, double stepNR){
  
  // extract parameters
  Ti        = DesignX.n_rows;
  len_mean  = DesignX.n_cols; 
  
  int iconv_nr=0; 
  double diff_nr=0.01;
  
  arma::vec delta0(K-1, fill::zeros);
  arma::vec deltatmp(K-1, fill::zeros);
  arma::vec delta1(K-1, fill::zeros);
  // arma::vec yvec(K, fill::zeros);
  arma::vec f(K, fill::zeros);;
  
  mat pc_nr(K,K,fill::zeros);
  // mat delta(Ti,K-1,fill::zeros);
  mat pc_tmp(K-1,K,fill::zeros);
  mat dfsum(K-1,K-1,fill::zeros);
  
  // t = 1
  while((iconv_nr<1000) & (diff_nr>conv)){
    pc_nr.zeros();
    pc_tmp.zeros();
    for(j=0;j<K-1;j++){
      pc_tmp.col(j) = delta0+xi.col(j);
    }
    pc_tmp.col(K-1) = delta0;
    for(j=0;j<K;j++){
      pc_nr.col(j).subvec(0,K-2) = exp(pc_tmp.col(j))/(1+sum(exp(pc_tmp.col(j))));
      pc_nr(K-1,j) = 1/(1+sum(exp(pc_tmp.col(j))));
    }
    dfsum.zeros();
    for(k=0;k<K-1;k++){
      for(j=k;j<K-1;j++){
        if(j==k){
          dfsum(k,j) = sum(pc_nr.row(k)%(1-pc_nr.row(k))%pm.row(0));
        }else{
          dfsum(k,j) = -sum(pc_nr.row(k)%pc_nr.row(j)%pm.row(0));
        }
      }
    }
    dfsum += trans(dfsum);
    //Since double count the diagonal elements
    dfsum.diag() = dfsum.diag()/2;
    f = pc_nr*trans(pm.row(0)) - trans(pm.row(1));
    deltatmp = delta0 - solve(dfsum,f.subvec(0,K-2),solve_opts::fast)*stepNR;
    // deltatmp = delta0 - pinv(dfsum)*f.subvec(0,K-2)*stepNR;
    
    //-------------------------------------------------------------------------------------//
    deltatmp = (deltatmp + delta0)/2;
    pc_nr.zeros();
    pc_tmp.zeros();
    for(j=0;j<K-1;j++){
      pc_tmp.col(j) = deltatmp+xi.col(j);
    }
    pc_tmp.col(K-1) = deltatmp;
    for(j=0;j<K;j++){
      pc_nr.col(j).subvec(0,K-2) = exp(pc_tmp.col(j))/(1+sum(exp(pc_tmp.col(j))));
      pc_nr(K-1,j) = 1/(1+sum(exp(pc_tmp.col(j))));
    }
    dfsum.zeros();
    for(k=0;k<K-1;k++){
      for(j=k;j<K-1;j++){
        if(j==k){
          dfsum(k,j) = sum(pc_nr.row(k)%(1-pc_nr.row(k))%pm.row(0));
        }else{
          dfsum(k,j) = -sum(pc_nr.row(k)%pc_nr.row(j)%pm.row(0));
        }
      }
    }
    dfsum += trans(dfsum);
    //Since double count the diagonal elements
    dfsum.diag() = dfsum.diag()/2;
    f = pc_nr*trans(pm.row(0)) - trans(pm.row(1));
    delta1 = delta0 - solve(dfsum,f.subvec(0,K-2),solve_opts::fast)*stepNR;
    // delta1 = delta0 - pinv(dfsum)*f.subvec(0,K-2)*stepNR;
    
    //-------------------------------------------------------------------------------------//
    diff_nr = sqrt(sum((delta1-delta0)%(delta1-delta0)));
    delta0 = delta1;
    iconv_nr = iconv_nr + 1;
  }
  pc.slice(1) = pc_nr;
  delta.row(1) = trans(delta0);
  
  // t > 1
  for(t=2;t<Ti;t++){
    iconv_nr = 0;
    diff_nr  = 0.01;
    
    while((iconv_nr<1000) & (diff_nr>conv)){
      // Conditional Probability: current(row)*previous(col)
      pc_nr.zeros();
      pc_tmp.zeros();
      for(j=0;j<K-1;j++){
        pc_tmp.col(j) = delta0+xi.col(j);
      }
      pc_tmp.col(K-1) = delta0;
      
      // log link
      for(j=0;j<K;j++){
        pc_nr.col(j).subvec(0,K-2) = exp(pc_tmp.col(j))/(1+sum(exp(pc_tmp.col(j))));
        pc_nr(K-1,j) = 1/(1+sum(exp(pc_tmp.col(j))));
      }
      
      dfsum.zeros();
      for(k=0;k<K-1;k++){
        for(j=k;j<K-1;j++){
          if(j==k){
            dfsum(k,j) = sum(pc_nr.row(k)%(1-pc_nr.row(k))%pm.row(t-1));
          }else{
            dfsum(k,j) = -sum(pc_nr.row(k)%pc_nr.row(j)%pm.row(t-1));
          }
        }
      }
      dfsum += trans(dfsum);
      dfsum.diag() = dfsum.diag()/2;
      
      f = pc_nr*trans(pm.row(t-1)) - trans(pm.row(t));
      delta1 = delta0 - solve(dfsum,f.subvec(0,K-2),solve_opts::fast)*stepNR;
      // delta1 = delta0 - pinv(dfsum)*f.subvec(0,K-2)*stepNR;
      diff_nr = sqrt(sum((delta1-delta0)%(delta1-delta0)));
      delta0 = delta1;
      iconv_nr = iconv_nr + 1;
    }
    pc.slice(t) = pc_nr;
    delta.row(t) = trans(delta0);
  }
}

// [[Rcpp::export]]
arma::cube CondProbCpp(arma::mat xi, int Ti, arma::mat& delta){
  
  K = xi.n_cols+1;
  
  vec yvec(K, fill::zeros);
  arma::mat pc_tmp(K-1, K, fill::zeros);
  arma::cube pc(K, K, Ti, fill::zeros);
  
  for(t=1;t<Ti;t++){
    pc_tmp.zeros();
    for(k=0;k<K-1;k++){
      for(j=0;j<K;j++){
        // indicator vector in conditional model
        yvec.zeros();
        yvec(j) = 1;
        pc_tmp(k,j) = delta(t,k) + \
          sum(trans(xi.row(k))%yvec.subvec(0,K-2));
      }
    }
    
    
    for(j=0;j<K;j++){
      for(k=0;k<K-1;k++){
        pc(k,j,t) = exp(pc_tmp(k,j))/(1+sum(exp(pc_tmp.col(j))));
      }
      pc(K-1,j,t) = 1/(1+sum(exp(pc_tmp.col(j))));
    }
  }
  return(pc);
}

// [[Rcpp::export]]
double LoglikCalcCpp(arma::vec par, arma::mat& DesignX, arma::mat& y_log, arma::vec& mm, double conv, double stepNR){
  
  // extract parameters
  Ti        = y_log.n_rows;
  K         = y_log.n_cols;
  len_mean  = DesignX.n_cols; 
  len_par   = par.size();
  
  // double conv = 1e-05;
  
  //-----------  Allocation Memory  -----------//
  arma::mat theta(K-1, len_mean, fill::zeros);
  arma::mat xi(K-1, K-1, fill::zeros);
  arma::mat eta(Ti, K-1, fill::zeros);
  arma::mat gamm(Ti, K-1, fill::zeros);
  arma::mat pm(Ti, K, fill::zeros);
  arma::mat delta(Ti, K-1, fill::zeros);
  
  arma::cube pc(K, K, Ti, fill::zeros);
  
  for(i=0;i<len_mean;i++){
    theta.col(i) = par.subvec(i*(K-1),(i+1)*(K-1)-1);
  }
  for(i=0;i<K-1;i++){
    xi.col(i) = par.subvec(len_mean*(K-1)+i*(K-1), \
           len_mean*(K-1)+i*(K-1)+(K-1)-1);
  }
  
  //----------- Marginal Parameter -----------//
  eta = DesignX*trans(theta);
  gamm = exp(eta)/(1 + exp(eta));
  
  // pm.zeros();
  for(i=0;i<K-1;i++){
    if(i==0){
      pm.col(i) = gamm.col(i);
    }else{
      pm.col(i) = gamm.col(i) - gamm.col(i-1);
    }
  }
  pm.col(K-1) = 1- gamm.col(K-2);
  
  // Check marginal probability positive
  while(bool positive_prob=any(vectorise(pm)<0)) {
    Rcout<<"Marginal Probability less than 0 !!!"<<std::endl;
    break;
  }
  
  //----- Newton-Raphson: Calculation of delta_tk values -----//
  //-----    Re-Calculation on Conditional Probability   -----//
  NRDeltaCpp(xi, pm, DesignX, delta, pc, K, conv, stepNR);
  
  //--------------------   Loglikelihood  --------------------//
  double tmpsum = 0.0;
  for(k=1;k<K-1;k++){
    tmpsum += y_log(0,k)*log(gamm(0,k)-gamm(0,k-1));
  }
  double logL1 = y_log(0,0)*log(gamm(0,0)) + tmpsum + \
    y_log(0,K-1)*log(1-gamm(0,K-2));
  
  double logL2 = 0.0;
  for(t=1;t<Ti;t++){
    tmpsum = 0;
    for(k=0;k<K-1;k++){
      tmpsum += y_log(t,k)*log(pc(k,mm(t-1),t)/pc(K-1,mm(t-1),t));
    }
    logL2 += tmpsum + log(pc(K-1,mm(t-1),t));
  }
  double logL = logL1 + logL2;
  
  return(logL);
}

// [[Rcpp::export]]
List GradHessCpp(arma::vec par, arma::mat DesignX, arma::mat y_log, arma::vec m, double conv, double stepNR){
  
  // extract parameters
  Ti        = y_log.n_rows;
  K         = y_log.n_cols;
  len_mean  = DesignX.n_cols; 
  len_par   = par.size();
  len_Marg  = len_mean*(K-1);
  len_Dep   = pow(K-1,2.0);
  
  // double conv = 1e-05;
  
  //-----------  Allocation Memory  -----------//
  double tmpsum1, tmpsum2;
  arma::mat theta(K-1, len_mean, fill::zeros);
  arma::mat xi(K-1, K-1, fill::zeros);
  arma::mat eta(Ti, K-1, fill::zeros);
  arma::mat gamm(Ti, K-1, fill::zeros);
  arma::mat pm(Ti, K, fill::zeros);
  arma::mat delta(Ti, K-1, fill::zeros);
  arma::mat H(K-1, K-1, fill::zeros);
  arma::mat s1(len_mean*(K-1), K-1, fill::zeros);
  arma::mat s2(K-1, len_Dep, fill::zeros);
  
  arma::cube pc(K, K, Ti, fill::zeros);
  
  for(i=0;i<len_mean;i++){
    theta.col(i) = par.subvec(i*(K-1),(i+1)*(K-1)-1);
  }
  for(i=0;i<K-1;i++){
    xi.col(i) = par.subvec(len_mean*(K-1)+i*(K-1), \
           len_mean*(K-1)+i*(K-1)+(K-1)-1);
  }
  
  Timer timer;
  timer.step("start"); 
  
  //----------- Marginal Parameter -----------//
  eta = DesignX*trans(theta);
  gamm = exp(eta)/(1 + exp(eta));
  pm.col(0) = gamm.col(0);
  pm.submat(0,1,Ti-1,K-2) = gamm.submat(0,1,Ti-1,K-2) - \
    gamm.submat(0,0,Ti-1,K-3);
  pm.col(K-1) = 1- gamm.col(K-2);
  
  // Check marginal probability positive
  while(bool positive_prob=any(vectorise(pm)<0)) {
    Rcout<<"Marginal Probability less than 0 !!!"<<std::endl;
    break;
  }
  
  timer.step("Marginal Prob"); 
  
  //----- Newton-Raphson: Calculation of delta_tk values -----//
  NRDeltaCpp(xi, pm, DesignX, delta, pc, K, conv, stepNR);
  
  timer.step("NR & Conditional Prob"); 
  
  //-----    Re-Calculation on Conditional Probability   -----//
  // pc = CondProbCpp(xi, Ti, delta);
  
  //---------- Derivatives of eta_tk w.r.t theta -------------//
  // (index of theta(i), eta_tk(k), eta_tk(t))
  arma::cube deta_dtheta(len_mean*(K-1), K-1, Ti, fill::zeros);   
  for(t=0;t<Ti;t++){
    for(k=0;k<K-1;k++){
      for(i=0;i<len_mean;i++){
        deta_dtheta(i*(K-1)+k,k,t) = DesignX(t,i);
      }
    }
  }
  
  //--------- Derivatives of gamm_tk w.r.t theta -------------//
  // (index of theta(kp), gamma_tk(k), gamma_tk(tt))
  arma::mat invgamm = 1-gamm;
  arma::cube dgamm_dtheta(len_Marg, K-1, Ti, fill::zeros);  
  for(t=0;t<Ti;t++){
    for(k=0;k<K-1;k++){
      dgamm_dtheta.slice(t).col(k) = gamm(t,k)*invgamm(t,k)*\
        deta_dtheta.slice(t).col(k);
    }
  }
  
  //--------- Derivatives of pm_tk w.r.t theta --------------//
  //index of theta(kp), pm_tk(k), pm_tk(tt)
  arma::cube dpm_dtheta(len_Marg, K, Ti, fill::zeros);
  // k = 1
  dpm_dtheta.subcube(0,0,0, len_Marg-1,0,Ti-1) = \
    dgamm_dtheta.subcube(0,0,0, len_Marg-1,0,Ti-1);
  // k = 2:(K-1)
  dpm_dtheta.subcube(0,1,0, len_Marg-1,K-2,Ti-1) =     \
    dgamm_dtheta.subcube(0,1,0, len_Marg-1,K-2,Ti-1) - \
    dgamm_dtheta.subcube(0,0,0, len_Marg-1,K-3,Ti-1);
  // k = K
  dpm_dtheta.subcube(0,K-1,0, len_Marg-1,K-1,Ti-1) = \
    -dgamm_dtheta.subcube(0,K-2,0, len_Marg-1,K-2,Ti-1);
    
    //--------- Derivatives of pc_tjk w.r.t delta -------------//
    // delta_tk(kp), pc_tkj(k), pc_tkj(j)*pc_tkj(tt)
    arma::cube dpc_ddelta(K-1, K, K*Ti, fill::zeros);
    for(t=0;t<Ti;t++){
      for(j=0;j<K;j++){
        for(k=0;k<K;k++){
          for(kp=0;kp<K-1;kp++){
            if(k==kp){
              dpc_ddelta(kp,k,j*Ti+t) = pc(k,j,t)*(1-pc(k,j,t));
            }else{
              dpc_ddelta(kp,k,j*Ti+t) = -pc(k,j,t)*pc(kp,j,t);
            }
          }
        }
      }
    }
    
    //--------- Derivatives of pc_tjk w.r.t xi   -------------//
    // xi(a,b), pc_tkj(k), pc_tkj(j)*pc_tkj(tt)
    arma::cube dpc_dxi(len_Dep, K, K*Ti, fill::zeros);
    for(t=0;t<Ti;t++){
      for(a=0;a<K-1;a++){
        for(b=0;b<K-1;b++){
          for(j=0;j<K;j++){
            if(b==j){
              for(k=0;k<K-1;k++){
                if(k==a){
                  dpc_dxi(b*(K-1)+a,k,j*Ti+t) = pc(k,j,t)*(1-pc(k,j,t));
                }else{
                  dpc_dxi(b*(K-1)+a,k,j*Ti+t) = -pc(k,j,t)*pc(a,j,t);
                }
              }
            }
            // k = K
            dpc_dxi(b*(K-1)+a,K-1,j*Ti+t) = -pc(K-1,j,t)*pc(a,j,t);
          }
        }
      }
    }
    
    timer.step("Basic Derivative"); 
    
    //----------     System Equation Solve     ---------------//
    
    //-------  Derivatives of Delta_tk w.r.t theta  ----------//
    //theta(kp), delta_tk(k), delta_tk(tt)
    arma::cube ddelta_dtheta(len_mean*(K-1), K-1, Ti, fill::zeros);
    //--------  Derivatives of Delta_tk w.r.t xi   -----------//
    //xi(a,b), delta_tk(k), delta_tk(tt)
    arma::cube ddelta_dxi(len_Dep, K-1, Ti, fill::zeros);
    
    for(t=1;t<Ti;t++){
      // Left Hand Side of eq1 & eq2
      H.zeros();
      for(j=0;j<K;j++){
        H += trans(dpc_ddelta.slice(j*Ti+t).submat(0,0,K-2,K-2)*pm(t-1,j));
      }
      // Right Hand Side of eq1 (s1: theta(kp), pm_tk(k))
      s1 = dpm_dtheta.slice(t).submat(0,0,len_Marg-1,K-2)-\
        trans(pc.slice(t).submat(0,0,K-2,K-1)*            \
        trans(dpm_dtheta.slice(t-1)));
      ddelta_dtheta.slice(t) = trans(solve(H,trans(s1),solve_opts::fast));
      // Left Hand Side of eq2 (s2: pc_tjk(k) & xi_kj(a,b))
      s2.zeros();
      for(j=0;j<K;j++){
        s2 -= trans(dpc_dxi.slice(j*Ti+t).submat(0,0,len_Dep-1,K-2)*pm(t-1,j));
      }
      ddelta_dxi.slice(t) = trans(solve(H,s2,solve_opts::fast));
    }
    
    timer.step("System Equation Solve"); 
    //+++++++++++++++++++++++++++++++++++++++++++//
    //               Gradient                    //
    //+++++++++++++++++++++++++++++++++++++++++++//
    
    //------------- theta -------------//
    // dlogL1_dtheta: len_mean*(K-1) * 1
    arma::vec tmpsumdL1(len_mean*(K-1), fill::zeros);
    arma::vec dlogL1_dtheta(len_mean*(K-1), fill::zeros);
    for(k=1;k<K-1;k++){
      tmpsumdL1 += y_log(0,k)/(gamm(0,k)-gamm(0,k-1))*\
        (dgamm_dtheta.slice(0).col(k)-dgamm_dtheta.slice(0).col(k-1));
    }
    dlogL1_dtheta = y_log(0,0)/gamm(0,0)*dgamm_dtheta.slice(0).col(0) + \
      tmpsumdL1 + y_log(0,K-1)/(1-gamm(0,K-2))*                         \
      (0-dgamm_dtheta.slice(0).col(K-2));
    
    // dlogL2_dtheta: len_mean*(K-1) * 1
    arma::vec dlogL2_dtheta(len_mean*(K-1), fill::zeros);
    for(t=1;t<Ti;t++){
      dlogL2_dtheta += ddelta_dtheta.slice(t)*\
        (trans(y_log.row(t).subvec(0,K-2))-   \
        pc.slice(t).col(m(t-1)).subvec(0,K-2));
    }
    
    //-------------  xi  -------------//
    // dlogL1_dxi: (K-1)*(K-1) * 1 (have to be zero)
    // dlogL2_dxi: (K-1)*(K-1) * 1
    arma::vec dlogL2_dxi(len_Dep, fill::zeros);
    for(a=0;a<K-1;a++){
      for(b=0;b<K-1;b++){
        for(t=1;t<Ti;t++){
          dlogL2_dxi(b*(K-1)+a) +=                        \
            sum((y_log.row(t).subvec(0,K-2)-              \
            trans(pc.slice(t).col(m(t-1)).subvec(0,K-2)))%\
            ddelta_dxi.slice(t).row(b*(K-1)+a))+          \
            (y_log(t,a)-pc(a,m(t-1),t))*y_log(t-1,b);
        }
      }
    }
    
    arma::vec dlogL_dphi(len_mean*(K-1)+len_Dep,fill::zeros);
    dlogL_dphi.subvec(0,len_mean*(K-1)-1) = dlogL1_dtheta + dlogL2_dtheta;
    dlogL_dphi.subvec(len_mean*(K-1),len_mean*(K-1)+(K-1)*(K-1)-1)=dlogL2_dxi;
    
    timer.step("Gradient"); 
    
    //+++++++++++++++++++++++++++++++++++++++++++//
    //                Hessian                    //
    //+++++++++++++++++++++++++++++++++++++++++++// 
    arma::mat ddlogL_ddphi(len_par, len_par, fill::zeros);
    
    // ddgamm_dthetatheta: (len_mean*(K-1)) * (len_mean*(K-1)) * (K-1)
    // index of , gamm_tk(k), theta_star(p), theta (q)
    arma::cube ddgamm_dthetatheta(K-1, len_mean*(K-1), \
                                  len_mean*(K-1), fill::zeros);
    for(p=0;p<len_mean*(K-1);p++){
      for(q=0;q<len_mean*(K-1);q++){
        ddgamm_dthetatheta.slice(q).col(p) = trans((1-2*gamm.row(0))%  \
          gamm.row(0)%(1-gamm.row(0))%deta_dtheta.slice(0).row(p)%     \
          deta_dtheta.slice(0).row(q));
      }
    }
    
    
    //-------------  theta by theta  -------------//
    // ddlogL1_ddtheta: len_mean*(K-1) * len_mean*(K-1)
    arma::mat ddlogL1_ddtheta(len_mean*(K-1), len_mean*(K-1), fill::zeros);
    
    for(p=0;p<len_mean*(K-1);p++){
      for(q=0;q<len_mean*(K-1);q++){
        tmpsum1=0.0;
        for(k=1;k<K-1;k++){
          tmpsum1 += y_log(0,k)*(1/pow((gamm(0,k)-gamm(0,k-1)),2.0)*\
            (dgamm_dtheta(p,k,0)-dgamm_dtheta(p,k-1,0))*            \
            (dgamm_dtheta(q,k,0)-dgamm_dtheta(q,k-1,0)) -           \
            1/(gamm(0,k)-gamm(0,k-1))*                              \
            (ddgamm_dthetatheta(k,p,q)-ddgamm_dthetatheta(k-1,p,q)));
        }
        ddlogL1_ddtheta(p,q) = y_log(0,0)*(1/(pow(gamm(0,0),2.0))*\
          dgamm_dtheta(p,0,0)*                                    \
          dgamm_dtheta(q,0,0)-                                    \
          1/gamm(0,0)*                                            \
          ddgamm_dthetatheta(0,p,q)) +                            \
          tmpsum1 +                                               \
          y_log(0,K-1)*(1/(pow(gamm(0,K-1-1),2.0))*               \
          dgamm_dtheta(p,0,K-1-1)*                                \
          dgamm_dtheta(q,0,K-1-1)+                                \
          1/gamm(0,K-1-1)*                                        \
          ddgamm_dthetatheta(K-1-1,p,q));
      }
    }
    
    // ddlogL1_ddtheta: len_mean*(K-1) * len_mean*(K-1)
    arma::mat ddlogL2_ddtheta(len_mean*(K-1), len_mean*(K-1), fill::zeros);
    for(t=1;t<Ti;t++){
      ddlogL2_ddtheta += ddelta_dtheta.slice(t).submat(0,0,len_Marg-1,K-2)*  \
        trans(ddelta_dtheta.slice(t)*                                        \
        dpc_ddelta.slice(m(t-1)*Ti+t).submat(0,0,K-2,K-2));
    }
    timer.step("Hessian theta by theta"); 
    
    //-------------    xi by theta   -------------//
    arma::mat ddlogL2_dxidtheta(len_Dep, len_mean*(K-1), fill::zeros);
    // ddlogL2_dxidtheta: (K-1)*(K-1) * len_mean*(K-1)
    for(a=0;a<K-1;a++){
      for(b=0;b<K-1;b++){
        for(q=0;q<len_mean*(K-1);q++){
          for(t=1;t<Ti;t++){
            for(k=0;k<K-1;k++){
              tmpsum1 = 0.0;
              for(g=0;g<K-1;g++){
                tmpsum1 += dpc_ddelta(g,k,m(t-1)*Ti+t)*\
                  ddelta_dxi(b*(K-1)+a,g,t)*           \
                  ddelta_dtheta(q,k,t);
              }
              ddlogL2_dxidtheta(b*(K-1)+a,q) += tmpsum1 + \
                dpc_dxi(b*(K-1)+a,k,m(t-1)*Ti+t)*ddelta_dtheta(q,k,t);
            }
          }
        }
      }
    }
    timer.step("Hessian xi by theta");
    
    //-------------      xi by xi    -------------// 
    arma::mat ddlogL2_ddxi(len_Dep, len_Dep, fill::zeros);
    int tmpindex;
    for(t=1;t<Ti;t++){
      tmpindex = m(t-1)*Ti+t;
      for(p=0;p<len_Dep;p++){
        for(q=p;q<len_Dep;q++){
          a = q%(K-1); // get the remainder
          b = ceil(q/(K-1));
          
          ddlogL2_ddxi(q,p) += 
            sum(ddelta_dxi.slice(t).row(p)*
            dpc_ddelta.slice(tmpindex).submat(0,0,K-2,K-2)%  \
            ddelta_dxi.slice(t).row(q))+                     \
            sum(dpc_dxi.slice(tmpindex).row(p).subvec(0,K-2)%\
            ddelta_dxi.slice(t).row(q))+                     \
            (sum(ddelta_dxi.slice(t).row(p)*                 \
            dpc_ddelta.slice(tmpindex).col(a))+              \
            dpc_dxi(p,a,tmpindex))*                          \
            y_log(t-1,b);
        }
      }
    }
    
    ddlogL2_ddxi += trans(ddlogL2_ddxi);
    //Since double count the diagonal elements
    ddlogL2_ddxi.diag() = ddlogL2_ddxi.diag()/2;
    timer.step("Hessian xi by xi"); 
    
    // ddtheta_ddtheta
    ddlogL_ddphi.submat(0,0,len_mean*(K-1)-1,len_mean*(K-1)-1) = \
      -(ddlogL1_ddtheta + ddlogL2_ddtheta);
      // ddxi_ddtheta
      ddlogL_ddphi.submat(len_mean*(K-1),0,len_par-1,len_mean*(K-1)-1) = \
      ddlogL2_dxidtheta;
      // ddtheta_ddxi
      ddlogL_ddphi.submat(0,len_mean*(K-1),len_mean*(K-1)-1,len_par-1) = \
        trans(ddlogL2_dxidtheta);
      // ddxi_ddxi
      ddlogL_ddphi.submat(len_mean*(K-1),len_mean*(K-1),len_par-1,len_par-1) = \
        ddlogL2_ddxi;
      
      NumericVector rest(timer);   // 
      int n = 1000000;
      for (int i=0; i<rest.size(); i++) {
        rest[i] = rest[i] / n;
      }
      
      List res;
      res["myG"] = dlogL_dphi;
      res["myH"] = ddlogL_ddphi;
      res["tim"] = rest;
      return(res);
}

// [[Rcpp::export]]
void GradCppnewVec21(arma::vec& par, arma::mat& DesignX, arma::mat& delta, arma::mat& y_log, arma::vec& m, arma::vec& Gradtheta, arma::mat& Hesstheta, double conv, double stepNR){
  
  // extract parameters
  Ti        = y_log.n_rows;
  K         = y_log.n_cols;
  len_mean  = DesignX.n_cols; 
  len_par   = par.size();
  len_Marg  = len_mean*(K-1);
  len_Dep   = pow(K-1,2.0);
  
  // double conv = 1e-05;
  
  //-----------  Allocation Memory  -----------//
  double tmpsum1, tmpsum2, tmpsum3;
  arma::mat theta(K-1, len_mean, fill::zeros);
  arma::mat xi(K-1, K-1, fill::zeros);
  arma::mat eta(Ti, K-1, fill::zeros);
  arma::mat gamm(Ti, K-1, fill::zeros);
  arma::mat pm(Ti, K, fill::zeros);
  arma::mat H(K-1, K-1, fill::zeros);
  arma::mat s1(len_Marg, K-1, fill::zeros);
  
  arma::cube pc(K, K, Ti, fill::zeros);
  
  //-----------  Allocation Memory  -----------//
  for(i=0;i<len_mean;i++){
    theta.col(i) = par.subvec(i*(K-1),(i+1)*(K-1)-1);
  }
  for(i=0;i<K-1;i++){
    xi.col(i) = par.subvec(len_Marg+i*(K-1), \
           len_Marg+i*(K-1)+(K-1)-1);
  }
  
  //----------- Marginal Parameter -----------//
  eta = DesignX*trans(theta);
  gamm = exp(eta)/(1 + exp(eta));
  pm.col(0) = gamm.col(0);
  pm.submat(0,1,Ti-1,K-2) = gamm.submat(0,1,Ti-1,K-2) - \
    gamm.submat(0,0,Ti-1,K-3);
  pm.col(K-1) = 1- gamm.col(K-2);
  
  // Check marginal probability positive
  while(bool positive_prob=any(vectorise(pm)<0)) {
    Rcout<<"Marginal Probability less than 0 !!!"<<std::endl;
    break;
  }
  
  //----- Newton-Raphson: Calculation of delta_tk values -----//
  //------    Calculation on Conditional Probability   -------//
  NRDeltaCpp(xi, pm, DesignX, delta, pc, K, conv, stepNR);
  
  //---------- Derivatives of eta_tk w.r.t theta -------------//
  // (index of theta(i), eta_tk(k), eta_tk(t))
  arma::cube deta_dtheta(len_Marg, K-1, Ti, fill::zeros);   
  for(t=0;t<Ti;t++){
    for(k=0;k<K-1;k++){
      for(i=0;i<len_mean;i++){
        deta_dtheta(i*(K-1)+k,k,t) = DesignX(t,i);
      }
    }
  }
  
  //--------- Derivatives of gamm_tk w.r.t theta -------------//
  // (index of theta(kp), gamma_tk(k), gamma_tk(tt))
  arma::mat invgamm = 1-gamm;
  arma::cube dgamm_dtheta(len_Marg, K-1, Ti, fill::zeros);  
  for(t=0;t<Ti;t++){
    for(k=0;k<K-1;k++){
      dgamm_dtheta.slice(t).col(k) = gamm(t,k)*invgamm(t,k)*\
        deta_dtheta.slice(t).col(k);
    }
  }
  
  //--------- Derivatives of pm_tk w.r.t theta --------------//
  //index of theta(kp), pm_tk(k), pm_tk(tt)
  arma::cube dpm_dtheta(len_Marg, K, Ti, fill::zeros);
  // k = 1
  dpm_dtheta.subcube(0,0,0, len_Marg-1,0,Ti-1) = \
    dgamm_dtheta.subcube(0,0,0, len_Marg-1,0,Ti-1);
  // k = 2:(K-1)
  dpm_dtheta.subcube(0,1,0, len_Marg-1,K-2,Ti-1) =     \
    dgamm_dtheta.subcube(0,1,0, len_Marg-1,K-2,Ti-1) - \
    dgamm_dtheta.subcube(0,0,0, len_Marg-1,K-3,Ti-1);
  // k = K
  dpm_dtheta.subcube(0,K-1,0, len_Marg-1,K-1,Ti-1) = \
    -dgamm_dtheta.subcube(0,K-2,0, len_Marg-1,K-2,Ti-1);
    
    //--------- Derivatives of pc_tjk w.r.t delta -------------//
    // delta_tk(kp), pc_tkj(k), pc_tkj(j)*pc_tkj(tt)
    arma::cube dpc_ddelta(K-1, K, K*Ti, fill::zeros);
    for(t=0;t<Ti;t++){
      for(j=0;j<K;j++){
        for(k=0;k<K;k++){
          for(kp=0;kp<K-1;kp++){
            if(k==kp){
              dpc_ddelta(kp,k,j*Ti+t) = pc(k,j,t)*(1-pc(k,j,t));
            }else{
              dpc_ddelta(kp,k,j*Ti+t) = -pc(k,j,t)*pc(kp,j,t);
            }
          }
        }
      }
    }
    
    //----------     System Equation Solve     ---------------//
    
    //-------  Derivatives of Delta_tk w.r.t theta  ----------//
    //theta(kp), delta_tk(k), delta_tk(tt)
    arma::cube ddelta_dtheta(len_Marg, K-1, Ti, fill::zeros);
    
    for(t=1;t<Ti;t++){
      // Left Hand Side of eq1 & eq2
      H.zeros();
      for(j=0;j<K;j++){
        H += trans(dpc_ddelta.slice(j*Ti+t).submat(0,0,K-2,K-2)*pm(t-1,j));
      }
      
      // Right Hand Side of eq1 (s1: theta(kp), pm_tk(k))
      s1 = dpm_dtheta.slice(t).submat(0,0,len_Marg-1,K-2)-\
        trans(pc.slice(t).submat(0,0,K-2,K-1)*            \
        trans(dpm_dtheta.slice(t-1)));
      ddelta_dtheta.slice(t) = trans(solve(H,trans(s1),solve_opts::fast));
    }
    
    //+++++++++++++++++++++++++++++++++++++++++++//
    //               Gradient                    //
    //+++++++++++++++++++++++++++++++++++++++++++//
    
    //------------- theta -------------//
    arma::mat dlogL_dtheta(len_Marg, Ti, fill::zeros);
    // dlogL1_dtheta: len_Marg * 1
    arma::vec tmpsumdL1(len_Marg, fill::zeros);
    for(k=1;k<K-1;k++){
      tmpsumdL1 += y_log(0,k)/(gamm(0,k)-gamm(0,k-1))*\
        (dgamm_dtheta.slice(0).col(k)-dgamm_dtheta.slice(0).col(k-1));
    }
    dlogL_dtheta.col(0) = y_log(0,0)/gamm(0,0)*dgamm_dtheta.slice(0).col(0) +  \
      tmpsumdL1 + y_log(0,K-1)/(1-gamm(0,K-2))*                                \
      (0-dgamm_dtheta.slice(0).col(K-2));
    // dlogL2_dtheta: len_Marg * 1
    for(t=1;t<Ti;t++){
      dlogL_dtheta.col(t) = ddelta_dtheta.slice(t)*\
        (trans(y_log.row(t).subvec(0,K-2))-        \
        pc.slice(t).col(m(t-1)).subvec(0,K-2));
    }
    Gradtheta = sum(dlogL_dtheta,1);
    
    //+++++++++++++++++++++++++++++++++++++++++++//
    //         Hessian Approximation             //
    //+++++++++++++++++++++++++++++++++++++++++++//
    Hesstheta = dlogL_dtheta*trans(dlogL_dtheta);
    
}

// [[Rcpp::export]]
void GradCppnewVec22(arma::vec& par, arma::mat& DesignX, arma::mat& delta, arma::mat& y_log, arma::vec& m, arma::vec& Gradxi, arma::mat& Hessxi, double conv, double stepNR){
  
  // extract parameters
  Ti        = y_log.n_rows;
  K         = y_log.n_cols;
  len_mean  = DesignX.n_cols;
  len_par   = par.size();
  len_Marg  = len_mean*(K-1);
  len_Dep   = pow(K-1,2.0);
  
  // double conv = 1e-05;
  
  //-----------  Allocation Memory  -----------//
  double tmpsum1, tmpsum2, tmpsum3;
  arma::mat theta(K-1, len_mean, fill::zeros);
  arma::mat xi(K-1, K-1, fill::zeros);
  arma::mat eta(Ti, K-1, fill::zeros);
  arma::mat gamm(Ti, K-1, fill::zeros);
  arma::mat pm(Ti, K, fill::zeros);
  arma::mat H(K-1, K-1, fill::zeros);
  arma::mat s1(len_Marg, K-1, fill::zeros);
  arma::mat s2(K-1, len_Dep, fill::zeros);
  
  arma::cube pc(K, K, Ti, fill::zeros);
  
  //-----------  Allocation Memory  -----------//
  for(i=0;i<len_mean;i++){
    theta.col(i) = par.subvec(i*(K-1),(i+1)*(K-1)-1);
  }
  for(i=0;i<K-1;i++){
    xi.col(i) = par.subvec(len_Marg+i*(K-1), \
           len_Marg+i*(K-1)+(K-1)-1);
  }
  
  //----------- Marginal Parameter -----------//
  eta = DesignX*trans(theta);
  gamm = exp(eta)/(1 + exp(eta));
  pm.col(0) = gamm.col(0);
  pm.submat(0,1,Ti-1,K-2) = gamm.submat(0,1,Ti-1,K-2) - \
    gamm.submat(0,0,Ti-1,K-3);
  pm.col(K-1) = 1- gamm.col(K-2);
  
  // Check marginal probability positive
  while(bool positive_prob=any(vectorise(pm)<0)) {
    Rcout<<"Marginal Probability less than 0 !!!"<<std::endl;
    break;
  }
  
  //----- Newton-Raphson: Calculation of delta_tk values -----//
  //-----    Re-Calculation on Conditional Probability   -----//
  NRDeltaCpp(xi, pm, DesignX, delta, pc, K, conv, stepNR);
  
  //---------- Derivatives of eta_tk w.r.t theta -------------//
  // (index of theta(i), eta_tk(k), eta_tk(t))
  arma::cube deta_dtheta(len_Marg, K-1, Ti, fill::zeros);
  for(t=0;t<Ti;t++){
    for(k=0;k<K-1;k++){
      for(i=0;i<len_mean;i++){
        deta_dtheta(i*(K-1)+k,k,t) = DesignX(t,i);
      }
    }
  }
  
  //--------- Derivatives of gamm_tk w.r.t theta -------------//
  // (index of theta(kp), gamma_tk(k), gamma_tk(tt))
  arma::mat invgamm = 1-gamm;
  arma::cube dgamm_dtheta(len_Marg, K-1, Ti, fill::zeros);
  for(t=0;t<Ti;t++){
    for(k=0;k<K-1;k++){
      dgamm_dtheta.slice(t).col(k) = gamm(t,k)*invgamm(t,k)*\
        deta_dtheta.slice(t).col(k);
    }
  }
  
  //--------- Derivatives of pm_tk w.r.t theta --------------//
  //index of theta(kp), pm_tk(k), pm_tk(tt)
  arma::cube dpm_dtheta(len_Marg, K, Ti, fill::zeros);
  // k = 1
  dpm_dtheta.subcube(0,0,0, len_Marg-1,0,Ti-1) = \
    dgamm_dtheta.subcube(0,0,0, len_Marg-1,0,Ti-1);
  // k = 2:(K-1)
  dpm_dtheta.subcube(0,1,0, len_Marg-1,K-2,Ti-1) =     \
    dgamm_dtheta.subcube(0,1,0, len_Marg-1,K-2,Ti-1) - \
    dgamm_dtheta.subcube(0,0,0, len_Marg-1,K-3,Ti-1);
  // k = K
  dpm_dtheta.subcube(0,K-1,0, len_Marg-1,K-1,Ti-1) = \
    -dgamm_dtheta.subcube(0,K-2,0, len_Marg-1,K-2,Ti-1);
    
    //--------- Derivatives of pc_tjk w.r.t delta -------------//
    // delta_tk(kp), pc_tkj(k), pc_tkj(j)*pc_tkj(tt)
    arma::cube dpc_ddelta(K-1, K, K*Ti, fill::zeros);
    for(t=0;t<Ti;t++){
      for(j=0;j<K;j++){
        for(k=0;k<K;k++){
          for(kp=0;kp<K-1;kp++){
            if(k==kp){
              dpc_ddelta(kp,k,j*Ti+t) = pc(k,j,t)*(1-pc(k,j,t));
            }else{
              dpc_ddelta(kp,k,j*Ti+t) = -pc(k,j,t)*pc(kp,j,t);
            }
          }
        }
      }
    }
    
    //--------- Derivatives of pc_tjk w.r.t xi   -------------//
    // xi(a,b), pc_tkj(k), pc_tkj(j)*pc_tkj(tt)
    arma::cube dpc_dxi(len_Dep, K, K*Ti, fill::zeros);
    for(t=0;t<Ti;t++){
      for(a=0;a<K-1;a++){
        for(b=0;b<K-1;b++){
          for(j=0;j<K;j++){
            if(b==j){
              for(k=0;k<K-1;k++){
                if(k==a){
                  dpc_dxi(b*(K-1)+a,k,j*Ti+t) = pc(k,j,t)*(1-pc(k,j,t));
                }else{
                  dpc_dxi(b*(K-1)+a,k,j*Ti+t) = -pc(k,j,t)*pc(a,j,t);
                }
              }
            }
            // k = K
            dpc_dxi(b*(K-1)+a,K-1,j*Ti+t) = -pc(K-1,j,t)*pc(a,j,t);
          }
        }
      }
    }
    
    //----------     System Equation Solve     ---------------//
    
    //-------  Derivatives of Delta_tk w.r.t theta  ----------//
    //theta(kp), delta_tk(k), delta_tk(tt)
    arma::cube ddelta_dtheta(len_Marg, K-1, Ti, fill::zeros);
    //--------  Derivatives of Delta_tk w.r.t xi   -----------//
    //xi(a,b), delta_tk(k), delta_tk(tt)
    arma::cube ddelta_dxi(len_Dep, K-1, Ti, fill::zeros);
    
    for(t=1;t<Ti;t++){
      // Left Hand Side of eq1 & eq2
      H.zeros();
      for(j=0;j<K;j++){
        H += trans(dpc_ddelta.slice(j*Ti+t).submat(0,0,K-2,K-2)*pm(t-1,j));
      }
      // Right Hand Side of eq1 (s1: theta(kp), pm_tk(k))
      s1 = dpm_dtheta.slice(t).submat(0,0,len_Marg-1,K-2)-\
        trans(pc.slice(t).submat(0,0,K-2,K-1)*            \
        trans(dpm_dtheta.slice(t-1)));
      ddelta_dtheta.slice(t) = trans(solve(H,trans(s1),solve_opts::fast));
      
      // Left Hand Side of eq2 (s2: pc_tjk(k) & xi_kj(a,b))
      s2.zeros();
      for(j=0;j<K;j++){
        s2 -= trans(dpc_dxi.slice(j*Ti+t).submat(0,0,len_Dep-1,K-2)*pm(t-1,j));
      }
      ddelta_dxi.slice(t) = trans(solve(H,s2,solve_opts::fast));
    }
    
    //+++++++++++++++++++++++++++++++++++++++++++//
    //               Gradient                    //
    //+++++++++++++++++++++++++++++++++++++++++++//
    
    arma::mat dlogL_dxi(len_Dep, Ti, fill::zeros);
    //-------------  xi  -------------//
    // dlogL1_dxi: (K-1)*(K-1) * 1 (have to be zero)
    // dlogL2_dxi: (K-1)*(K-1) * 1
    for(a=0;a<K-1;a++){
      for(b=0;b<K-1;b++){
        for(t=1;t<Ti;t++){
          dlogL_dxi(b*(K-1)+a,t) =                          \
            sum((y_log.row(t).subvec(0,K-2)-                \
            trans(pc.slice(t).col(m(t-1)).subvec(0,K-2)))%  \
            ddelta_dxi.slice(t).row(b*(K-1)+a))+            \
            (y_log(t,a)-pc(a,m(t-1),t))*y_log(t-1,b);
        }
      }
    }
    Gradxi = sum(dlogL_dxi,1);
    
    //+++++++++++++++++++++++++++++++++++++++++++//
    //         Hessian Approximation             //
    //+++++++++++++++++++++++++++++++++++++++++++//
    Hessxi = dlogL_dxi*trans(dlogL_dxi);
    
}

// [[Rcpp::export]]
List SeqQuasiNWCpp(arma::vec par, arma::mat& y_log, arma::vec& mm, arma::mat& DesignX, arma::vec stepsize, int maxit, double mytol){
  
  //--------------------------- extract parameters ----------------------//
  Ti        = y_log.n_rows;
  K         = y_log.n_cols;
  len_mean  = DesignX.n_cols;
  len_par   = par.size();
  len_Marg  = (K-1)*len_mean;
  len_Dep   = (K-1)*(K-1);
  
  //--------------------------- Declar identifiers ----------------------//
  arma::vec phi0 = par;
  arma::vec phi1 = par;
  arma::vec phi0_out = par;
  arma::vec phi1_out = par;
  
  int iconv_NW  = 0;
  int iconv_NW1 = 0;
  int iconv_NW2 = 0;
  int convergence;
  
  double dif_NW  = 0.01;
  double dif_NW1 = 0.01;
  double dif_NW2 = 0.01;
  
  double stepNW = stepsize(0);
  double stepNR = stepsize(1);
  
  arma::uvec upId1 = arma::linspace<arma::uvec>(0, len_Marg-1, len_Marg);
  arma::uvec upId2 = arma::linspace<arma::uvec>(len_Marg, len_par-1, len_Dep);
  
  arma::mat delta(Ti,K-1,fill::zeros);
  
  // Gradient
  arma::vec Gradtheta(len_Marg, fill::zeros);
  arma::vec Gradxi(len_Dep, fill::zeros);
  
  // Hessian
  arma::mat Hesstheta(len_Marg, len_Marg, fill::zeros);
  arma::mat Hessxi(len_Dep, len_Dep, fill::zeros);
  
  //-------------------- Sequential Newton's Methods -------------------//
  while((iconv_NW<maxit) & (dif_NW>mytol)){
    
    // Updating Mean Parameter
    iconv_NW1 = 0; 
    dif_NW1   = 0.01;
    while((iconv_NW1<maxit) & (dif_NW1>mytol)){
      Gradtheta.zeros();
      Hesstheta.zeros();
      GradCppnewVec21(phi0, DesignX, delta, y_log, mm, Gradtheta, Hesstheta, mytol, stepNR);
      phi1(upId1) = phi0(upId1) + inv(Hesstheta)*Gradtheta*stepNW;
      dif_NW1 = pow(sum(square(phi1-phi0)), 0.5);
      iconv_NW1 +=1;
      phi0 = phi1;
    }
    // Check Convergence
    phi1_out = phi0;
    iconv_NW += iconv_NW1;
    dif_NW = pow(sum(square(phi0_out - phi1_out)), 0.5);
    phi0_out = phi1_out;
    if((iconv_NW>maxit) | (dif_NW<mytol)){
      break;
    }
    
    // Updating Dependence Parameter
    iconv_NW2 = 0; 
    dif_NW2   = 0.01;
    while((iconv_NW2<maxit) & (dif_NW2>mytol)){
      Gradxi.zeros();
      Hessxi.zeros();
      GradCppnewVec22(phi0, DesignX, delta, y_log, mm, Gradxi, Hessxi, mytol, stepNR);
      phi1(upId2) = phi0(upId2) + inv(Hessxi)*Gradxi*stepNW;
      dif_NW2 = pow(sum(square(phi1-phi0)), 0.5);
      iconv_NW2 +=1; 
      phi0 = phi1;
    }
    // Check Convergence
    phi1_out = phi0;
    iconv_NW += iconv_NW2;
    dif_NW = pow(sum(square(phi0_out - phi1_out)), 0.5);
    phi0_out = phi1_out;
  }
  
  //--------------------------- Return Result -------------------------//
  List res;
  res["parEst"] = phi0_out;
  res["counts"] = iconv_NW;
  res["Gradtheta"] = Gradtheta;
  res["Gradxi"] = Gradxi;
  res["Hesstheta"] = Hesstheta;
  res["Hessxi"] = Hessxi;
  if(iconv_NW < maxit){
    res["convergence"] = 0;
  }else{
    res["convergence"] = 1;
  }
  
  return(res);
}

// [[Rcpp::export]]
void GradHessCppnewVec21(arma::vec& par, arma::mat& DesignX, arma::mat& delta, arma::mat& y_log, arma::vec& m, arma::vec& dlogL_dtheta, arma::mat& ddlogL_ddtheta, double conv, double stepNR){
  
  // extract parameters
  Ti        = y_log.n_rows;
  K         = y_log.n_cols;
  len_mean  = DesignX.n_cols; 
  len_par   = par.size();
  len_Marg  = len_mean*(K-1);
  len_Dep   = pow(K-1,2.0);
  
  // double conv = 1e-05;
  
  //-----------  Allocation Memory  -----------//
  double tmpsum1, tmpsum2, tmpsum3;
  arma::mat theta(K-1, len_mean, fill::zeros);
  arma::mat xi(K-1, K-1, fill::zeros);
  arma::mat eta(Ti, K-1, fill::zeros);
  arma::mat gamm(Ti, K-1, fill::zeros);
  arma::mat pm(Ti, K, fill::zeros);
  arma::mat H(K-1, K-1, fill::zeros);
  arma::mat s1(len_Marg, K-1, fill::zeros);
  
  arma::cube pc(K, K, Ti, fill::zeros);
  
  //-----------  Allocation Memory  -----------//
  for(i=0;i<len_mean;i++){
    theta.col(i) = par.subvec(i*(K-1),(i+1)*(K-1)-1);
  }
  for(i=0;i<K-1;i++){
    xi.col(i) = par.subvec(len_Marg+i*(K-1), \
           len_Marg+i*(K-1)+(K-1)-1);
  }
  
  //----------- Marginal Parameter -----------//
  eta = DesignX*trans(theta);
  gamm = exp(eta)/(1 + exp(eta));
  pm.col(0) = gamm.col(0);
  pm.submat(0,1,Ti-1,K-2) = gamm.submat(0,1,Ti-1,K-2) - \
    gamm.submat(0,0,Ti-1,K-3);
  pm.col(K-1) = 1- gamm.col(K-2);
  
  // Check marginal probability positive
  while(bool positive_prob=any(vectorise(pm)<0)) {
    Rcout<<"Marginal Probability less than 0 !!!"<<std::endl;
    break;
  }
  
  //----- Newton-Raphson: Calculation of delta_tk values -----//
  //------    Calculation on Conditional Probability   -------//
  NRDeltaCpp(xi, pm, DesignX, delta, pc, K, conv, stepNR);
  
  //---------- Derivatives of eta_tk w.r.t theta -------------//
  // (index of theta(i), eta_tk(k), eta_tk(t))
  arma::cube deta_dtheta(len_Marg, K-1, Ti, fill::zeros);   
  for(t=0;t<Ti;t++){
    for(k=0;k<K-1;k++){
      for(i=0;i<len_mean;i++){
        deta_dtheta(i*(K-1)+k,k,t) = DesignX(t,i);
      }
    }
  }
  
  //--------- Derivatives of gamm_tk w.r.t theta -------------//
  // (index of theta(kp), gamma_tk(k), gamma_tk(tt))
  arma::mat invgamm = 1-gamm;
  arma::cube dgamm_dtheta(len_Marg, K-1, Ti, fill::zeros);  
  for(t=0;t<Ti;t++){
    for(k=0;k<K-1;k++){
      dgamm_dtheta.slice(t).col(k) = gamm(t,k)*invgamm(t,k)*\
        deta_dtheta.slice(t).col(k);
    }
  }
  
  //--------- Derivatives of pm_tk w.r.t theta --------------//
  //index of theta(kp), pm_tk(k), pm_tk(tt)
  arma::cube dpm_dtheta(len_Marg, K, Ti, fill::zeros);
  // k = 1
  dpm_dtheta.subcube(0,0,0, len_Marg-1,0,Ti-1) = \
    dgamm_dtheta.subcube(0,0,0, len_Marg-1,0,Ti-1);
  // k = 2:(K-1)
  dpm_dtheta.subcube(0,1,0, len_Marg-1,K-2,Ti-1) =     \
    dgamm_dtheta.subcube(0,1,0, len_Marg-1,K-2,Ti-1) - \
    dgamm_dtheta.subcube(0,0,0, len_Marg-1,K-3,Ti-1);
  // k = K
  dpm_dtheta.subcube(0,K-1,0, len_Marg-1,K-1,Ti-1) = \
    -dgamm_dtheta.subcube(0,K-2,0, len_Marg-1,K-2,Ti-1);
    
    //--------- Derivatives of pc_tjk w.r.t delta -------------//
    // delta_tk(kp), pc_tkj(k), pc_tkj(j)*pc_tkj(tt)
    arma::cube dpc_ddelta(K-1, K, K*Ti, fill::zeros);
    for(t=0;t<Ti;t++){
      for(j=0;j<K;j++){
        for(k=0;k<K;k++){
          for(kp=0;kp<K-1;kp++){
            if(k==kp){
              dpc_ddelta(kp,k,j*Ti+t) = pc(k,j,t)*(1-pc(k,j,t));
            }else{
              dpc_ddelta(kp,k,j*Ti+t) = -pc(k,j,t)*pc(kp,j,t);
            }
          }
        }
      }
    }
    
    //----------     System Equation Solve     ---------------//
    
    //-------  Derivatives of Delta_tk w.r.t theta  ----------//
    //theta(kp), delta_tk(k), delta_tk(tt)
    arma::cube ddelta_dtheta(len_Marg, K-1, Ti, fill::zeros);
    
    for(t=1;t<Ti;t++){
      // Left Hand Side of eq1 & eq2
      H.zeros();
      for(j=0;j<K;j++){
        H += trans(dpc_ddelta.slice(j*Ti+t).submat(0,0,K-2,K-2)*pm(t-1,j));
      }
      
      // Right Hand Side of eq1 (s1: theta(kp), pm_tk(k))
      s1 = dpm_dtheta.slice(t).submat(0,0,len_Marg-1,K-2)-\
        trans(pc.slice(t).submat(0,0,K-2,K-1)*            \
        trans(dpm_dtheta.slice(t-1)));
      ddelta_dtheta.slice(t) = trans(solve(H,trans(s1),solve_opts::fast));
    }
    
    //+++++++++++++++++++++++++++++++++++++++++++//
    //               Gradient                    //
    //+++++++++++++++++++++++++++++++++++++++++++//
    
    //------------- theta -------------//
    // dlogL1_dtheta: len_Marg * 1
    arma::vec tmpsumdL1(len_Marg, fill::zeros);
    for(k=1;k<K-1;k++){
      tmpsumdL1 += y_log(0,k)/(gamm(0,k)-gamm(0,k-1))*\
        (dgamm_dtheta.slice(0).col(k)-dgamm_dtheta.slice(0).col(k-1));
    }
    dlogL_dtheta = y_log(0,0)/gamm(0,0)*dgamm_dtheta.slice(0).col(0) +  \
      tmpsumdL1 + y_log(0,K-1)/(1-gamm(0,K-2))*                         \
      (0-dgamm_dtheta.slice(0).col(K-2));
    
    // dlogL2_dtheta: len_Marg * 1
    for(t=1;t<Ti;t++){
      dlogL_dtheta += ddelta_dtheta.slice(t)*\
        (trans(y_log.row(t).subvec(0,K-2))-  \
        pc.slice(t).col(m(t-1)).subvec(0,K-2));
    }
    
    //+++++++++++++++++++++++++++++++++++++++++++//
    //                Hessian                    //
    //+++++++++++++++++++++++++++++++++++++++++++// 
    
    // ddgamm_dthetatheta: (len_Marg) * (len_Marg) * (K-1)
    // index of , gamm_tk(k), theta_star(p), theta (q)
    arma::cube ddgamm_dthetatheta(K-1,len_Marg,len_Marg,fill::zeros);
    for(p=0;p<len_Marg;p++){
      for(q=0;q<len_Marg;q++){
        ddgamm_dthetatheta.slice(q).col(p) = trans((1-2*gamm.row(0))%  \
          gamm.row(0)%(1-gamm.row(0))%deta_dtheta.slice(0).row(p)%     \
          deta_dtheta.slice(0).row(q));
      }
    }
    
    //-------------  theta by theta  -------------//
    // ddlogL1_ddtheta: len_Marg * len_Marg
    for(p=0;p<len_Marg;p++){
      for(q=p;q<len_Marg;q++){
        tmpsum1=0.0;
        for(k=1;k<K-1;k++){
          tmpsum1 += y_log(0,k)*(1/pow((gamm(0,k)-gamm(0,k-1)),2.0)*\
            (dgamm_dtheta(p,k,0)-dgamm_dtheta(p,k-1,0))*            \
            (dgamm_dtheta(q,k,0)-dgamm_dtheta(q,k-1,0)) -           \
            1/(gamm(0,k)-gamm(0,k-1))*                              \
            (ddgamm_dthetatheta(k,p,q)-ddgamm_dthetatheta(k-1,p,q)));
        }
        ddlogL_ddtheta(p,q) = y_log(0,0)*(1/(pow(gamm(0,0),2.0))* \
          dgamm_dtheta(p,0,0)*                                    \
          dgamm_dtheta(q,0,0)-                                    \
          1/gamm(0,0)*                                            \
          ddgamm_dthetatheta(0,p,q)) +                            \
          tmpsum1 +                                               \
          y_log(0,K-1)*(1/(pow(gamm(0,K-1-1),2.0))*               \
          dgamm_dtheta(p,0,K-1-1)*                                \
          dgamm_dtheta(q,0,K-1-1)+                                \
          1/gamm(0,K-1-1)*                                        \
          ddgamm_dthetatheta(K-1-1,p,q));
      }
    }
    
    ddlogL_ddtheta += trans(ddlogL_ddtheta);
    //Since double count the diagonal elements
    ddlogL_ddtheta.diag() = ddlogL_ddtheta.diag()/2;
    
    // ddlogL2_ddtheta: len_Marg * len_Marg
    for(t=1;t<Ti;t++){
      ddlogL_ddtheta += ddelta_dtheta.slice(t).submat(0,0,len_Marg-1,K-2)*  \
        trans(ddelta_dtheta.slice(t)*                                       \
        dpc_ddelta.slice(m(t-1)*Ti+t).submat(0,0,K-2,K-2));
    }
}

// [[Rcpp::export]]
void GradHessCppnewVec22(arma::vec& par, arma::mat& DesignX, arma::mat& delta, arma::mat& y_log, arma::vec& m, arma::vec& dlogL_dxi, arma::mat& ddlogL_ddxi, double conv, double stepNR){
  
  // extract parameters
  Ti        = y_log.n_rows;
  K         = y_log.n_cols;
  len_mean  = DesignX.n_cols; 
  len_par   = par.size();
  len_Marg  = len_mean*(K-1);
  len_Dep   = pow(K-1,2.0);
  
  // double conv = 1e-05;
  
  //-----------  Allocation Memory  -----------//
  double tmpsum1, tmpsum2, tmpsum3;
  arma::mat theta(K-1, len_mean, fill::zeros);
  arma::mat xi(K-1, K-1, fill::zeros);
  arma::mat eta(Ti, K-1, fill::zeros);
  arma::mat gamm(Ti, K-1, fill::zeros);
  arma::mat pm(Ti, K, fill::zeros);
  arma::mat H(K-1, K-1, fill::zeros);
  arma::mat s1(len_Marg, K-1, fill::zeros);
  arma::mat s2(K-1, len_Dep, fill::zeros);
  
  arma::cube pc(K, K, Ti, fill::zeros);
  
  //-----------  Allocation Memory  -----------//
  for(i=0;i<len_mean;i++){
    theta.col(i) = par.subvec(i*(K-1),(i+1)*(K-1)-1);
  }
  for(i=0;i<K-1;i++){
    xi.col(i) = par.subvec(len_Marg+i*(K-1), \
           len_Marg+i*(K-1)+(K-1)-1);
  }
  
  //----------- Marginal Parameter -----------//
  eta = DesignX*trans(theta);
  gamm = exp(eta)/(1 + exp(eta));
  pm.col(0) = gamm.col(0);
  pm.submat(0,1,Ti-1,K-2) = gamm.submat(0,1,Ti-1,K-2) - \
    gamm.submat(0,0,Ti-1,K-3);
  pm.col(K-1) = 1- gamm.col(K-2);
  
  // Check marginal probability positive
  while(bool positive_prob=any(vectorise(pm)<0)) {
    Rcout<<"Marginal Probability less than 0 !!!"<<std::endl;
    break;
  }
  
  //----- Newton-Raphson: Calculation of delta_tk values -----//
  //-----    Re-Calculation on Conditional Probability   -----//
  NRDeltaCpp(xi, pm, DesignX, delta, pc, K, conv, stepNR);
  
  //---------- Derivatives of eta_tk w.r.t theta -------------//
  // (index of theta(i), eta_tk(k), eta_tk(t))
  arma::cube deta_dtheta(len_Marg, K-1, Ti, fill::zeros);   
  for(t=0;t<Ti;t++){
    for(k=0;k<K-1;k++){
      for(i=0;i<len_mean;i++){
        deta_dtheta(i*(K-1)+k,k,t) = DesignX(t,i);
      }
    }
  }
  
  //--------- Derivatives of gamm_tk w.r.t theta -------------//
  // (index of theta(kp), gamma_tk(k), gamma_tk(tt))
  arma::mat invgamm = 1-gamm;
  arma::cube dgamm_dtheta(len_Marg, K-1, Ti, fill::zeros);  
  for(t=0;t<Ti;t++){
    for(k=0;k<K-1;k++){
      dgamm_dtheta.slice(t).col(k) = gamm(t,k)*invgamm(t,k)*\
        deta_dtheta.slice(t).col(k);
    }
  }
  
  //--------- Derivatives of pm_tk w.r.t theta --------------//
  //index of theta(kp), pm_tk(k), pm_tk(tt)
  arma::cube dpm_dtheta(len_Marg, K, Ti, fill::zeros);
  // k = 1
  dpm_dtheta.subcube(0,0,0, len_Marg-1,0,Ti-1) = \
    dgamm_dtheta.subcube(0,0,0, len_Marg-1,0,Ti-1);
  // k = 2:(K-1)
  dpm_dtheta.subcube(0,1,0, len_Marg-1,K-2,Ti-1) =     \
    dgamm_dtheta.subcube(0,1,0, len_Marg-1,K-2,Ti-1) - \
    dgamm_dtheta.subcube(0,0,0, len_Marg-1,K-3,Ti-1);
  // k = K
  dpm_dtheta.subcube(0,K-1,0, len_Marg-1,K-1,Ti-1) = \
    -dgamm_dtheta.subcube(0,K-2,0, len_Marg-1,K-2,Ti-1);
    
    //--------- Derivatives of pc_tjk w.r.t delta -------------//
    // delta_tk(kp), pc_tkj(k), pc_tkj(j)*pc_tkj(tt)
    arma::cube dpc_ddelta(K-1, K, K*Ti, fill::zeros);
    for(t=0;t<Ti;t++){
      for(j=0;j<K;j++){
        for(k=0;k<K;k++){
          for(kp=0;kp<K-1;kp++){
            if(k==kp){
              dpc_ddelta(kp,k,j*Ti+t) = pc(k,j,t)*(1-pc(k,j,t));
            }else{
              dpc_ddelta(kp,k,j*Ti+t) = -pc(k,j,t)*pc(kp,j,t);
            }
          }
        }
      }
    }
    
    //--------- Derivatives of pc_tjk w.r.t xi   -------------//
    // xi(a,b), pc_tkj(k), pc_tkj(j)*pc_tkj(tt)
    arma::cube dpc_dxi(len_Dep, K, K*Ti, fill::zeros);
    for(t=0;t<Ti;t++){
      for(a=0;a<K-1;a++){
        for(b=0;b<K-1;b++){
          for(j=0;j<K;j++){
            if(b==j){
              for(k=0;k<K-1;k++){
                if(k==a){
                  dpc_dxi(b*(K-1)+a,k,j*Ti+t) = pc(k,j,t)*(1-pc(k,j,t));
                }else{
                  dpc_dxi(b*(K-1)+a,k,j*Ti+t) = -pc(k,j,t)*pc(a,j,t);
                }
              }
            }
            // k = K
            dpc_dxi(b*(K-1)+a,K-1,j*Ti+t) = -pc(K-1,j,t)*pc(a,j,t);
          }
        }
      }
    }
    
    //----------     System Equation Solve     ---------------//
    
    //-------  Derivatives of Delta_tk w.r.t theta  ----------//
    //theta(kp), delta_tk(k), delta_tk(tt)
    arma::cube ddelta_dtheta(len_Marg, K-1, Ti, fill::zeros);
    //--------  Derivatives of Delta_tk w.r.t xi   -----------//
    //xi(a,b), delta_tk(k), delta_tk(tt)
    arma::cube ddelta_dxi(len_Dep, K-1, Ti, fill::zeros);
    
    for(t=1;t<Ti;t++){
      // Left Hand Side of eq1 & eq2
      H.zeros();
      for(j=0;j<K;j++){
        H += trans(dpc_ddelta.slice(j*Ti+t).submat(0,0,K-2,K-2)*pm(t-1,j));
      }
      // Right Hand Side of eq1 (s1: theta(kp), pm_tk(k))
      s1 = dpm_dtheta.slice(t).submat(0,0,len_Marg-1,K-2)-\
        trans(pc.slice(t).submat(0,0,K-2,K-1)*            \
        trans(dpm_dtheta.slice(t-1)));
      ddelta_dtheta.slice(t) = trans(solve(H,trans(s1),solve_opts::fast));
      
      // Left Hand Side of eq2 (s2: pc_tjk(k) & xi_kj(a,b))
      s2.zeros();
      for(j=0;j<K;j++){
        s2 -= trans(dpc_dxi.slice(j*Ti+t).submat(0,0,len_Dep-1,K-2)*pm(t-1,j));
      }
      ddelta_dxi.slice(t) = trans(solve(H,s2,solve_opts::fast));
    }
    
    //+++++++++++++++++++++++++++++++++++++++++++//
    //               Gradient                    //
    //+++++++++++++++++++++++++++++++++++++++++++//
    
    //-------------  xi  -------------//
    // dlogL1_dxi: (K-1)*(K-1) * 1 (have to be zero)
    // dlogL2_dxi: (K-1)*(K-1) * 1
    for(a=0;a<K-1;a++){
      for(b=0;b<K-1;b++){
        for(t=1;t<Ti;t++){
          dlogL_dxi(b*(K-1)+a) +=                         \
            sum((y_log.row(t).subvec(0,K-2)-              \
            trans(pc.slice(t).col(m(t-1)).subvec(0,K-2)))%\
            ddelta_dxi.slice(t).row(b*(K-1)+a))+          \
            (y_log(t,a)-pc(a,m(t-1),t))*y_log(t-1,b);
        }
      }
    }
    
    //+++++++++++++++++++++++++++++++++++++++++++//
    //                Hessian                    //
    //+++++++++++++++++++++++++++++++++++++++++++// 
    //-------------      xi by xi    -------------// 
    // ddlogL2_ddxi: (K-1)*(K-1) * (K-1)*(K-1)
    int tmpindex;
    for(t=1;t<Ti;t++){
      tmpindex = m(t-1)*Ti+t;
      for(p=0;p<len_Dep;p++){
        for(q=p;q<len_Dep;q++){
          a = q%(K-1); // get the remainder
          b = ceil(q/(K-1));
          
          ddlogL_ddxi(q,p) += 
            sum(ddelta_dxi.slice(t).row(p)*
            dpc_ddelta.slice(tmpindex).submat(0,0,K-2,K-2)%  \
            ddelta_dxi.slice(t).row(q))+                     \
            sum(dpc_dxi.slice(tmpindex).row(p).subvec(0,K-2)%\
            ddelta_dxi.slice(t).row(q))+                     \
            (sum(ddelta_dxi.slice(t).row(p)*                 \
            dpc_ddelta.slice(tmpindex).col(a))+              \
            dpc_dxi(p,a,tmpindex))*                          \
            y_log(t-1,b);
        }
      }
    }
    ddlogL_ddxi += trans(ddlogL_ddxi);
    //Since double count the diagonal elements
    ddlogL_ddxi.diag() = ddlogL_ddxi.diag()/2;
}

// [[Rcpp::export]]
List SeqNewtonCpp(arma::vec par, arma::mat& y_log, arma::vec& mm, arma::mat& DesignX, arma::vec stepsize, int maxit, double mytol){
  
  //--------------------------- extract parameters ----------------------//
  Ti        = y_log.n_rows;
  K         = y_log.n_cols;
  len_mean  = DesignX.n_cols;
  len_par   = par.size();
  len_Marg  = (K-1)*len_mean;
  len_Dep   = (K-1)*(K-1);
  
  //--------------------------- Declar identifiers ----------------------//
  arma::vec phi0 = par;
  arma::vec phi1 = par;
  arma::vec phi0_out = par;
  arma::vec phi1_out = par;
  
  int iconv_NW  = 0;
  int iconv_NW1 = 0;
  int iconv_NW2 = 0;
  int convergence;
  
  double dif_NW  = 0.01;
  double dif_NW1 = 0.01;
  double dif_NW2 = 0.01;
  
  double stepNW = stepsize(0);
  double stepNR = stepsize(1);
  
  arma::uvec upId1 = arma::linspace<arma::uvec>(0, len_Marg-1, len_Marg);
  arma::uvec upId2 = arma::linspace<arma::uvec>(len_Marg, len_par-1, len_Dep);
  
  arma::mat delta(Ti,K-1,fill::zeros);
  // Gradient
  arma::vec dlogL_dtheta(len_Marg, fill::zeros);
  arma::mat ddlogL_ddtheta(len_Marg, len_Marg, fill::zeros);
  // Hessian
  arma::vec dlogL_dxi(len_Dep, fill::zeros);
  arma::mat ddlogL_ddxi(len_Dep, len_Dep, fill::zeros);
  
  //-------------------- Sequential Newton's Methods -------------------//
  while((iconv_NW<maxit) & (dif_NW>mytol)){
    
    // Updating Mean Parameter
    iconv_NW1 = 0; 
    dif_NW1   = 0.01;
    while((iconv_NW1<maxit) & (dif_NW1>mytol)){
      dlogL_dtheta.zeros();
      ddlogL_ddtheta.zeros();
      GradHessCppnewVec21(phi0, DesignX, delta, y_log, mm, dlogL_dtheta, ddlogL_ddtheta, mytol, stepNR);
      phi1(upId1) = phi0(upId1) + inv(ddlogL_ddtheta)*dlogL_dtheta*stepNW;
      dif_NW1 = pow(sum(square(phi1-phi0)), 0.5);
      iconv_NW1 +=1;
      phi0 = phi1;
    }
    // Check Convergence
    phi1_out = phi0;
    iconv_NW += iconv_NW1;
    dif_NW = pow(sum(square(phi0_out - phi1_out)), 0.5);
    phi0_out = phi1_out;
    if((iconv_NW>maxit) | (dif_NW<mytol)){
      break;
    }
    
    // Updating Dependence Parameter
    iconv_NW2 = 0; 
    dif_NW2   = 0.01;
    while((iconv_NW2<maxit) & (dif_NW2>mytol)){
      dlogL_dxi.zeros();
      ddlogL_ddxi.zeros();
      GradHessCppnewVec22(phi0, DesignX, delta, y_log, mm, dlogL_dxi, ddlogL_ddxi, mytol, stepNR);
      phi1(upId2) = phi0(upId2) + inv(ddlogL_ddxi)*dlogL_dxi*stepNW;
      dif_NW2 = pow(sum(square(phi1-phi0)), 0.5);
      iconv_NW2 +=1; 
      phi0 = phi1;
    }
    // Check Convergence
    phi1_out = phi0;
    iconv_NW += iconv_NW2;
    dif_NW = pow(sum(square(phi0_out - phi1_out)), 0.5);
    phi0_out = phi1_out;
  }
  
  //--------------------------- Return Result -------------------------//
  List res;
  res["parEst"] = phi0_out;
  res["counts"] = iconv_NW;
  if(iconv_NW < maxit){
    res["convergence"] = 0;
  }else{
    res["convergence"] = 1;
  }
  
  return(res);
}

// [[Rcpp::export]]
List MTMResidCpp(arma::vec parEst, arma::mat& y_log, arma::vec& mm, arma::mat& DesignX, int K, int Ti, double mytol, double stepNR){
  
  // extract parameters
  len_mean  = DesignX.n_cols; 
  len_par   = parEst.size();
  len_Marg  = len_mean*(K-1);
  len_Dep   = pow(K-1,2.0);
  
  //-----------  Allocation Memory  -----------//
  arma::mat theta(K-1, len_mean, fill::zeros);
  arma::mat xi(K-1, K-1, fill::zeros);
  arma::mat eta(Ti, K-1, fill::zeros);
  arma::mat gamm(Ti, K-1, fill::zeros);
  arma::mat pm(Ti, K, fill::zeros);
  arma::mat delta(Ti, K-1, fill::zeros);
  
  arma::cube pc(K, K, Ti, fill::zeros);
  
  for(i=0;i<len_mean;i++){
    theta.col(i) = parEst.subvec(i*(K-1),(i+1)*(K-1)-1);
  }
  for(i=0;i<K-1;i++){
    xi.col(i) = parEst.subvec(len_mean*(K-1)+i*(K-1), \
           len_mean*(K-1)+i*(K-1)+(K-1)-1);
  }
  
  //----------- Marginal Parameter -----------//
  eta = DesignX*trans(theta);
  gamm = exp(eta)/(1 + exp(eta));
  
  // pm.zeros();
  for(i=0;i<K-1;i++){
    if(i==0){
      pm.col(i) = gamm.col(i);
    }else{
      pm.col(i) = gamm.col(i) - gamm.col(i-1);
    }
  }
  pm.col(K-1) = 1- gamm.col(K-2);
  
  // Check marginal probability positive
  while(bool positive_prob=any(vectorise(pm)<0)) {
    Rcout<<"Marginal Probability less than 0 !!!"<<std::endl;
    break;
  }
  
  //----- Newton-Raphson: Calculation of delta_tk values -----//
  //-----    Re-Calculation on Conditional Probability   -----//
  NRDeltaCpp(xi, pm, DesignX, delta, pc, K, mytol, stepNR);
  
  //--------  Multivariate One-step-ahead Prediction  -------//
  
  // One-step-ahead Prediction: conditional expectation on passing values
  arma::mat EXW(Ti, K, fill::zeros);
  EXW.row(0) = pm.row(0);
  for(t=1;t<Ti;t++){
    for(k=0;k<K;k++){
      EXW(t,k) = pc(k,mm(t-1),t);
    }
  }
  
  // normal residuals
  arma::mat Err = y_log - pm;
  
  // One-step-ahead prediction residuals
  arma::mat ErrW = y_log - EXW;
  
  // Weighted One-step-ahead prediction residuals
  arma::mat tmpmycov(K, K, fill::zeros);
  arma::mat ResidW(Ti, K-1, fill::zeros);
  for(t=1;t<Ti;t++){
    tmpmycov.zeros();
    // Covariance matrix
    for(k=0;k<K;k++){
      for(kp=0;kp<K;kp++){
        if(k==kp){
          tmpmycov(k,kp) = pm(t,k) - sum(pm.row(t-1)%pc.slice(t).row(k)%pc.slice(t).row(k));
        }else{
          tmpmycov(k,kp) = -pm(t,k)*pm(t,kp)                                                            \
          - (pm(t,kp)*sum(pc.slice(t).row(k)%pm.row(t-1))-pm(t,kp)*sum(pc.slice(t).row(k)%pm.row(t-1))) \
          - (pm(t,k)*sum(pc.slice(t).row(kp)%pm.row(t-1))-pm(t,k)*sum(pc.slice(t).row(kp)%pm.row(t-1))) \
          + (sum(pc.slice(t).row(k)%pc.slice(t).row(kp)%pm.row(t-1))                                    \
          - sum(pm.row(t-1)%pc.slice(t).row(k))*sum(pm.row(t-1)%pc.slice(t).row(kp)));
        }
      }
    }
    // eigenvalue decomposition to get square root of variance matrix
    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, tmpmycov.submat(0,0,K-2,K-2));
    arma::mat tmpsqrt(K-1, K-1, fill::zeros);
    tmpsqrt.diag() = sqrt(eigval);
    arma::mat rooti = arma::inv(eigvec*tmpsqrt*inv(eigvec));
    ResidW.row(t) = ErrW.row(t).subvec(0,K-2)*rooti;
  }
  
  
  //----------           Return Result         ------------//
  List res;
  res["Pred"]         = pm;
  res["Resid"]        = Err;
  res["OneStepPred"]  = EXW;
  res["OneStepResid"] = ErrW;
  res["OneStepResidW"]= ResidW;
  
  return(res);
  
}

// [[Rcpp::export]]
List MTMResidCpp1(arma::vec parEst, arma::mat& y_log, arma::vec& mm, arma::mat& DesignX, int K, int Ti, double mytol, double stepNR){
  
  // extract parameters
  len_mean  = DesignX.n_cols; 
  len_par   = parEst.size();
  len_Marg  = len_mean*(K-1);
  len_Dep   = pow(K-1,2.0);
  
  //-----------  Allocation Memory  -----------//
  arma::mat theta(K-1, len_mean, fill::zeros);
  arma::mat xi(K-1, K-1, fill::zeros);
  arma::mat eta(Ti, K-1, fill::zeros);
  arma::mat gamm(Ti, K-1, fill::zeros);
  arma::mat pm(Ti, K, fill::zeros);
  arma::mat delta(Ti, K-1, fill::zeros);
  
  arma::cube pc(K, K, Ti, fill::zeros);
  
  for(i=0;i<len_mean;i++){
    theta.col(i) = parEst.subvec(i*(K-1),(i+1)*(K-1)-1);
  }
  for(i=0;i<K-1;i++){
    xi.col(i) = parEst.subvec(len_mean*(K-1)+i*(K-1), \
           len_mean*(K-1)+i*(K-1)+(K-1)-1);
  }
  
  //----------- Marginal Parameter -----------//
  eta = DesignX*trans(theta);
  gamm = exp(eta)/(1 + exp(eta));
  
  // pm.zeros();
  for(i=0;i<K-1;i++){
    if(i==0){
      pm.col(i) = gamm.col(i);
    }else{
      pm.col(i) = gamm.col(i) - gamm.col(i-1);
    }
  }
  pm.col(K-1) = 1- gamm.col(K-2);
  
  // Check marginal probability positive
  while(bool positive_prob=any(vectorise(pm)<0)) {
    Rcout<<"Marginal Probability less than 0 !!!"<<std::endl;
    break;
  }
  
  //----- Newton-Raphson: Calculation of delta_tk values -----//
  //-----    Re-Calculation on Conditional Probability   -----//
  NRDeltaCpp(xi, pm, DesignX, delta, pc, K, mytol, stepNR);
  
  //--------  Multivariate One-step-ahead Prediction  -------//
  
  // One-step-ahead Prediction: conditional expectation on passing values
  arma::mat EXW(Ti, K, fill::zeros);
  EXW.row(0) = pm.row(0);
  for(t=1;t<Ti;t++){
    for(k=0;k<K;k++){
      EXW(t,k) = pc(k,mm(t-1),t);
    }
  }
  
  // normal residuals
  arma::mat Err = y_log - pm;
  
  // One-step-ahead prediction residuals
  arma::mat ErrW = y_log - EXW;
  
  // Weighted One-step-ahead prediction residuals
  arma::mat tmpmycov(K, K, fill::zeros);
  arma::mat ResidW(Ti, K-1, fill::zeros);
  for(t=1;t<Ti;t++){
    tmpmycov.zeros();
    // Covariance matrix
    for(k=0;k<K;k++){
      for(kp=0;kp<K;kp++){
        if(k==kp){
          tmpmycov(k,kp) = pm(t,k) - 2*sum(pm.row(t)%pc.slice(t).row(k)%pm.row(t-1))\
          + sum(pm.row(t-1)%pc.slice(t).row(k)%pc.slice(t).row(k));
        }else{
          tmpmycov(k,kp) = -pm(t,k)*pm(t,kp)                                                            \
          - (pm(t,kp)*sum(pc.slice(t).row(k)%pm.row(t-1))-pm(t,kp)*sum(pc.slice(t).row(k)%pm.row(t-1))) \
          - (pm(t,k)*sum(pc.slice(t).row(kp)%pm.row(t-1))-pm(t,k)*sum(pc.slice(t).row(kp)%pm.row(t-1))) \
          + (sum(pc.slice(t).row(k)%pc.slice(t).row(kp)%pm.row(t-1))                                    \
          - sum(pm.row(t-1)%pc.slice(t).row(k))*sum(pm.row(t-1)%pc.slice(t).row(kp)));
        }
      }
    }
    // eigenvalue decomposition to get square root of variance matrix
    vec eigval;
    mat eigvec;
    eig_sym(eigval, eigvec, tmpmycov.submat(0,0,K-2,K-2));
    arma::mat tmpsqrt(K-1, K-1, fill::zeros);
    tmpsqrt.diag() = sqrt(eigval);
    arma::mat rooti = arma::inv(eigvec*tmpsqrt*inv(eigvec));
    ResidW.row(t) = ErrW.row(t).subvec(0,K-2)*rooti;
  }
  
  
  //----------           Return Result         ------------//
  List res;
  res["Pred"]         = pm;
  res["Resid"]        = Err;
  res["OneStepPred"]  = EXW;
  res["OneStepResid"] = ErrW;
  res["OneStepResidW"]= ResidW;
  
  return(res);
  
}
