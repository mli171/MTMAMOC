// Title: Univariate CUSUM Hyppothesis testing for Changepoint Detection
// Model: Gaussian Latent Copula Model
// Author: Mo Li
// Create Date: 05/13/2019
// Revise Date: No.8:    02/14/2020 (Speedup by Rcpp on almost everything)
//              No.9:    02/22/2020 (Correct Memory Issue related to Innovation
//                                  Algorithm when parallel)
//              No.10:   03/16/2020 (Correct Parameter Estimation under stationary Case)
//              No.11:   03/22/2020 Add Clip Prediction Function
//              No.12:   03/27/2020 Optimize Code to improve computing efficiency
//
//
// Scripts Contents:
//==============================================================================
// 1. Sequential Parameter Estimation
//    1.1 Marginal parameter estimation in vector form LSE;
//        Gradient + Hessian: ----------------------- "gr_hs"
//        Newton's Method: -------------------------- "NW_cpp"
//    1.2 Dependence parameter Estimation by MOM;
//        Sample 2nd moments: -----------------------  In R
//        Theoretical 2nd moments: ------------------ "TheorCovCal"
//        Optimize function for matching: ----------- "diffMOM"
// 2. One Step Ahead Prediction
//    2.1 Innovation Algorithm ---------------------- "UniInnovRcpp"
//    2.2 Prediction -------------------------------- "ClipPred"
//
// 3. CUSUM Test Statistic
//    3.1 Standized one-step-ahead prediction -------  "ClipPred"
//    3.2 CUSUM calculation -------------------------  "MyCusum"
//
//
//
//
// Support Function Lists:
//==============================================================================
// 1. LSE_cpp():          Sum of square calculation;
// 2. gr_hs():            Gradient and Hessian matrix calculation;
// 3. NW_cpp():           Newton's method in parameter estimation;
// 4. triangl():          Extract the lower trangular matrix element and store in
//                        an array;
// 5. pmvnorm_cpp():      Modified Rcpp version of multivaraite normal probability
//                        calculation (Require "mvtnormAPI.h");
// 6. diffMOM():          Calculate the difference between the sample 2nd moment and
//                        theoretical 2nd moment used for dependence parameter
//                        estimation;
// 7. TheorCovCal():      Theoretical 2nd moment of Clip model calculation'
// 8. KappLongRcpp():     Kappa function defined for multinomial data in univariate
//                        form;
// 9. UniInnovRcpp():     Univaraite version of Innovation algorithm, return the
//                        Innovation coefficients,
//                        Predicted mean square error,
//                        How many lags are used for prediction;
// 10. ClipPred():        Input the innovation coefficients from UniInnovRcpp(),
//                        calculate the one-step-ahead prediction on the
//                        residuals and the Observed categories;
// 11. MyCusum():         Input the innovation MSE from UniInnovRcpp() and
//                        one-step-ahead prediction residuals from ClipPred(),
//                        calculate the CUSUM test statistic and return the
//                        estimated changepoint location and CUSUM statistic
//                        value;
//
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::depends(mvtnorm)]]
#include <mvtnormAPI.h>


using namespace Rcpp;
using namespace arma;

int len_par, len_mean;
int Ts, Th, Td, K, t, i, j, k, p, n;


// [[Rcpp::export]]
double LSE_cpp(arma::vec par, arma::mat X, arma::mat DesignX){

  Th = X.n_rows;
  Td = Th/8;
  Ts = Td*24;
  K = X.n_cols;
  len_par = par.size();

  //--------------------------- Declar identifiers ----------------------//
  double h;
  arma::vec ci(K-1, fill::zeros);
  arma::vec seg(K+1, fill::zeros);
  arma::vec mst(Th, fill::zeros);
  arma::vec norm_quan;
  arma::mat EX(Th, K, fill::zeros);

  //-------------------------- Parameter Extraction ---------------------//
  ci      = par.subvec(0,K-2);
  seg(0) = R_NegInf;
  seg.subvec(1, K-1) = ci;
  seg(K) = R_PosInf;

  // mst = DesignX*par.subvec(K-1, len_par-1);
  if(len_par>K-1){
    mst = DesignX*par.subvec(K-1, len_par-1);
  }
  for(h=0;h<Th;h++){
    norm_quan = normcdf(seg - mst(h));
    EX.row(h) = trans(norm_quan.subvec(1,K) - norm_quan.subvec(0,K-1));
  }

  //-------------------------    Results Output   ----------------------//
  double LSE = accu(pow(X-EX, 2.0));
  return(LSE);
}

// [[Rcpp::export]]
double LSE_cpp1(arma::vec par, arma::mat X, arma::mat DesignX){
  
  Th = X.n_rows;
  K = X.n_cols;
  len_par = par.size();
  
  //--------------------------- Declar identifiers ----------------------//
  double h;
  double phi, beta0, beta1, beta2, mst;
  arma::vec ci(K-1, fill::zeros);
  arma::vec seg(K+1, fill::zeros);
  arma::vec norm_quan;
  arma::vec tmp(Th, fill::zeros);
  arma::mat EX(Th, K, fill::zeros);
  
  //-------------------------- Parameter Extraction ---------------------//
  ci(0)   = 0;
  ci.subvec(1,K-2) = par.subvec(0,K-3);
  seg(0)  = R_NegInf;
  seg.subvec(1, K-1) = ci;
  seg(K)  = R_PosInf;
  beta0   = par(K-2);
  beta1   = par(K-1);
  beta2   = par(K);
  phi     = par(K+1);
  
  for(t=0;t<Th;t++){
    tmp(t) = pow(phi, t);
  }
  double tmp1, tmp2, tmp3;
  for(h=0;h<Th;h++){
    tmp1 = 0;
    tmp2 = 0;
    tmp3 = 0;
    for(t=0;t<h;t++){
      tmp1 += tmp(t);
      tmp2 += tmp(t)*DesignX(h-t,1);
      tmp3 += tmp(t)*DesignX(h-t,2);
    }
    mst = beta0*tmp1 + beta1*tmp2 + beta2*tmp3;
    // norm_quan = normcdf((seg-mst)/pow(1-phi*phi, 0.5));
    norm_quan = normcdf((seg-mst)/(1/pow(1-phi*phi, 0.5)));
    EX.row(h) = trans(norm_quan.subvec(1,K) - norm_quan.subvec(0,K-1));
  }
  
  //-------------------------    Results Output   ----------------------//
  // Rcout << "par = " << par << std::endl;
  double LSE = accu(pow(X-EX, 2.0));
  return(LSE);
}

// [[Rcpp::export]]
mat gr_hs(vec par, mat X, mat DesignX){

  Th = X.n_rows;
  Td = Th/8;
  Ts = Td*24;
  K = X.n_cols;
  len_par = par.size();
  len_mean = len_par - K + 1;

  //--------------------------- Declar identifiers ----------------------//
  double h;

  vec ci, seg(K+1, fill::zeros), mst(Th, fill::zeros);
  vec norm_quan;

  mat EX(Th, K, fill::zeros);
  mat dSS_dtheta(1, len_par, fill::zeros);
  mat ddSS_ddtheta(len_par, len_par, fill::zeros);


  //--------------------------- Extract Parameters ----------------------//
  ci     = par.subvec(0,K-2);
  seg(0) = R_NegInf;
  seg.subvec(1, K-1) = ci;
  seg(K) = R_PosInf;

  //-------------------------- Expectation & Mean ----------------------//
  if(len_par>K-1){
    mst = DesignX*par.subvec(K-1, len_par-1);
  }
  // mst = DesignX*par.subvec(K-1, len_par-1);
  for(h=0;h<Th;h++){
    norm_quan = normcdf(seg - mst(h));
    EX.row(h) = trans(norm_quan.subvec(1,K) - norm_quan.subvec(0,K-1));
  }

  //------------------------- Gradient Calculation ---------------------//
  for(h=0;h<Th;h++){
    for(k=0;k<K-1;k++){
      // w.r.t ci's
      dSS_dtheta(0,k) += ((X(h,k)-EX(h,k)) - (X(h,k+1)-EX(h,k+1)))*normpdf(seg(k+1)-mst(h));
    }
    for(j=0;j<K;j++){
      // w.r.t mean parameters
      for(p=0;p<len_mean;p++){
        dSS_dtheta(0,K-1+p) += (X(h,j)-EX(h,j))*(-DesignX(h,p))*(normpdf(seg(j+1)-mst(h)) - \
          normpdf(seg(j)-mst(h)));
      }
    }
  }
  dSS_dtheta = -2*dSS_dtheta;

  //-------------------------  Hessian Calculation ---------------------//
  //------ Diag and Upper trangle ------//
  for(h=0;h<Th;h++){
    // w.r.t ci(i is 1st) by ci(j is 2nd)
    for(i=0;i<K-1;i++){
      for(j=0;j<i+1;j++){
        if(i==j){
          ddSS_ddtheta(j,i) += 2*(2*normpdf(seg(i+1)-mst(h))*normpdf(seg(i+1)-mst(h)) + \
            ((X(h,i+1)-EX(h,i+1)) - (X(h,i)-EX(h,i)))*(seg(i+1)-mst(h))*normpdf(seg(i+1)-mst(h)));
        }else if(abs(i-j)<2){
          ddSS_ddtheta(j,i) -= normpdf(seg(i+1)-mst(h))*normpdf(seg(j+1)-mst(h));
        }else{
          ddSS_ddtheta(j,i) = 0;
        }
      }
    }
    // w.r.t ci(i is 1st) by mst(j is 2nd)
    for(i=0;i<K-1;i++){
      for(j=K-1;j<len_par;j++){
        ddSS_ddtheta(i,j) += 2*(DesignX(h,j-(K-1))*normpdf(seg(i+1)-mst(h))* \
          (normpdf(seg(i+2)-mst(h))-2*normpdf(seg(i+1)-mst(h)) +             \
          normpdf(seg(i)-mst(h)) +                                           \
          ((X(h,i)-EX(h,i)) - (X(h,i+1)-EX(h,i+1)))*(seg(i+1)-mst(h))));
      }
    }
    // w.r.t mst(i is 1st) by mst(j is 2nd)
    for(i=K-1;i<len_par;i++){
      for(j=K-1;j<i+1;j++){
        for(k=0;k<K;k++){
          if(k==0){
            ddSS_ddtheta(j,i) += 2*(DesignX(h,i-(K-1))*DesignX(h,j-(K-1))*                                           \
              ((normpdf(seg(k+1)-mst(h))-0)*(normpdf(seg(k+1)-mst(h))-0) + \
              (X(h,k)-EX(h,k))*(normpdf(seg(k+1)-mst(h))*(seg(k+1)-mst(h))-0)));
          }else if(k==K-1){
            ddSS_ddtheta(j,i) += 2*(DesignX(h,i-(K-1))*DesignX(h,j-(K-1))*                                           \
              (0-normpdf(seg(k)-mst(h)))*(0-normpdf(seg(k)-mst(h))) + \
              (X(h,k)-EX(h,k))*(0-normpdf(seg(k)-mst(h))*(seg(k)-mst(h))));
          }else{
            ddSS_ddtheta(j,i) += 2*(DesignX(h,i-(K-1))*DesignX(h,j-(K-1))*                                           \
              ((normpdf(seg(k+1)-mst(h))-normpdf(seg(k)-mst(h)))*(normpdf(seg(k+1)-mst(h))-normpdf(seg(k)-mst(h))) + \
              (X(h,k)-EX(h,k))*(normpdf(seg(k+1)-mst(h))*(seg(k+1)-mst(h))-normpdf(seg(k)-mst(h))*(seg(k)-mst(h)))));
          }
        }
      }
    }
  }

  //---- lower trangular matrix element ----//
  ddSS_ddtheta += trans(ddSS_ddtheta);
  ddSS_ddtheta.diag() = ddSS_ddtheta.diag()/2;   //Since double count the diagonal elements


  //-------------------------    Results Output   ----------------------//
  mat der(1+len_par, len_par, fill::zeros);
  der.submat(0, 0, 0, len_par-1) = dSS_dtheta;
  der.submat(1, 0, len_par, len_par-1) = ddSS_ddtheta;

  return(der);
}

// [[Rcpp::export]]
List NW_cpp(vec par, mat X, mat DesignX, double stepsize, double conv){

  len_par = par.size();

  //--------------------------- Declar identifiers ----------------------//
  vec par0 = par;
  vec par1 = par0;
  int iconv_nr = 0;
  int conv_messg;

  double diff_nr = 0.01;
  // double conv = 1e-07;
  // double LSE;

  mat der(1+len_par, len_par, fill::zeros);
  mat dSS_dtheta(1, len_par, fill::zeros);
  mat ddSS_ddtheta(len_par, len_par, fill::zeros);

  while((iconv_nr<1000) & (diff_nr>conv)){
    der = gr_hs(par0, X, DesignX);
    dSS_dtheta = der.submat(0, 0, 0, len_par-1);
    ddSS_ddtheta = der.submat(1, 0, len_par, len_par-1);
    par1 = par0 - trans(dSS_dtheta*inv(ddSS_ddtheta))*stepsize;
    diff_nr = pow(sum(square(par0 - par1)), 0.5);
    // Rcout << "par0 = " << trans(par0) << std::endl;
    par0 = par1;
    // Rcout << "par1 = " << trans(par1) << std::endl;
    // LSE = LSE_cpp(par0, X, DesignX);
    // Rcout << "LSE = " << LSE << std::endl;
    iconv_nr += 1;
    //iconv_nr = 100000;
    // Rcout << "diff_nr = " << diff_nr << std::endl;
    // Rcout << "===============================" << std::endl;
  }

  //-------------------------    Results Output   ----------------------//
  List res;
  res["par"]=par0;
  // res["value"]=LSE;
  if(iconv_nr < 1000){
    conv_messg = 0;
  }else{
    conv_messg = 1;
  }
  res["convergence"]=conv_messg;
  res["gr"] = dSS_dtheta;
  res["hs"] = ddSS_ddtheta;

  return(res);
}

//[[Rcpp::export]]
vec triangl(const arma::mat& X){
  int n = X.n_cols;
  arma::vec res(n * (n-1) / 2);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < i; ++j) {
      res(j + i * (i-1) / 2) = X(i, j);
    }
  }
  return res;
}

// [[Rcpp::export]]
double pmvnorm_cpp(arma::vec& lb, arma::vec& ub, arma::vec& mu, arma::vec& lowertrivec, double abseps = 1e-3){

  int n = lb.n_elem;
  int nu = 0;
  int maxpts = 25000;     // default in mvtnorm: 25000
  double releps = 0;      // default in mvtnorm: 0
  int rnd = 1;            // Get/PutRNGstate

  double* lb_ = lb.memptr();                        // lower bound;
  double* ub_ = ub.memptr();                        // upper bound;
  double* correlationMatrix = lowertrivec.memptr(); // array of correlation coefficients;
  int* infin = new int[n];                          // Integer, array of integration limits flag;
  double* mu_ = mu.memptr();                        // array of non-centrality parameters;

  // if INFIN(I) < 0, Ith limits are (-infinity, infinity);
  // if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
  // if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
  // if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)];

  for (int i = 0; i < n; ++i) {
    if(lb(i) == R_NegInf){
      infin[i] = 0;
    }else if (ub(i) == R_PosInf){
      infin[i] = 1;
    }else{
      infin[i] = 2;
    }
  }

  // return values
  double error;
  double value;
  int inform;

  // Fortran function details check:
  // https://github.com/cran/mvtnorm/blob/master/src/mvt.f
  mvtnorm_C_mvtdst(&n, &nu, lb_, ub_,
                   infin, correlationMatrix, mu_,
                   &maxpts, &abseps, &releps,
                   &error, &value, &inform, &rnd);
                   delete[] (infin);

                   return value;
}

// [[Rcpp::export]]
double diffMOM(double par, vec Marg_est, int K, mat DesignX, double MySampleCov){

  Ts = DesignX.n_rows;

  double MyTheorCov=0.0, mydiff;
  vec ci, seg(K+1, fill::zeros), mst(Ts, fill::zeros);
  vec mymu(2), lb(2), ub(2), correle(1);
  correle(0) = par;

  //-------------------------- Parameter Extraction ---------------------//
  ci     = Marg_est.subvec(0,K-2);
  seg(0) = R_NegInf;
  seg.subvec(1, K-1) = ci;
  seg(K) = R_PosInf;
  if(len_par>K-1){
    mst = DesignX*Marg_est.subvec(K-1, len_par-1);
  }
  for(t=0;t<Ts-1;t++){
    mymu = mst.subvec(t,t+1);
    for(k=0;k<K;k++){
      lb = {seg(k), seg(k)};
      // Rcout << "Mo 1 = " << lb << std::endl;
      ub = {seg(k+1), seg(k+1)};
      // Rcout << "Mo 2 = " << ub << std::endl;
      // Rcout << "Mo 3 = " << pmvnorm_cpp(lb, ub, mymu, corrvec, 1e-3)<< std::endl;
      MyTheorCov += pmvnorm_cpp(lb, ub, mymu, correle, 1e-3);
    }
  }
  // mydiff = abs(MyTheorCov/Ts - MySampleCov);
  mydiff = pow(MyTheorCov/Ts - MySampleCov, 2);
  // Rcout << "par = " << par << std::endl;
  // Rcout << "mydiff = " << mydiff << std::endl;
  return(mydiff);
}

// [[Rcpp::export]]
double TheorCovCal(double phi, vec Marg_est, mat DesignX, int K, int Ts){

  double MyTheorCov=0.0;
  vec ci, seg(K+1, fill::zeros), mst(Ts, fill::zeros);
  NumericVector mymu(2), lb(2), ub(2);
  mat mycorr = { {1, phi}, {phi, 1} };

  SEXP JointProb;
  Environment pkg = Environment::namespace_env("mvtnorm");
  Function f = pkg["pmvnorm"];

  //-------------------------- Parameter Extraction ---------------------//
  ci     = Marg_est.subvec(0,K-2);
  seg(0) = R_NegInf;
  seg.subvec(1, K-1) = ci;
  seg(K) = R_PosInf;

  if(len_par>K-1){
    mst = DesignX*Marg_est.subvec(K-1, len_par-1);
  }
  for(t=0;t<Ts-1;t++){
    mymu = mst.subvec(t,t+1);
    for(k=0;k<K;k++){
      lb = {seg(k), seg(k)};
      ub = {seg(k+1), seg(k+1)};
      JointProb = f(lb, ub, mymu, mycorr);
      Rcpp::NumericVector vector_curr(JointProb);
      MyTheorCov += vector_curr(0);
    }
  }
  return(MyTheorCov/Ts);
}

// [[Rcpp::export]]
double KappLongRcpp(int ti, int tj, double mui, double muj, double EXi, double EXj, int K, vec seg, double phi){

  double res=0.0, prob=0.0;
  vec mymu(2), lb(2), ub(2), correle(1);
  mymu = {mui, muj};
  correle(0) = pow(phi, abs(ti-tj));

  for(i=0;i<K;i++){
    for(j=0;j<K;j++){
      lb = {seg(i), seg(j)};
      ub = {seg(i+1), seg(j+1)};
      prob = pmvnorm_cpp(lb, ub, mymu, correle, 1e-3);
      res += (i+1)*(j+1)*prob;
    }
  }
  return(res - EXi*EXj);
}

// [[Rcpp::export]]
List UniInnovRcpp(vec EXL, vec mst, vec Marg_est, double phi_est, int K, int numCoef){

  Ts = EXL.size();
  int nn;
  vec ci, seg(K+1, fill::zeros);
  ci     = Marg_est.subvec(0,K-2);
  seg(0) = R_NegInf;
  seg.subvec(1, K-1) = ci;
  seg(K) = R_PosInf;

  double tmp;
  vec V(Ts, fill::zeros);
  mat HQ(Ts, numCoef, fill::zeros);

  //---------- n = 0 ----------//
  V(0) = KappLongRcpp(0, 0, mst(0), mst(0), EXL(0), EXL(0), K, seg, phi_est);

  //---------- n = 1 ----------//
  HQ(0,0) = KappLongRcpp(1, 0, mst(1), mst(0), EXL(1), EXL(0), K, seg, phi_est)/V(0);
  tmp = V(0)*pow(HQ(0,0), 2.0);
  V(1) = KappLongRcpp(1, 1, mst(1), mst(1), EXL(1), EXL(1), K, seg, phi_est) - tmp;

  //---------- n >= 2 ----------// n=nn+1
  int mylag=0;
  double tmpcheck=100;
  for(nn=1;nn<Ts-1;nn++){
  //for(nn=1;nn<100;nn++){
    if((tmpcheck>0.00001) & (mylag==0)){
      // k=0
      HQ(nn,nn) = KappLongRcpp(nn+1, 0, mst(nn+1), mst(0), EXL(nn+1), EXL(0), K, seg, phi_est)/V(0);
      // k>0
      for(k=1;k<nn+1;k++){
        tmp = 0.0;
        for(j=0;j<k;j++){
          tmp += HQ(k-1,k-1-j)*HQ(nn,nn-j)*V(j);
        }
        HQ(nn,nn-k) = (KappLongRcpp(nn+1+1, k+1, mst(nn+1), mst(k), EXL(nn+1), EXL(k), K, seg, phi_est)-tmp)/V(k);
      }
      // Prediction MSE //
      tmp = 0.0;
      for(j=0;j<nn+1;j++){
        tmp += HQ(nn,nn-j)*HQ(nn,nn-j)*V(j);
      }
      V(nn+1) = KappLongRcpp(nn+1, nn+1, mst(nn+1), mst(nn+1), EXL(nn+1), EXL(nn+1), K, seg, phi_est) - tmp;

      tmpcheck = fabs(HQ(nn-1,nn-1));
      // set the lag=n-1 if the innovation coefficient small enough //
      if(fabs(HQ(nn,nn)) < 0.00001){
        // set lags to save computation time
        mylag = nn;
        // erase therest lags to save memory since they are all 0
        HQ.shed_cols(mylag+1, numCoef-1);
      }
    }else{
      // Only to cover the effective coefficients (lag), the rest are less than 1e-05
      for(k=nn+1-mylag;k<nn+1;k++){
        tmp = 0.0;
        for(j=std::max(k-mylag,0);j<k;j++){
            if(nn-j<mylag+1){
              tmp += HQ(k-1,k-1-j)*HQ(nn,nn-j)*V(j);
            }
        }
        HQ(nn,nn-k) = (KappLongRcpp(nn+1+1, k+1, mst(nn+1), mst(k), EXL(nn+1), EXL(k), K, seg, phi_est)-tmp)/V(k);
      }
      // Prediction MSE //
      tmp = 0.0;
      for(j=nn+1-mylag;j<nn+1;j++){
        if(nn-j<mylag+1){
          tmp += HQ(nn,nn-j)*HQ(nn,nn-j)*V(j);
        }
      }
      V(nn+1) = KappLongRcpp(nn+1, nn+1, mst(nn+1), mst(nn+1), EXL(nn+1), EXL(nn+1), K, seg, phi_est) - tmp;

      tmpcheck = fabs(HQ(nn-1,mylag));
    }
  }

  List res;
  res["HQ"]=HQ;
  res["V"]=V;
  res["lag"]=mylag;

  return(res);
}

// [[Rcpp::export]]
List ClipPred(vec EXL, vec X_hour, int mylag, mat HQ){

  Ts = EXL.size();

  int predlag;
  vec Err = X_hour - EXL;
  vec PredX(Ts, fill::zeros);
  vec PredErr(Ts, fill::zeros), ErrErr(Ts, fill::zeros);

  for(n=0;n<Ts-1;n++){
    if(mylag>0){
      predlag = mylag;
    }else{
      predlag = n;
    }
    for(j=0;j<std::min(n+1,predlag);j++){
      PredErr(n+1) += HQ(n,j)*(Err(n+1-j-1)-PredErr(n+1-j-1));
    }
  }
  ErrErr = Err - PredErr;
  PredX = PredErr + EXL;

  //---------- Return Result ----------//
  List res;
  res["Innovation"] = ErrErr;
  res["Prediction"] = PredX;
  return(res);
}

// [[Rcpp::export]]
List MyCusum(vec ErrErr, vec V){

  Ts = ErrErr.size();

  int tau;
  vec It = ErrErr/sqrt(V);
  vec Csm(Ts, fill::zeros);

  vec tmpsum1 = cumsum(It);
  double tmpsum2 = sum(It);
  double tmpsqrt = sqrt(Ts);

  for(tau=1;tau<Ts-1;tau++){
    Csm(tau) = fabs((tmpsum1(tau) - (tau+1)*tmpsum2/Ts)/tmpsqrt);
    // Rcout << "tau " << tau << std::endl;
  }
  uword CpLoc = Csm.index_max();
  double CpVal = Csm(CpLoc);

  //---------- Return Result ----------//
  List res;
  res["Csm"]      = Csm;
  res["Location"] = CpLoc + 1;
  res["Maximum"]  = CpVal;
  return(res);
}

// // [[Rcpp::export]]
// This is the original version of CUSUM calculation for backup
// List MyCusum(vec ErrErr, vec V){
//
//   Ts = ErrErr.size();
//
//   vec It = ErrErr/sqrt(V);
//   vec Csm(Ts, fill::zeros);
//
//   for(tau=1;tau<Ts-1;tau++){
//     Csm(tau) = fabs((sum(It.subvec(0,tau)) - (tau+1)*sum(It)/Ts)/sqrt(Ts));
//     // Rcout << "tau " << tau << std::endl;
//   }
//   uword CpLoc = Csm.index_max();
//   double CpVal = Csm(CpLoc);
//
//   //---------- Return Result ----------//
//   List res;
//   res["Csm"]      = Csm;
//   res["Location"] = CpLoc + 1;
//   res["Maximum"]  = CpVal;
//   return(res);
// }
