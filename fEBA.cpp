#include <RcppArmadillo.h>

#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <sstream>
#include <cmath>

//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
arma::cx_cube ghat(arma::mat X, int N, int K, bool stdz){
  
  //initialization
  int Ts=X.n_rows;
  int R=X.n_cols;
  int B=floor(Ts/N);
  int Fs=floor(N/2)+1;
  int drp=floor(Ts%N/2+0.5);
  
  arma::mat xmat=join_horiz(ones(N,1),linspace(1,N,N));
  arma::mat vec(N,1);
  arma::mat linfit(2,1);
  arma::mat tmp(N,K);
  arma::mat tapers=zeros(N,K);
  
  arma::cx_mat tmp2(N,R);
  arma::cx_mat tmp3(N,1);
  arma::cx_mat tmp4(N,pow(R,2));
  arma::cx_mat tmpfft(N,K);

  arma::cx_cube lmtse(B,R,R);
  arma::cx_cube fftcb(N,K,R);
  arma::cx_cube tmp5(N,pow(R,2),K);
  arma::cx_cube mtspec(N,pow(R,2),B);
  
  //sine tapers
  for(int i=0;i<N;i++){
    for(int j=0;j<K;j++){
      tapers(i,j)=sqrt(2/double(N+1))*sin(M_PI*(i+1)*(j+1)/double(N+1));
    }
  }
  
  for(int i=0;i<B;i++){
    for (int j=0;j<R;j++){
      //linearly detrend
      vec=X(span(i*N+drp,(i+1)*N+drp-1),j);
      linfit=inv(xmat.t()*xmat)*xmat.t()*vec;
      vec=vec-xmat*linfit;
      //standard to unit variance
      if(stdz){vec=vec/repelem(stddev(vec),N,1);}
      //fft
      tmp=tapers;
      tmp.each_col() %= vec;
      fftcb.slice(j)=fft(tmp);
    }
    
    //multitaper spectral estimator
    for(int k=0;k<K;k++){
      tmp2=fftcb.subcube(0,k,0,N-1,k,R-1);
      for(int j=0;j<R;j++){
        tmp3=conj(fftcb.subcube(0,k,j,N-1,k,j));
        tmp4.submat(0,j*R,N-1,(j+1)*R-1)=sqrt(repmat(tmp3,1,R) % tmp2);
      }
    tmp5.slice(k)=tmp4;
    }
    mtspec.slice(i)=mean(tmp5,2);
  }
  
  //throw warning if observations discarded
  ostringstream text;
  if(Ts%N!=0){
    if(Ts%N%2!=0){
      text << "Warning: T is not a multiple of N. First " 
           << drp << " and last " << drp+1 << " observations have been discarded.";
      warning(text.str());
    }else if(Ts%N%2==0){
      text << "Warning: T is not a multiple of N. First " 
           << drp << " and last " << drp << " observations have been discarded.";
      warning(text.str());
    }
  }
  return mtspec;
}

// WORKING CODE AS OF 5/1/2020
// // [[Rcpp::export]]
// arma::cx_cube partNcpp(arma::mat X, int N, int K, bool stdz){
//   int Ts=X.n_rows;
//   int R=X.n_cols;
//   int B=floor(Ts/N);
//   int Fs=floor(N/2)+1;
//   
//   arma::cx_cube lmtse(B,R,R);
//   arma::cx_cube tmpmt(N,K,R);
//   arma::cx_cube tmpmtc(N,pow(R,2),K);
//   arma::cx_cube tmpmtc2(N,pow(R,2),B);
//   arma::cx_mat tmpmt2(N,R);
//   arma::cx_mat tmpmt4(N,pow(R,2));
//   arma::cx_mat mtspec(N,pow(R,2));
//   
//   //partition time series and detrend/standardize
//   Rcpp::List Xarr(B);
//   arma::mat xmat=join_horiz(ones(N,1),linspace(1,N,N));
//   arma::mat vec(N,1);
//   arma::mat linfit(2,1);
//   arma::mat tmp(N,K);
//   arma::cx_mat tmpfft(N,K);
//   int drp=floor(Ts%N/2+0.5);
//   arma::vec freq=regspace(0,1/N,floor(Ts/2)/Ts);
//   
//   //sine tapers
//   arma::mat tapers=zeros(N,K);
//   for(int i=0;i<N;i++){
//     for(int j=0;j<K;j++){
//       tapers(i,j)=sqrt(2/double(N+1))*sin(M_PI*(i+1)*(j+1)/double(N+1));
//     }
//   }
//   
//   for(int i=0;i<B;i++){
//     for (int j=0;j<R;j++){
//       //linearly detrend
//       vec=X(span(i*N+drp,(i+1)*N+drp-1),j);
//       linfit=inv(xmat.t()*xmat)*xmat.t()*vec;
//       vec=vec-xmat*linfit;
//       //standard to unit variance
//       if(stdz){vec=vec/repelem(stddev(vec),N,1);}
//       //fft
//       tmp=tapers;
//       tmp.each_col() %= vec;
//       tmpmt.slice(j)=fft(tmp);
//     }
//     //multitaper estimate
//     for(int i=0;i<K;i++){
//       tmpmt2=tmpmt.subcube(0,i,0,N-1,i,R-1);
//       for(int j=0;j<R;j++){
//         arma::cx_mat tmpmt3=conj(tmpmt.subcube(0,i,j,N-1,i,j));
//         tmpmt4.submat(0,j*R,N-1,(j+1)*R-1)=sqrt(repmat(tmpmt3,1,R) % tmpmt2);
//       }
//       tmpmtc.slice(i)=tmpmt4;
//     }
//     mtspec=mean(tmpmtc,2);
//     tmpmtc2.slice(i)=mtspec;
//     
//     // tmpmt
//     Xarr[i]=tmp;
//   }
//   
//   //Throw warning if observations discarded
//   ostringstream text;
//   if(Ts%N!=0){
//     if(Ts%N%2!=0){
//       text << "Warning: T is not a multiple of N. First " 
//            << drp << " and last " << drp+1 << " observations have been discarded.";
//       warning(text.str());
//     }else if(Ts%N%2==0){
//       text << "Warning: T is not a multiple of N. First " 
//            << drp << " and last " << drp << " observations have been discarded.";
//       warning(text.str());
//     }
//   }
//   return tmpmtc2;
// }


/*** R
ghat(X,N,K,FALSE)
*/
