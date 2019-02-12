#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericMatrix LogLike_m4_Pr(NumericVector d,NumericMatrix Y, NumericMatrix para, double Sigma) {
  int n = d.size(), B = para.nrow();
  NumericVector fy1(n),fy2(n),u2(n),log_dy1(n);
  NumericMatrix log_likelihood(B,n);
 
for(int j = 0; j < B; j++){
 NumericVector theta = para(j,_);  
  for(int i = 0; i < n; i++) {
    fy1[i] = exp(theta[0])+(exp(theta[1])*d[i])/(exp(theta[2])+d[i]);  
    
	fy2[i] = theta[3]+(theta[4]*d[i]);
    u2[i] = 1/(1+exp(-fy2[i]));
	
    
    log_dy1[i]=R::dnorm(Y(i,0),fy1[i],Sigma,1);
    
    if(Y(i,1)==1)
      log_likelihood(j,i)= log_dy1[i]+ log(u2[i]);
    if(Y(i,1)==0)
      log_likelihood(j,i)= log_dy1[i]+ log(1-u2[i]);

    if(log_likelihood(j,i) < (-5e+2))
      log_likelihood(j,i) = (-5e+2);
	} 
  }
  return log_likelihood;   
}