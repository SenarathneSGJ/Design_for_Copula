#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericMatrix LogLike_m3_Cl(NumericVector d,NumericMatrix Y, NumericMatrix para, double Sigma) {
  int n = d.size(), B = para.nrow();
  NumericVector fy1(n),fy2(n),u1(n),u2(n),ClCopula(n),Log_ClCopula(n),C_cop(n),log_dy1(n);
  NumericMatrix log_likelihood(B,n);
  double alpha;
 
for(int j = 0; j < B; j++){
 NumericVector theta = para(j,_);  
  for(int i = 0; i < n; i++) {
    fy1[i] = exp(theta[0])+(exp(theta[1])*d[i])/(exp(theta[2])+d[i]);  
    u1[i] = R::pnorm(Y(i,0),fy1[i],Sigma,1,0);
    
	fy2[i] = theta[3]+(theta[4]*d[i]);
    u2[i] = 1/(1+exp(-fy2[i]));
	
    alpha=(2*theta[5]/(1-theta[5]));
	
    C_cop[i]=pow(((pow(u1[i],-alpha))+(pow(u2[i],-alpha))-1),(-1/alpha));
    if(C_cop[i] < 0)
      C_cop[i]= 0;
	if(C_cop[i] > 1)
      C_cop[i]= 1;
    ClCopula[i]=C_cop[i]*pow(u1[i],(-alpha-1))/(pow(u1[i],-alpha)+pow(u2[i],-alpha)-1);
	
	if(ISNAN(ClCopula[i]))
	  ClCopula[i]= 0;
	if(ClCopula[i] < 0)
      ClCopula[i]= 0;
	if(ClCopula[i] > 1)
      ClCopula[i]= 1;
	
    log_dy1[i]=R::dnorm(Y(i,0),fy1[i],Sigma,1);
    
    if(Y(i,1)==1)
      log_likelihood(j,i)= log_dy1[i]+ log(ClCopula[i]);
    if(Y(i,1)==0)
      log_likelihood(j,i)= log_dy1[i]+ log(1-ClCopula[i]);

    if(log_likelihood(j,i) < (-5e+2))
      log_likelihood(j,i) = (-5e+2);
	} 
  }
  return log_likelihood;   
}