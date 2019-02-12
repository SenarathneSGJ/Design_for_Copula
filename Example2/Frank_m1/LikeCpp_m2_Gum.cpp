#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]

NumericMatrix LogLike_m2_Gum(NumericVector d,NumericMatrix Y, NumericMatrix para, double Sigma) {
  int n = d.size(), B = para.nrow();
  NumericVector fy1(n),fy2(n),u1(n),u2(n),GuCopula(n),Log_GuCopula(n),G_cop(n),log_dy1(n);
  NumericMatrix log_likelihood(B,n);
  double alpha;
 
for(int j = 0; j < B; j++){
 NumericVector theta = para(j,_);  
  for(int i = 0; i < n; i++) {
    fy1[i] = exp(theta[0])+(exp(theta[1])*d[i])/(exp(theta[2])+d[i]);  
    u1[i] = R::pnorm(Y(i,0),fy1[i],Sigma,1,0);
    
	fy2[i] = theta[3]+(theta[4]*d[i]);
    u2[i] = 1/(1+exp(-fy2[i]));
	
    alpha=(1/(1-theta[5]));
    G_cop[i]=exp(-1*pow(((pow(-log(u1[i]),alpha))+(pow(-log(u2[i]),alpha))),(1/alpha)));
    if(G_cop[i] < 0)
      G_cop[i]= 0;
	if(G_cop[i] > 1)
      G_cop[i]= 1;
    GuCopula[i]=G_cop[i]*(1/u1[i])*pow(-log(u1[i]),(alpha-1))* pow(((pow(-log(u1[i]),alpha))+(pow(-log(u2[i]),alpha))),((1/alpha)-1)); 
	
	if(ISNAN(GuCopula[i]))
	  GuCopula[i]= 0;
	if(GuCopula[i] < 0)
      GuCopula[i]= 0;
	if(GuCopula[i] > 1)
      GuCopula[i]= 1;
	
    log_dy1[i]=R::dnorm(Y(i,0),fy1[i],Sigma,1);
    
    if(Y(i,1)==1)
      log_likelihood(j,i)= log_dy1[i]+ log(GuCopula[i]);
    if(Y(i,1)==0)
      log_likelihood(j,i)= log_dy1[i]+ log(1-GuCopula[i]);

    if(log_likelihood(j,i) < (-5e+2))
      log_likelihood(j,i) = (-5e+2);
	} 
  }
  return log_likelihood;   
}