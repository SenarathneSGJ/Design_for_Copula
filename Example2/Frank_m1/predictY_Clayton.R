PredictY_Cl=function(xdata, theta,Sigma)
{

  theta[,1:3]=exp(theta[,1:3])
  
  mu1= theta[,1]+ (theta[,2]*xdata[[1]]/(theta[,3]+xdata[[1]]))
  
  mu2= 1/(1 + exp(-theta[,4]-theta[,5]*xdata[[1]]))
  pi2=mu2
  
  alpha_Cl= 2*theta[,6]/(1-theta[,6])
  
  Y1=rnorm(length(mu1),mu1,Sigma)
  pi1=pnorm(Y1,mu1,Sigma)
    
  Cl_Copula= ((pi1^(-alpha_Cl))+(pi2^(-alpha_Cl))-1)^(-1/alpha_Cl)
  Cl_Copula[Cl_Copula<0]=0
	Cl_Copula[Cl_Copula>1]=1
    
  Cl_Copula=Cl_Copula*(pi1^(-alpha_Cl-1))/((pi1^(-alpha_Cl))+(pi2^(-alpha_Cl))-1) 
  Cl_Copula[is.na(Cl_Copula) | Cl_Copula<0]=0
  Cl_Copula[Cl_Copula>1]=1
    
  Y2=rbinom(length(Cl_Copula),1,Cl_Copula)
  y=data.frame(Y1,Y2)
  
  return(y)
  
}