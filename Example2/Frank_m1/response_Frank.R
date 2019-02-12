
response_Fr=function(xdata,Sigma)
{
  mean_theta=c(2.5, 14.5, 0.2, -2.123,3.728,14.14) #actual parameters
  #para=(Eo,Emax,ED50,Bo,B1,alpha)
  
  mu1= mean_theta[1]+ (mean_theta[2]*xdata[[1]]/(mean_theta[3]+xdata[[1]]))

  mu2= 1/(1 + exp(-mean_theta[4]-mean_theta[5]*xdata[[1]]))
  
  Y1=rnorm(1,mu1,Sigma)
  pi1=pnorm(Y1,mu1,Sigma)
  pi2=mu2
  
  F_Copula=(exp(-mean_theta[6]*pi1))*(exp(-mean_theta[6]*pi2)-1)/((exp(-mean_theta[6])-1)+(exp(-mean_theta[6]*pi1)-1)*(exp(-mean_theta[6]*pi2)-1))
  F_Copula[F_Copula>1]=1
  F_Copula[F_Copula<0]=0
  
  Y2=rbinom(1,1,F_Copula)
  
  y=data.frame(Y1,Y2)
  return(y)
  
}