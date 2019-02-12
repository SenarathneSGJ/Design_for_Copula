PredictY_Gum=function(xdata, theta,Sigma)
{

  theta[,1:3]=exp(theta[,1:3])
  
  mu1= theta[,1]+ (theta[,2]*xdata[[1]]/(theta[,3]+xdata[[1]]))
  
  mu2= 1/(1 + exp(-theta[,4]-theta[,5]*xdata[[1]]))
  pi2=mu2
 
  alpha_Gu=1/(1-theta[,6])
  
  Y1=rnorm(length(mu1),mu1,Sigma)
  pi1=pnorm(Y1,mu1,Sigma)
    
  G_Copula=exp(-1*(((-log(pi1))^alpha_Gu)+((-log(pi2))^alpha_Gu))^(1/alpha_Gu))
	G_Copula[G_Copula<0]=0
  G_Copula[G_Copula>1]=1
	
  G_Copula=G_Copula*(1/pi1)*((-log(pi1))^(alpha_Gu-1))*((((-log(pi1))^alpha_Gu)+((-log(pi2))^alpha_Gu))^((1/alpha_Gu)-1))
  G_Copula[G_Copula<0]=0
  G_Copula[G_Copula>1]=1
   
  Y2=rbinom(length(G_Copula),1,G_Copula)
  
  y=data.frame(Y1,Y2)
  return(y)
  
}