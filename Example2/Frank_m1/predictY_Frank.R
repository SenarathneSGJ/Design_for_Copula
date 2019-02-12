PredictY_Fr=function(xdata,theta,Sigma)
{
  alpha_Fr=c()
  
  for(i in 1:nrow(theta))
  {
    #alpha_Fr[i]=uniroot(fn2, lower=-1e+10, upper=1e+10,k=theta[i,6])$root
    alpha_Fr[i]=multiroot(fn2, start=.01,k=theta[i,6],positive=TRUE,useFortran = FALSE)$root
  }
  
  theta[,1:3]=exp(theta[,1:3])
  mu1= theta[,1]+(theta[,2]*xdata[[1]]/(theta[,3]+xdata[[1]]))
  
  mu2= 1/(1 + exp(-theta[,4]-theta[,5]*xdata[[1]]))
  pi2=mu2
 
  Y1=rnorm(length(mu1),mu1,Sigma)
  pi1=pnorm(Y1,mu1,Sigma)

  F_Copula= (exp(-alpha_Fr*pi1))*(exp(-alpha_Fr*pi2)-1)/((exp(-alpha_Fr)-1)+(exp(-alpha_Fr*pi1)-1)*(exp(-alpha_Fr*pi2)-1))
  F_Copula[F_Copula<0]=0
  F_Copula[F_Copula>1]=1
    
 Y2=rbinom(length(F_Copula),1,F_Copula)
 y=data.frame(Y1,Y2)
  
  return(y)
  
}