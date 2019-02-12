PredictY_Pr=function(xdata, theta,Sigma)
{

  theta[,1:3]=exp(theta[,1:3])
  
  mu1= theta[,1]+ (theta[,2]*xdata[[1]]/(theta[,3]+xdata[[1]]))
  
  mu2= 1/(1 + exp(-theta[,4]-theta[,5]*xdata[[1]]))
 
  Y1=rnorm(length(mu1),mu1,Sigma)
  Y2=rbinom(length(mu2),1,mu2)

  y=data.frame(Y1,Y2)
  return(y)
  
}