
likelihood_m3_Cl=function(data,para,Sigma) 
{
  log_likelihood=matrix(ncol=nrow(data),nrow=nrow(para))
  para[,1:3]=exp(para[,1:3])
  
  X1=matrix(rep(data[,1],nrow(para)),ncol=nrow(data),byrow=T)
  Para1= matrix(rep(para[,1],nrow(data)),ncol=nrow(data))
  Para2=data.matrix(para[,2])
  Para3= matrix(rep(para[,3],nrow(data)),ncol=nrow(data))
  
  mu1=Para1+((Para2%*%(t(data[,1])))/(Para3+X1))    
  mu1=matrix(mu1,ncol=nrow(data))
  
  Y1= matrix(rep(data[,2],nrow(para)),ncol=nrow(data),byrow=T)
  Y2= matrix(rep(data[,3],nrow(para)),ncol=nrow(data),byrow=T)
  pi1=pnorm(Y1,mu1,Sigma)
  
  datax=data.matrix(cbind(x= rep(1,nrow(data)),data[,1]))
  paranew=data.matrix(para[,4:5])
  Dy2=paranew%*%t(datax)
  
  pi2 = 1/(1 + exp(-Dy2)) 
  alpha_Cl= 2*para[,6]/(1-para[,6])
  
  Cl_Copula= ((pi1^(-alpha_Cl))+(pi2^(-alpha_Cl))-1)^(-1/alpha_Cl)
  Cl_Copula[Cl_Copula<0]=0
  Cl_Copula[Cl_Copula>1]=1
  
  Cl_Copula=Cl_Copula*(pi1^(-alpha_Cl-1))/((pi1^(-alpha_Cl))+(pi2^(-alpha_Cl))-1) 
  Cl_Copula[is.na(Cl_Copula) | Cl_Copula<0]=0
  Cl_Copula[Cl_Copula>1]=1
  
  log_likelihood = dnorm(Y1,mu1,Sigma,log = TRUE)
  
  log_likelihood[Y2 == 1] = log_likelihood[Y2 == 1] + log(Cl_Copula[Y2 == 1])
  log_likelihood[Y2 == 0] = log_likelihood[Y2 == 0] + log(1-Cl_Copula[Y2 ==0])
  
  log_likelihood[log_likelihood<(-5e+2)]<-(-5e+2)

  return(log_likelihood)
}


