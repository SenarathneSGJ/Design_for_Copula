
likelihood_m2_Gum=function(data,para,Sigma) 
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
  z1= (Y1-mu1)/Sigma
  pi1=pnorm(z1)
  
  datax=data.matrix(cbind(x= rep(1,nrow(data)),data[,1]))
  paranew=data.matrix(para[,4:5])
  Dy2=paranew%*%t(datax)
  
  pi2 = 1/(1 + exp(-Dy2)) 
  
  alpha_Gu=1/(1-para[,6])
  
  G_Copula=exp(-1*(((-log(pi1))^alpha_Gu)+((-log(pi2))^alpha_Gu))^(1/alpha_Gu))
  G_Copula[G_Copula<0]=0
  G_Copula[G_Copula>1]=1
  
  G_Copula=G_Copula*(1/pi1)*((-log(pi1))^(alpha_Gu-1))*((((-log(pi1))^alpha_Gu)+((-log(pi2))^alpha_Gu))^((1/alpha_Gu)-1))
  G_Copula[is.na(G_Copula) | G_Copula<0]=0
  G_Copula[G_Copula>1]=1
  
  log_likelihood = dnorm(Y1,mu1,Sigma,log = TRUE)
  
  log_likelihood[Y2 == 1] = log_likelihood[Y2 == 1] + log(G_Copula[Y2 == 1])
  log_likelihood[Y2 == 0] = log_likelihood[Y2 == 0] + log(1-G_Copula[Y2 == 0])
  
  log_likelihood[log_likelihood<(-5e+2)]<-(-5e+2)

  return(log_likelihood)
}


