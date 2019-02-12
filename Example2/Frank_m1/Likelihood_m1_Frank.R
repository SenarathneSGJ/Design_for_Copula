
likelihood_m1_Fr=function(data,para,Sigma) 
{
  alpha_Fr=c()
  
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
  
  for(i in 1:nrow(para))
  {
    #alpha_Fr[i]=uniroot(fn2, lower=-1e+10, upper=1e+10,k=para[i,6])$root
    alpha_Fr[i]=multiroot(fn2, start=.01,k=para[i,6],positive=TRUE,useFortran = FALSE)$root
  }
  
  F_Copula= (exp(-alpha_Fr*pi1))*(exp(-alpha_Fr*pi2)-1)/((exp(-alpha_Fr)-1)+(exp(-alpha_Fr*pi1)-1)*(exp(-alpha_Fr*pi2)-1))
  F_Copula[F_Copula>1]=1
  F_Copula[F_Copula<0]=0
  
  log_likelihood = dnorm(Y1,mu1,Sigma,log = TRUE)
  
  log_likelihood[Y2 == 1] = log_likelihood[Y2 == 1] + log(F_Copula[Y2 == 1])
  log_likelihood[Y2 == 0] = log_likelihood[Y2 == 0] + log(1-F_Copula[Y2 == 0])
 
  log_likelihood[log_likelihood<(-5e+2)]<-(-5e+2)
  
  
  return(log_likelihood)
}


