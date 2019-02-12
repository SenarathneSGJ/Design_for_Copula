
likelihood_m2=function(data,para) 
{
  alpha_Fr=c()
  log_likelihood=matrix(ncol=nrow(data),nrow=nrow(para))
  para1=data.matrix(para[,1:3])
  datax=data.matrix(cbind(x= rep(1,nrow(data)),data[1:2]))
  Dy1=para1%*%t(datax)
  
  para2=data.matrix(para[,4:6])
  Dy2=para2%*%t(datax)
  
  pi1 = 1/(1 + exp(-Dy1)) 
  pi2 = 1/(1 + exp(-Dy2))
  
  for(i in 1:nrow(para))
  {
    #alpha_Fr[i]=uniroot(fn2, lower=-1e+10, upper=1e+10,k=para[i,6])$root
    alpha_Fr[i]=multiroot(fn2, start=.1,k=para[i,7],maxiter=250)$root
  }
  
  frankC=-1*(alpha_Fr^(-1))*log(1+((exp(-alpha_Fr*pi1)-1)*(exp(-alpha_Fr*pi2)-1)/(exp(-alpha_Fr)-1)))
  frankC[frankC==Inf | frankC>pi1]<-pi1[frankC==Inf | frankC>pi1]
  frankC[frankC>pi2]<-pi2[frankC>pi2]
  
  p11= frankC
  p10=pi1-p11
  p01=pi2-p11
  p00=1-pi1-pi2+p11
  p00[p00<0]=0
  
  for(j in 1:nrow(data))
  {
    
    if(data[j,4]==1 & data[j,5]==1)
    {
      log_likelihood[,j]=log(p11[,j])
      
    }else if(data[j,4]==1 & data[j,5]==0)
    {
      log_likelihood[,j]=log(p10[,j])
      
    }else if(data[j,4]==0 & data[j,5]==1)
    {
      log_likelihood[,j]=log(p01[,j])
      
    }else
    {
      log_likelihood[,j]=log(p00[,j])
    }  
  }
  log_likelihood[log_likelihood<(-1e+100)]=(-1e+100)
  
  return(log_likelihood)
}


