
likelihood_m2=function(data,para) 
{
  log_likelihood=matrix(ncol=nrow(data),nrow=nrow(para))
  para1=data.matrix(para[,1:4])
  datax=data.matrix(cbind(x= rep(1,nrow(data)),data[1:3]))
  Dy1=para1%*%t(datax)
  
  para2=data.matrix(para[,5:8])
  Dy2=para2%*%t(datax)
  
  pi1 = 1/(1 + exp(-Dy1)) 
  pi2 = 1/(1 + exp(-Dy2))
  
  p11= pi1*pi2
  p10=pi1-p11
  p01=pi2-p11
  p00=1-pi1-pi2+p11
  
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


