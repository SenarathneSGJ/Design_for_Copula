mh1=function(theta,N,W,datat,R)
{
  theta=theta[sample(nrow(theta),N,prob=W,replace=T), ]   # resampling
  
  
  mu_eps=rep(0,ncol(theta))
  
  for (m in 1:R)
  {
    
    Sigma_eps=cov(theta)
    eps= rmvnorm(n = nrow(theta), mu_eps, Sigma_eps,method="svd")
    theta_prop=theta+eps 
    
    ind.alpha=c()
    ##### Find indicator for alpha
    ind.alpha=ifelse(((theta_prop[,9] > -0.99 & theta_prop[,9] < -0.05) | (theta_prop[,9] > 0.05 & theta_prop[,9] < 0.99)),1,0)
    vec=which(ind.alpha==0)
    theta_prop[vec,9]=theta[vec,9]
    
    Loglkhd_j= rowSums(likelihood_m1(data=datat,para=theta))
    Loglkhd_prop=rowSums(likelihood_m1(data=datat,para=theta_prop))
    
    # Calculate the likelihood of theta and theta* under the prior
    
    
    log_prior_theta_j <- dmvnorm((theta[,1:8]),mean=rep(0,8),sigma=diag(rep(SD^2,8)),log=TRUE) 
    log_prior_theta_prop<- dmvnorm((theta_prop[,1:8]),mean=rep(0,8),sigma=diag(rep(SD^2,8)),log=TRUE)  
    
    # Compute likelihood on the log-scale
    
    r <- exp(Loglkhd_prop - Loglkhd_j + log_prior_theta_prop - log_prior_theta_j)
    
    rxy=ifelse(r>1,1,r)
    r1=runif(N)
    replace=which(rxy>r1 & ind.alpha==1)
    theta[replace,]=theta_prop[replace,]    
    
  }
   
return(theta)  
}
