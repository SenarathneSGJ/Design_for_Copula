mh4=function(theta,N,W,datat,R)
{
  theta=theta[sample(nrow(theta),N,prob=W,replace=T), ]   # resampling
  
  
  mu_eps=rep(0,ncol(theta))
  
  for (m in 1:R)
  {
    
    Sigma_eps=cov(theta)
    eps= rmvnorm(n = nrow(theta), mu_eps, Sigma_eps,method="svd")
    theta_prop=theta+eps 
    
    #Loglkhd_j= rowSums(likelihood_m4_Pr(data=datat,para=theta,Sigma=Sigma))
    #Loglkhd_prop=rowSums(likelihood_m4_Pr(data=datat,para=theta_prop,Sigma=Sigma))
    
    Loglkhd_j= rowSums(LogLike_m4_Pr(d=datat[,1],Y=as.matrix(datat[,2:3]),para=as.matrix(theta),Sigma=Sigma))
    Loglkhd_prop=rowSums(LogLike_m4_Pr(d=datat[,1],Y=as.matrix(datat[,2:3]),para=as.matrix(theta_prop),Sigma=Sigma))
    
    # Calculate the likelihood of theta and theta* under the prior
    log_prior_theta_j <- dmvnorm((theta[,1:5]),mean=Exp_theta,sigma=diag(c(rep(SD^2,3),9,9)),log=TRUE) 
    log_prior_theta_prop<- dmvnorm((theta_prop[,1:5]),mean=Exp_theta,sigma=diag(c(rep(SD^2,3),9,9)),log=TRUE) 
    
    # Compute likelihood on the log-scale
    
    r <- exp(Loglkhd_prop - Loglkhd_j + log_prior_theta_prop - log_prior_theta_j)
       
    rxy=ifelse(r>1,1,r)
    r1=runif(N)
    replace=which(rxy>r1)
    theta[replace,]=theta_prop[replace,]  
    
  }
  
  return(theta)  
}