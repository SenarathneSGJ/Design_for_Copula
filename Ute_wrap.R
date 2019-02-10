Utility_wrap<-function(d,B)
{
  theta=B[[1]]
  W=B[[2]]
  post_model_probs=B[[3]]
  LogZs=B[[4]]
  
  theta1=theta[,1:9]
  theta2=theta[,10:17]
  theta_i=list(theta1,theta2)
  Likehood_function=c(likelihood_m1,likelihood_m2)
  ute_est=c()
  
  for(m in 1:2)
  {
    ute_est[m]=ute_kld(data=d,theta=theta_i[[m]],W=W[,m],likelihood_Func=Likehood_function[[m]])
  }
  
  utility_Est=sum(post_model_probs*ute_est)  
  
  utility_Disc=ute_mid(designT=d,theta=theta,W=W,post_model_probs=post_model_probs,LogZs=LogZs) 
  utility= utility_Disc+utility_Est
  #print(utility)
  
  return(utility)
  
}  