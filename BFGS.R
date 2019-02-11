BFGS=function(theta,W,post_model_probs,LogZs) # Next design point via BFGS
{
  
  optX=data.frame()
  utes=c()
  vec.choice <- seq(-1,1,by=1)
  num.v <- length(vec.choice)
  
  possible_designs=data.frame(x1=rep(vec.choice,each=(num.v^2)),X2=rep(vec.choice,each=num.v,times=num.v),x3=rep(vec.choice,times=(num.v^2)))
  
  for(k in 1: nrow(possible_designs))
    {
      utes[k]=Utility_wrap(X=possible_designs[k,],theta=theta,W=W,post_model_probs=post_model_probs,LogZs=LogZs)[[1]]
  }  
  
  Ute_data=cbind(possible_designs,utes)
  
  ind_max <- which.max(Ute_data$utes)
  opt_X= Ute_data[ind_max,]
    
  return(data.frame(opt_X))
  
}
