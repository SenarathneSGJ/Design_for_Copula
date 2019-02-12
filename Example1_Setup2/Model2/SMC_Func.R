
SMC_Func=function(T1,theta,W,R)
{  
  N<-nrow(theta)
  datat=data.frame()
  prop_data=data.frame()
  ute_t=data.frame()
  log_det=c()
  post_model_probs_all=data.frame()
  utility=matrix(ncol=3,nrow=T1)
  num_models=2
  LogZs = rep(0,num_models)  
  post_model_probs=data.frame(0.5,0.5)
  N1=500
  W3=data.frame(w31=c(rep((1/N1),N1)),w32=c(rep((1/N1),N1)))
    
  for(t1 in 0:(T1-1))
  {
    xdata=c()
    ydata=c()
    theta1s=(theta[,1:9])
    theta2s=(theta[,10:16])
    theta1s=theta1s[sample(nrow(theta1s),N1,prob=W[,1],replace=T), ]
    theta2s=theta2s[sample(nrow(theta2s),N1,prob=W[,2],replace=T), ]
    thetas=data.frame(theta1s,theta2s)
    
    B.list <- list(thetas,W3,post_model_probs,LogZs)
    
    #opt_data1=BFGS(theta=thetas,W=W3,post_model_probs=post_model_probs,LogZs=LogZs)
    opt_data=ace(Utility_wrap, start.d=matrix(c(0,0,0),1), B=list(B.list,B.list), Q = 12, N1 = 5, N2 = 0, lower=-1, upper = 1, deterministic = TRUE)
    xdata=data.frame(opt_data$phase1.d)
    ute_t=rbind(ute_t,t(opt_data$phase1.trace))
   
    # Generate data from true model
    ydata=data.frame(response2(xdata=xdata),row.names=c("y1","y2"))
    
    prop_data= merge(xdata,t(ydata))
    datat=rbind(datat,prop_data)
    
    w1=exp(likelihood_m1(data=datat[t1+1,],para=theta[,1:9]))       
    w2=exp(likelihood_m2(data=datat[t1+1,],para=theta[,10:16]))     
    w=data.frame(w1,w2)
    # Unnormalised importance weights:
	W<- w*W             
    
    #update the marginal likelihoods for each model
    LogZs = LogZs + log(colSums(W))
    
    #update the posterior model probabilities
    post_model_probs = exp(LogZs - max(LogZs))
    post_model_probs = post_model_probs/sum(post_model_probs)
    post_model_probs_all = rbind(post_model_probs_all,post_model_probs)
    # Normalised importance weights:
    W <- W/colSums(W) 
    
    ESS= 1/(colSums(W^2))
 
    if(ESS[1]<.75*N | t1==T1-1)
    {
      theta1=mh1(theta=theta[,1:9],N=N,W=W[,1],datat=datat,R=R)
      theta[,1:9] <- theta1   
      W[,1] <- rep(1,N)/N     
    }  
     
    if(ESS[2]<.75*N | t1==T1-1)
    {
      theta2=mh2(theta=theta[,10:16],N=N,W=W[,2],datat=datat,R=R)
      theta[,10:16] <- theta2 
      W[,2] <- rep(1,N)/N     
    }  
    
    print(t1)
    utility[t1+1,1:3]=Utility_wrap(d=xdata,B=list(theta,W,post_model_probs,LogZs))
    log_det[t1+1]=crit_cov(theta_vals=theta[,10:16],theta_w=W[,2])
  }

 datanew=cbind(theta,W) 
 names(ESS)=c("Model1_ESS","Model2_ESS")
 names(post_model_probs_all)=c("post_model_probs_m1","post_model_probs_m2")
 out=list(ESS,datanew,datat,post_model_probs_all,utility,log_det,ute_t)
 return(out)
 
}

