
SMC_Func=function(Starting_designs,theta,W,R,Sigma,T1)
{  
  N<-nrow(theta)
  datat=data.frame()
  ute_t=data.frame()
  prop_data=data.frame()
  num_models=4
  utility=data.frame()
  LogZs = rep(0,num_models)  
  post_model_probs=t(rep((1/num_models),num_models))
  post_model_probs_all=data.frame()
   Log_det=c()
  
  for(t1 in 0:(T1-1))
  {
    #ptm=Sys.time()   
    opt_data=BF_optX(theta=theta,W=W,post_model_probs=post_model_probs,LogZs=LogZs)
    xdata=opt_data[[1]]
    
    ute_t=rbind(ute_t,opt_data[[3]])
    ################################################################################
    # Generate data from true model (xdata, true parameter values), stochastic
    ydata=response_Fr(xdata=xdata,Sigma=Sigma)
    ################################################################################

    prop_data= merge(xdata,ydata)
    datat=rbind(datat,prop_data)
    
    w1=exp(likelihood_m1_Fr(data=datat[t1+1,],para=theta[,1:6],Sigma=Sigma))       
    #w2=exp(likelihood_m2_Gum(data=datat[t1+1,],para=theta[,7:12],Sigma=Sigma)) 
    w2=exp(LogLike_m2_Gum(d=datat[t1+1,1],Y=as.matrix(datat[t1+1,2:3]),para=as.matrix(theta[,7:12]),Sigma=Sigma))
    #w3=exp(likelihood_m3_Cl(data=datat[t1+1,],para=theta[,13:18],Sigma=Sigma))
    w3=exp(LogLike_m3_Cl(d=datat[t1+1,1],Y=as.matrix(datat[t1+1,2:3]),para=as.matrix(theta[,13:18]),Sigma=Sigma))
    #w4=exp(likelihood_m4_Pr(data=datat[t1+1,],para=theta[,19:23],Sigma=Sigma)) 
    w4=exp(LogLike_m4_Pr(d=datat[t1+1,1],Y=as.matrix(datat[t1+1,2:3]),para=as.matrix(theta[,19:23]),Sigma=Sigma)) 
     
    w=data.frame(w1,w2,w3,w4)
    W<- w*W    # Unnormalised importance weights
    
    #update the marginal likelihoods for each model
    W_temp=c(log(colSums(W)))
    LogZs = LogZs + (W_temp)
    
    #update the posterior model probabilities
    post_model_probs = exp(LogZs - max(LogZs))
    post_model_probs = post_model_probs/sum(post_model_probs)
    post_model_probs=data.matrix(t(post_model_probs))
    post_model_probs_all = rbind(post_model_probs_all,post_model_probs)
    
   
     # Normalised importance weights
    W[,1]=W[,1]/colSums(W)[1]
    W[,2]=W[,2]/colSums(W)[2]
    W[,3]=W[,3]/colSums(W)[3]
    W[,4]=W[,4]/colSums(W)[4]
    
    ESS= 1/(colSums(W^2))
    
    if(ESS[1]<.75*N | t1==T1-1)
    {
      theta1_new=mh1(theta=theta[,1:6],N=N,W=W[,1],datat=datat,R=R)
      theta[,1:6] <- theta1_new   
      W[,1] <- rep(1,N)/N     
    }  
     
    if(ESS[2]<.75*N | t1==T1-1)
    {
      theta2_new=mh2(theta=theta[,7:12],N=N,W=W[,2],datat=datat,R=R)
      theta[,7:12] <- theta2_new 
      W[,2] <- rep(1,N)/N     
    } 
    
    if(ESS[3]<.75*N | t1==T1-1)
    {
      theta3_new=mh3(theta=theta[,13:18],N=N,W=W[,3],datat=datat,R=R)
      theta[,13:18] <- theta3_new 
      W[,3] <- rep(1,N)/N     
    }  
    
    if(ESS[4]<.75*N | t1==T1-1)
    {
      theta4_new=mh4(theta=theta[,19:23],N=N,W=W[,4],datat=datat,R=R)
      theta[,19:23] <- theta4_new 
      W[,4] <- rep(1,N)/N     
    } 
    
    print(t1)
	print(post_model_probs)
    utility[t1+1,1]=u_Entropy(design_all=xdata, theta=theta, W=W,post_model_probs=post_model_probs,LogZs=LogZs,B=5000) 
    #print(utility[t1+1,1])
	  Log_det[t1+1]=crit_cov(theta_vals=theta[,1:6],theta_w=W[,1])
    #print(Log_det[t1+1])
  
  }

 datanew=cbind(theta,W) 
 names(ESS)=c("Model1_ESS","Model2_ESS","Model3_ESS","Model4_ESS")
 names(post_model_probs_all)=c("post_model_probs_m1","post_model_probs_m2","post_model_probs_m3","post_model_probs_m4")
 out=list(ESS,datanew,datat,post_model_probs_all,utility,ute_t,Log_det)
 return(out)
 
}

