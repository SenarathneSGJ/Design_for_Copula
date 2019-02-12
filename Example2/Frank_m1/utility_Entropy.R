u_Entropy=function(design_all, theta, W,post_model_probs,LogZs,B) 
{  
  crit=c()
  crit_Est=c()
  Ud_m=c()
  w=data.frame()
  np_ute=300
  
  theta1=theta[,1:6]
  theta2=theta[,7:12]
  theta3=theta[,13:18]
  theta4=theta[,19:23]
  
  theta11=theta1[sample(nrow(theta1),np_ute,prob=W[,1],replace=T), ]
  theta21=theta2[sample(nrow(theta2),np_ute,prob=W[,2],replace=T), ]
  theta31=theta3[sample(nrow(theta3),np_ute,prob=W[,3],replace=T), ]
  theta41=theta4[sample(nrow(theta4),np_ute,prob=W[,4],replace=T), ]
  
  Y1=PredictY_Fr(xdata=design_all,theta=theta11,Sigma=Sigma)
  Y2=PredictY_Gum(xdata=design_all,theta=theta21,Sigma=Sigma)
  Y3=PredictY_Cl(xdata=design_all,theta=theta31,Sigma=Sigma)
  Y4=PredictY_Pr(xdata=design_all,theta=theta41,Sigma=Sigma)
  
  theta1=theta1[sample(nrow(theta1),B,prob=W[,1],replace=T), ]
  theta2=theta2[sample(nrow(theta2),B,prob=W[,2],replace=T), ]
  theta3=theta3[sample(nrow(theta3),B,prob=W[,3],replace=T), ]
  theta4=theta4[sample(nrow(theta4),B,prob=W[,4],replace=T), ]
  
  Y= list(Y1,Y2,Y3,Y4)
  dataX=rep(design_all,np_ute)
  
  for(m in 1:4)
  {
	data_select=data.frame(X=dataX,Y[[m]])
    # dataframe with ncol=B(no of particles)  and nrow=np_ute(no of design points)
    log_w1=likelihood_m1_Fr(data=data_select,para=theta1,Sigma=Sigma)  
	#log_w2=likelihood_m2_Gum(data=data_select,para=theta2,Sigma=Sigma)
    log_w2=LogLike_m2_Gum(d=dataX,Y=as.matrix(Y[[m]]),para=as.matrix(theta2),Sigma=Sigma)
	#log_w3=likelihood_m3_Cl(data=data_select,para=theta3,Sigma=Sigma)  
    log_w3=LogLike_m3_Cl(d=dataX,Y=as.matrix(Y[[m]]),para=as.matrix(theta3),Sigma=Sigma)
	#log_w4=likelihood_m4_Pr(data=data_select,para=theta4,Sigma=Sigma)
    log_w4=LogLike_m4_Pr(d=dataX,Y=as.matrix(Y[[m]]),para=as.matrix(theta4),Sigma=Sigma)
	
    log_W_temp<- cbind(log_w1,log_w2,log_w3,log_w4)      #store all the log(weights) in a dataframe
    W_temp<-exp(log_W_temp)
    
    # Unnormalised importance weights
    W_temp<- W_temp*(1/B)
    probs_all= c(colSums(W_temp))
    
    W_temp1=matrix(rep(t(probs_all),B),nrow=B,byrow=T)
    W_temp2 <- W_temp/W_temp1     # Normalised importance weights
    
    crit1 = colSums(W_temp2*log_W_temp)-log(colSums(W_temp)) 
    crit1[probs_all==0]=0
    crit1= crit1[(np_ute*(m-1)+1):(np_ute*m)]
    crit_Est[m]=mean(crit1)
    
    for(j in 1:np_ute)
    {
      probs=c(probs_all[j],probs_all[j+np_ute],probs_all[j+2*np_ute],probs_all[j+3*np_ute])
      log_Z = LogZs + log(probs)
      p = log_Z - max(log_Z)
      p = exp(p)          
      
      p = p/sum(p)
      crit[j] = log(p[m])
    }
    
    Ud_m[m] = mean(crit)
  }
  
    utility_Disc=sum(post_model_probs*Ud_m)        # discrimination utility
    utility_Est=sum(post_model_probs*crit_Est)     # Estimation utility
    utility= utility_Disc+utility_Est
	return(utility)
  
}