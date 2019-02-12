ute_kld=function(data, theta, W,likelihood_Func) #data= current design point , theta= particle set, W=weights
{  
  
  crit=c()
  probs=c()
  
  Y=data.frame(y1=c(0,0,1,1),y2=c(0,1,0,1)) # select four outcomes (Y1,Y2)=(00,01,10,11)
  X=data[rep(seq_len(nrow(data)), each=4),]
  
  data_select=cbind(X,Y)
  
  
  for(j in 1:4) # run for four outcomes (00,01,10,11)
  {
    
    log_w <- (likelihood_Func(data=data_select[j,],para=theta))  # likelihood_Func returns log_likelihood values
    w=exp(log_w)
    
    W_temp1= w*W                     # Unnormalised importance weights
    probs[j]=sum(W_temp1)
    W_temp2 = W_temp1/sum(W_temp1)   # Normalised importance weights
    crit[j] = sum(W_temp2*log_w)-log(sum(W_temp1))       #U(d,y|.)

  }
  
  
  ud=sum(crit*probs)    #U(d|.)
  return(ud)
  
}

