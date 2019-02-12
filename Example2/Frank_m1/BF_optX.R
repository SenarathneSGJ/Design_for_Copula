BF_optX=function(theta,W,post_model_probs,LogZs) # Next design point via ACE
{

  tempX <- data.frame(X=seq(0.05,2,by=0.05)) 
  print(tempX)

  utes=foreach(i = 1:nrow(tempX),.packages= c("mvtnorm","rootSolve","stats","Rcpp"),.noexport = c("LogLike_m2_Gum","LogLike_m3_Cl","LogLike_m4_Pr"),.export=c("Sigma","u_Entropy","PredictY_Fr","PredictY_Gum","fn2","integrand","PredictY_Cl","PredictY_Pr","likelihood_m1_Fr"),.verbose=TRUE,.combine = c) %dopar% 
  {
    u_Entropy(design_all=tempX[i,],theta=theta,W=W,post_model_probs=post_model_probs,LogZs=LogZs,B=5000)
  }
  print(utes)
  item=which.max(utes)
  OptX=tempX[item,]   
  
  out=list(X=OptX,utility=utes[item],utes_t=utes)

  return(out)
  
}