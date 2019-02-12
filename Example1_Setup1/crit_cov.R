crit_cov=function(theta_vals,theta_w)
{
  
  var_theta <- cov.wt(theta_vals,c(theta_w))
  
  var_theta=var_theta$cov
    
  crit = log(det(var_theta))
  return(crit)
  
}
