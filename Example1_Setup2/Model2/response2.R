
response2=function(xdata)
{
  mean_theta=c(1,4,1,0,1,-0.5,1,0,14.14) #actual parameters
  ydata=c()
  
    y1= mean_theta[1]+ (mean_theta[2]*xdata[1])+ (mean_theta[3]*xdata[2])+ (mean_theta[4]*xdata[3])
    y2= mean_theta[5]+ (mean_theta[6]*xdata[1])+ (mean_theta[7]*xdata[2])+ (mean_theta[8]*xdata[3])
    
    pi1 = 1/(1 + exp(-y1))
    pi2 = 1/(1 + exp(-y2))
    
    frankC=-1*(mean_theta[9]^(-1))*log(1+((exp(-mean_theta[9]*pi1)-1)*(exp(-mean_theta[9]*pi2)-1)/(exp(-mean_theta[9])-1)))
    
    p11= frankC
    p10=pi1-p11
    p01=pi2-p11
    p00=1-pi1-pi2+p11
    
  
    u1=runif(1)
    
    if(u1<p00)
    {
      y=c(0,0)
    }else if(u1<p00+p01)
    {
      y=c(0,1)
      
    }else if(u1<p00+p01+p10)
    {
      y=c(1,0)
      
    }else
    {
      y=c(1,1)
      
    }
    
    ydata=c(ydata,y)
    return(y)
  
}