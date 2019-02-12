
library("stats")
library("mvtnorm")
library("rootSolve")
library("Rcpp")
library("doParallel")

source('utility_Entropy.R')
source('BF_optX.R')
source('Likelihood_m1_Frank.R')
sourceCpp("LikeCpp_m2_Gum.cpp")
sourceCpp("LikeCpp_m3_Cl.cpp")
sourceCpp("LikeCpp_m4_Pr.cpp")

source('response_Frank.R')
source('SMC_Func.R')
source('mh1.R')
source('mh2.R')
source('mh3.R')
source('mh4.R')
source('predictY_Gumbel.R')
source('predictY_Frank.R')
source('predictY_Clayton.R')
source('predictY_Product.R')
source('crit_cov.R')
N=5000
SD=0.5
Sigma=0.4
Exp_theta=c(2,13,0.25,0,0)
Exp_theta[1:3]=log(Exp_theta[1:3])

integrand <- function(t) {t/(exp(t)-1)}
fn2=function(alpha,k){1-(4/alpha)+(4/alpha^2)*(integrate(integrand,lower=0,upper=alpha))$value-k}

Log_Eo=rnorm(N,0,SD)+Exp_theta[1]
Log_Emax=rnorm(N,0,SD)+Exp_theta[2]
Log_ED50=rnorm(N,0,SD)+Exp_theta[3]
b0=rnorm(N,0,3)
b1=rnorm(N,0,3)
Tau=c(runif(N,0.01,0.99))

theta1=data.frame(Log_Eo,Log_Emax,Log_ED50,b0,b1,Tau,row.names=NULL)  
theta2=data.frame(Log_Eo,Log_Emax,Log_ED50,b0,b1,Tau,row.names=NULL) 
theta3=data.frame(Log_Eo,Log_Emax,Log_ED50,b0,b1,Tau,row.names=NULL) 
theta4=data.frame(Log_Eo,Log_Emax,Log_ED50,b0,b1,row.names=NULL) 

# Particle Set
theta=cbind(theta1,theta2,theta3,theta4)  

W1=c(rep((1/N),N))
W2=c(rep((1/N),N))
W3=c(rep((1/N),N))
W4=c(rep((1/N),N))

W=data.frame(W1,W2,W3,W4)

vec=SMC_Func(theta=theta, W=W,R=30,Sigma=Sigma,T1=600)

v_name=paste("vec","RData",sep=".")
save(vec, file = v_name)
