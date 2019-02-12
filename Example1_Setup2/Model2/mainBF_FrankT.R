
source('Ute_wrap.R')
source('Ute_kld.R')
source('Ute_mid.R')
source('crit_cov.R')
source('Likelihood_m1.R')
source('Likelihood_m2.R')
source('response2.R')
source('SMC_Func.R')
source('mh1.R')
source('mh2.R')

library("stats")
library("mvtnorm")
library("acebayes")
library("rootSolve")

N=5000
SD=4

# Prior 
b01=rnorm(N,0,SD)
b11=rnorm(N,0,SD)
b21=rnorm(N,0,SD)
b31=rnorm(N,0,SD)
b02=rnorm(N,0,SD)
b12=rnorm(N,0,SD)
b22=rnorm(N,0,SD)
b32=rnorm(N,0,SD)
Tau=c(runif(N,-0.99,0.99))
theta1=data.frame(b01,b11,b21,b31,b02,b12,b22,b32,Tau,row.names=NULL)  

r01=rnorm(N,0,SD)
r11=rnorm(N,0,SD)
r21=rnorm(N,0,SD)
r02=rnorm(N,0,SD)
r12=rnorm(N,0,SD)
r22=rnorm(N,0,SD)
Tau2=c(runif(N,-0.99,0.99))
theta2=data.frame(r01,r11,r21,r02,r12,r22,Tau2,row.names=NULL) 

#paricle set
theta=data.frame(theta1,theta2)

W1=c(rep((1/N),N))
W2=c(rep((1/N),N))
W=data.frame(W1,W2)

integrand <- function(t) {t/(exp(t)-1)}
fn2=function(alpha,k){1-(4/alpha)+(4/alpha^2)*(integrate(integrand,lower=0,upper=alpha))$value-k}

vec=SMC_Func(T1=250,theta=theta,W=W,R=30) 

v_name=paste("vec","RData",sep=".")
save(vec, file = v_name)


