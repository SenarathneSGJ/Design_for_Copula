source('Ute_wrap.R')
source('Ute_kld.R')
source('Ute_mid.R')
source('crit_cov.R')
source('Likelihood_m1.R')
source('Likelihood_m2.R')
source('response1.R')
source('SMC_Func.R')
source('mh1.R')
source('mh2.R')

library("stats")
library("mvtnorm")
library("acebayes")
library("rootSolve")

N=5000
SD=4

AI <- Sys.getenv("casenumber")
NAI <- as.numeric(AI)
NAI=NAI
set.seed(NAI*1000)

# Prior for model 1
b01=rnorm(N,0,SD)
b11=rnorm(N,0,SD)
b21=rnorm(N,0,SD)
b31=rnorm(N,0,SD)
b02=rnorm(N,0,SD)
b12=rnorm(N,0,SD)
b22=rnorm(N,0,SD)
b32=rnorm(N,0,SD)
Tau=c(runif(N/2,-0.99,-0.05),runif(N/2,0.05,0.99))

theta1=data.frame(b01,b11,b21,b31,b02,b12,b22,b32,Tau,row.names=NULL)  #paricle set

#Prior for model 2
b01=rnorm(N,0,SD)
b11=rnorm(N,0,SD)
b21=rnorm(N,0,SD)
b31=rnorm(N,0,SD)
b02=rnorm(N,0,SD)
b12=rnorm(N,0,SD)
b22=rnorm(N,0,SD)
b32=rnorm(N,0,SD)

theta2=data.frame(b01,b11,b21,b31,b02,b12,b22,b32,row.names=NULL) 

theta=data.frame(theta1,theta2)

Starting_designs=data.frame(X1=runif(250,-1,1),X2=runif(250,-1,1),X3=runif(250,-1,1)) #Starting design points in each iteration

seeds=5
ptm <- proc.time()
W1=c(rep((1/N),N))
W2=c(rep((1/N),N))
W=data.frame(W1,W2)

integrand <- function(t) {t/(exp(t)-1)}
fn2=function(alpha,k){1-(4/alpha)+(4/alpha^2)*(integrate(integrand,lower=0,upper=alpha))$value-k}

vec=SMC_Func(Starting_designs=Starting_designs,theta=theta,W=W,R=30) 

ptm2=proc.time() 
ptm2-ptm

NAII=NAI
v_name=paste("vec",NAII,sep="")
v_name2=paste(v_name,"RData",sep=".")

save(vec, file = v_name2)


