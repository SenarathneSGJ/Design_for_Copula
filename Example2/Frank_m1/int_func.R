out=function(alpha){
  fn=list()
  for(i in 1:length(k))
  {
    alpha1 <- alpha[i]
    k1=k[i]
    fn[i]=function(alpha1){1-(4/alpha1)+(4/alpha1^2)*(integrate(integrand,lower=0,upper=alpha1))$value-k1}
  }
}


makefunc <- function(k) { k; function(alpha) 1-(4/alpha)+(4/alpha^2)*(integrate(integrand,lower=0,upper=alpha))$value-k}

funclist <- list()
for (i in 1:length(kk)) funclist[[i]] <- makefunc(kk[i])
makefunc <- function(k,i) { k;i; function(alpha) 1-(4/alpha[i])+(4/alpha[i]^2)*(integrate(integrand,lower=0,upper=alpha[i]))$value-k}




funclist <- list()
for (i in 1:10) funclist[[i]] <- makefunc(para[1:10,6],i)

testFun <- function(alpha){
  
  sapply(1:length(alpha),function(i){funclist[[i]](alpha)})
  
}

multiroot(testFun, start=rep(.01,nrow(para)),useFortran = FALSE)$root