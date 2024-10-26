gdata<-function(alpha1,beta1,theta,n,m,c1,c2){
  #this program aims to generate length biased sample {X,A,delta} 
  #and right censored sample {Y,d} 
  #input:
  #  c1,c2--- are censoring rate for X0 and tY separately
  
  X0 = rgamma(n,shape = alpha1+1, rate = beta1)
  A = NULL
  for(i in 1:n){
    A[i] = runif(n=1,max=X0[i])
  }
  C = rexp(n, rate = c1)
  X = pmin(X0,A+C)
  delta = (A+C>=X0)*1          # 1 dead;0 censor
  tY = rexp(m,rate = theta)
  V = rexp(m,rate = c2)
  Y = pmin(tY,V)
  d = (V>=tY)*1
  list(X=X,A=A,delta=delta,Y=Y,d=d)
}
