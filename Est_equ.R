###GMM estimation when theta=theta_MLE
###GMM1 function is estimating function based on weight W1=XSc(X-A)
###GMM2 function is estimating function based on weight W2=\int_0^t Sc(u)du

GMM1<-function(Delta,Rd,prob,h){
  X = Rd$X
  A = Rd$A
  delta = Rd$delta
  n = length(X)
  
  SC = eKM(X-A,X-A,delta)
  W1 = X * SC
  a = delta / W1
  
  y = Rd$Y
  d = Rd$d
  #ty = Rd[,6]
  mle_theta = sum(d)/sum(y)  # exponential distribution MLE
  #ordy = t(t(sort(ty)))
  #dfy = apply(ordy,2,Dist_y,theta=mle_theta)
  #eta = fquant(dfy,prob,ordy)
  eta = qexp(prob,rate = mle_theta)
  ss = rep(eta+Delta,n)-X
  phi = pnorm(ss/h)-rep(prob,n)
  
  g1 = sum(phi * a)
  return(g1)
} 

GMM2<-function(Delta,Rd,prob,h){
  X = Rd$X
  A = Rd$A
  delta = Rd$delta
  n = length(X)
  
  W2 = paiw(X-A,X,delta)
  a2 = delta/W2
  
  y = Rd$Y
  d = Rd$d
  mle_theta = sum(d)/sum(y)
  eta = qexp(prob,rate = mle_theta)
  ss = rep(eta+Delta,n)-X
  phi = pnorm(ss/h)-rep(prob,n)
  
  g2 = sum(phi * a2)
  return(g2)
}


eKM<-function(wc,obs,delta){
#This is a program to calculate Kaplan-meier estimator
# input:
#   wc-- the points at which the KM-estimator of censored df is evalued
#   obs--- the censored data, i.e. X=min(T,C); T survival data, C censored data
#   delta--- The indicator of censored data, i.e. delta=I(T<=C);
# output:
#   censdf--denote the Kaplan-Meier of censoring survival function
  m = length(obs)
  UV= unique(obs)  # or unique(wc)?????
  dc = NULL
  S = NULL
  for(j in 1:length(UV)){
    s = sort(UV)[j]
    z = sum((obs>=rep(s,m)))
    dc[j] = 1-sum((obs==rep(s,m))*(1-delta))/(z+0.0001)
  }
  S0 = cumprod(dc)   #calculate SC for every observed data
  for(i in 1:m){
    id = sum(wc[i]>=sort(UV))
    S[i] = S0[id]
  }
  return(S)
}


paiw<-function(Y0,t,delta){
  #this function aims to compute the \int_0^t Sc(u)du
  #Y0 is the censored data
  paiw = NULL
  ind = order(Y0)
  Y1 = Y0[ind]
  delta1 = delta[ind]
  SC = eKM(Y1,Y1,delta1)
  l = length(Y1)
  for(i in 1:l){
    ss = Y1[Y1<t[i]]
    temp = SC[Y1<t[i]]
    h = diff(c(ss,t[i]))
    paiw[i] = sum(temp*h)
  }
  return(paiw)
}
