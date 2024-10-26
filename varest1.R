#this program aims to calculate the weight1 of adjusted kafang distri.
#
# Inputs:
# theta -- the parameter
# Delta -- the parameter
# Rd    -- the observed data
# prob  -- the probabltity of quantile
# h     -- the bandwidth
    
# Outputs:
# Var_weight1 -- the weight of adjusted kafang distri.


Var_weight1<-function(theta,Delta,Rd,prob,h){
  d = Rd$d
  m = length(d)
  
  z = G1.fun(theta,Delta,Rd,prob,h)
  sigma01 = mean(z^2)
  sigma11 = Sigma11.fun(theta,Delta,Rd,prob,h)
  Inform = (sum(d)/theta^2)/m
  beta01 = mean(dG1.fun(theta,Delta,Rd,prob,h))
  
  vv1 = beta01^2/Inform 
  vv2 = sigma01 + vv1
  vv3 = sigma11 + vv1
  var_weight1 = vv2/vv3
  return(var_weight1)
}

#Var_weight1(theta,Delta=Delta0,Rd=rd,prob=p,h=h)

Sigma11.fun<-function(theta,Delta,Rd,prob,h){
  X = Rd$X
  A = Rd$A
  delta = Rd$delta
  n = length(X)
  tilV = X-A
 
  temp10 = outer(tilV,tilV,'-')   #outer(X,tilV,'-')
  temp20 = outer(tilV,tilV,'-')
  temp1 = (temp10>=0)*1
  temp2 = (temp20>=0)*1
  
  part0 = G1.fun(theta,Delta,Rd,prob,h)

  B0 = part0 * temp1
  B1 = colSums(B0)/n
  phi = colSums(temp2)/n
  part1 = B1 * (1-delta) /phi
  
  ybar = colSums(temp2)
  part20 = part1/ybar
  part21 = part20 * t(temp2)
  part2 = colSums(part21)
  
  ss = part0 + part1 - part2
  sigma11 = mean(ss^2)
  return(sigma11)
}

Gamma01.fun<-function(theta,Delta,Rd,prob,h){
  X = Rd$X
  A = Rd$A
  delta = Rd$delta
  n = length(X)
  tilV = X-A
  SC = eKM(X-A,X-A,delta)
  W = X * SC
  a = delta / W
  
  eta = qexp(prob,rate = theta)
  ss = rep(eta+Delta,n)-X 
  phi = dnorm(ss/h)/h
  
  r01 = mean(a*phi)
  return(r01)
}


Sigma1.fun<-function(theta,Delta,Rd,prob,h){
  d = Rd$d
  m = length(d)
  # theta_MLE = sum(d)/sum(y)
  sigma11 = Sigma11.fun(theta,Delta,Rd,prob,h)
  Inform = (sum(d)/theta^2)/m
  beta01 = mean(dG1.fun(theta,Delta,Rd,prob,h))
  gamma01 = Gamma01.fun(theta,Delta,Rd,prob,h)
  
  sigma = (sigma11 + beta01^2/Inform)/(gamma01^2)
  return(sigma)
}



BGMM1<-function(Delta,delta,eta_mle,X,W,h){
  a = delta / W
  ss = rep(eta_mle+Delta,n)-X
  phi = pnorm(ss/h)-rep(p,n)
  
  g1 = sum(phi * a)
  return(g1)
} 

