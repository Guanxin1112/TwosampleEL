#this program aims to calculate the weight2 of adjusted kafang distri.
#
# Inputs:
# theta -- the parameter
# Delta -- the parameter
# Rd    -- the observed data
# prob  -- the probabltity of quantile
# h     -- the bandwidth

# Outputs:
# Var_weight2 -- the weight of adjusted kafang distri.


Var_weight2<-function(theta,Delta,Rd,prob,h){
  d = Rd$d
  m = length(d)
  
  z = G2.fun(theta,Delta,Rd,prob,h)
  sigma02 = mean(z^2)
  sigma12 = Sigma12.fun(theta,Delta,Rd,prob,h)
  Inform = (sum(d)/theta^2)/m
  beta02 = mean(dG2.fun(theta,Delta,Rd,prob,h))
  
  vv1 = beta02^2/Inform 
  vv2 = sigma02 + vv1
  vv3 = sigma12 + vv1
  var_weight2 = vv2/vv3
  return(var_weight2)
}


Sigma12.fun<-function(theta,Delta,Rd,prob,h){
  X = Rd$X
  A = Rd$A
  delta = Rd$delta
  n = length(X)
  tilV = X-A
  WW1 = paiw(tilV,tilV,delta)
  WW2 = paiw(tilV,X,delta)
  
  temp10 = outer(X,tilV,'-')
  temp20 = outer(tilV,tilV,'-')
  temp30 = outer(1/WW2,WW1,'*')
  temp1 = (temp10>=0)*1
  temp2 = (temp20>=0)*1
  temp3 = (1 - temp30)
  g2 = temp3 * temp1
  
  part0 = G2.fun(theta,Delta,Rd,prob,h)
  
  B0 = part0 * g2
  B2 = colSums(B0)/n
  phi = colSums(temp2)/n
  part1 = B2 * (1-delta) /phi
  
  ybar = colSums(temp2)
  part20 = part1/ybar
  part21 = part20 * t(temp2)
  part2 = colSums(part21)
  
  ss = part0 + part1 - part2
  sigma12 = mean(ss^2)
  return(sigma12)
}

Gamma02.fun<-function(theta,Delta,Rd,prob,h){
  X = Rd$X
  A = Rd$A
  delta = Rd$delta
  n = length(X)
  tilV = X-A
  W2 = paiw(tilV,X,delta)
  a = delta / W2
  
  eta = qexp(prob,rate = theta)
  ss = rep(eta+Delta,n)-X 
  phi = dnorm(ss/h)/h
  
  r02 = mean(a*phi)
  return(r02)
}


Sigma2.fun<-function(theta,Delta,Rd,prob,h){
  # y = Rd[,4]
  d = Rd$d
  m = length(d)
  # theta_MLE = sum(d)/sum(y)
  sigma12 = Sigma12.fun(theta,Delta,Rd,prob,h)
  Inform = (sum(d)/theta^2)/m
  beta02 = mean(dG2.fun(theta,Delta,Rd,prob,h))
  gamma02 = Gamma02.fun(theta,Delta,Rd,prob,h)
  
  sigma = (sigma12 + beta02^2/Inform)/(gamma02^2)
  return(sigma)
}

