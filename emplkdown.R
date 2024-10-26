#emplkdown1 is Newton down-hill method using relmo to calculate lambda and 
#use the true value of Delta0 to estimate theta_EL 
#  by maximaize empirical likelihood using emplk1 function
#then estimate Delta using estimating equation G1 with theta_EL 

emplkdown1<-function(theta,Delta,Rd,prob,h){
  theta0 = theta
  #likhd0 = -1e7
  step = 0
  max_step = 100
  y = Rd$Y
  d = Rd$d
  z =  G1.fun(theta,Delta,Rd,prob,h)
  EL = eqelm(z,maxit=25, gradtol=1e-7,svdtol = 1e-9, itertrace=F)
  nlam = EL$lam
  gradient_wts = c(1.0, 0.37267800, 0.22395160, 0.11938701, 0.06250000, 0.03125000, 0.01562500,
    0.00781250, 0.00390625, 0.00195313, 0.00097656, 0.00004883, 0.00000244, 0.00000012, 0.00000001, 0.0)
  
  while(step < max_step){
    step = step + 1
    z =  G1.fun(theta,Delta,Rd,prob,h)
    EL = eqelm(z,maxit=25, gradtol=1e-7,svdtol = 1e-9, itertrace=F)
    nlam = EL$lam
    temp = 1 + nlam * z
    dlog = -llog(temp,1/n)
    likhd0 = sum(dlog) + sum(d)*log(theta) - theta*sum(y)
    
    NWT = DDtheta1(theta,Delta,nlam,Rd,prob,h)
    grad = NWT[1]
    hess = NWT[2]
    ninner = 0
    for(i in 1:length(gradient_wts)){
      thetanew = theta - gradient_wts[i]*grad/hess
      if(thetanew<=0) {thetanew=theta0;likhd1=likhd0;nlamnew = nlam;break} #thetanew=theta;likhd0=likhd0;?
      znew =  G1.fun(thetanew,Delta,Rd,prob,h)
      ELnew = eqelm(znew,maxit=25, gradtol=1e-7,svdtol = 1e-9, itertrace=F)
      nlamnew = ELnew$lam
      tempnew = 1 + nlamnew * znew
      dlognew = -llog(tempnew,1/n)
      likhd1 = sum(dlognew) + sum(d)*log(thetanew) - thetanew*sum(y)
      ninner = ninner+1
      if(likhd1 > likhd0) break
      #if(ninner == length(gradient_wts)) break
    }
    #cat("ninner=",ninner,"\n")

    #if(ninner == length(gradient_wts)) break; #leave out the situation
    #######when ninner=16 but likh1>=lihd0. So we change to next row
    
    if(likhd1 < likhd0) break   #exclude the situation that when ninner = 1 to 
    ########length(gradient_wts), likhd1 always less than likhd0
    else{
      theta = thetanew;
      nlam = nlamnew;
      if(abs((likhd0-likhd1)/likhd0) <1e-8) break
      likhd0 = likhd1; # unnecessary since we have updated theta
      }
  }
  #cat("thetanew=",thetanew,"\n")
  return(cbind(likhd0,theta,nlam,step)) #not thetanew,nlamnew since new is not good than old
}


##################################################################
#estimating function for Delta=Delta0
G01.fun<-function(theta,Delta,Rd,prob,h){
  sum(G1.fun(theta,Delta=Delta0,Rd,prob,h))
} 

#hession matrix
DDtheta1<-function(theta,Delta,lambda,Rd,prob,h){
  n = length(Rd$X)
  y = Rd$Y
  d = Rd$d
  g1 = G1.fun(theta,Delta,Rd,prob,h)
  dg1 = dG1.fun(theta,Delta,Rd,prob,h)
  ddg1 = ddG1.fun(theta,Delta,Rd,prob,h)
  
  arg = rep(1,n) + lambda*g1
  Ma =  t(t(arg))
  b0 = apply(Ma,2,llogp,1/n)
  c0 = apply(Ma,2,llogpp,1/n)
  b1 = b0 * dg1
  
  grad = -sum(b1)*lambda + sum(d)/theta-sum(y)
  hess = -sum(c0*dg1^2)*lambda^2-sum(b0*ddg1)*lambda-sum(d)/(theta^2)
  return(cbind(grad,hess))
}


G1.fun<-function(theta,Delta,Rd,prob,h){
  X = Rd$X
  A = Rd$A
  delta = Rd$delta
  n = length(X)
  
  SC = eKM(X-A,X-A,delta)
  W = X * SC
  a = delta / W
  #ty = Rd[,6]
  #ordy = t(t(sort(ty)))
  #dfy = apply(ordy,2,Dist_y,theta)
  #eta = fquant(dfy,prob,ordy)
  eta = qexp(prob,rate = theta)
  ss = rep(eta+Delta,n)-X 
  phi = pnorm(ss/h)-rep(prob,n)
  g1 = phi * a
  return(g1)
} 

dG1.fun<-function(theta,Delta,Rd,prob,h){
  X = Rd$X
  A = Rd$A
  delta = Rd$delta
  n = length(X)
  #Y = Rd[,4]
  SC = eKM(X-A,X-A,delta)
  W = X * SC
  a = delta / W
  #ty = Rd[,6]
  #ordy = t(t(sort(ty)))
  #dfy = apply(ordy,2,Dist_y,theta)
  #eta = fquant(dfy,prob,ordy)
  eta = qexp(prob,rate = theta)
  feta = dexp(eta,rate = theta)
  temp = -eta/theta
  
  ss = rep(eta+Delta,n)-X
  dphi = temp*(dnorm(ss/h)/h-rep(feta,n))
  
  dg1 = dphi * a
  return(dg1)
}

ddnorm<-function(y) -y*dnorm(y)

ddG1.fun<-function(theta,Delta,Rd,prob,h){
  X = Rd$X
  A = Rd$A
  delta = Rd$delta
  n = length(X)
  #Y = Rd[,4]
  SC = eKM(X-A,X-A,delta)
  W = X * SC
  a = delta/W

  eta = qexp(prob,rate = theta)
  feta = dexp(eta,rate = theta)
  
  temp1 = -eta/theta
  temp2 = 2*eta/(theta^2)
  temp3 = theta*eta*exp(-theta*eta)
  
  ss = rep(eta+Delta,n)-X
  d1 = (temp1^2)*ddnorm(ss/h)/(h^2)
  d2 = temp2*dnorm(ss/h)/h
  d3 = temp3*temp1
  d4 =  feta*temp2 
  
  ddg1 = (d1+d2-d3-d4)*a
  return(ddg1)
}
