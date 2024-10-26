#two_emplkdown1 is Newton hill-down method using relmo to calculate lambda and 
#estimate theta_EL,Delta_EL at the same time
two_emplkdown2<-function(beta,Rd,prob,h){
  svdtol = 1e-9
  thetap = beta[1]
  Deltap = beta[2]
  y = Rd$Y
  d = Rd$d
  #likhd0 = -1e7
  step = 0
  max_step = 200
  gradient_wts = c(1.0, 0.37267800, 0.22395160, 0.11938701, 0.06250000, 0.03125000, 0.01562500,
                   0.00781250, 0.00390625, 0.00195313, 0.00097656, 0.00004883, 0.00000244, 0.00000012, 0.00000001, 0.0)
  
  while(step < max_step){
    step = step + 1
    theta = beta[1]
    Delta = beta[2]
    z =  G2.fun(theta,Delta,Rd,prob,h)
    EL = eqelm(z,maxit=25, gradtol=1e-7,svdtol = 1e-9, itertrace=F)
    nlam = EL$lam
    temp = 1 + nlam * z
    dlog = -llog(temp,1/n)
    likhd0 = sum(dlog) + sum(d)*log(theta) - theta*sum(y)
    
    NWT = Two_Grad2(beta,nlam,Rd,prob,h)
    grad = NWT$grad
    hess = NWT$hess
    
    svdh = svd(hess)
    if( min(svdh$d) < max(svdh$d)*svdtol + 1e-128)
      svdh$d <- svdh$d + max(svdh$d)*svdtol + 1e-128
    ##make it unsingular
    shess = svdh$v %*% diag(1/svdh$d) %*% t(svdh$u)
    
    ninner = 0
    for(j in 1:length(gradient_wts)){
      betanew = beta - gradient_wts[j]*shess %*% grad
      thetanew = betanew[1]
      if(thetanew<=0) {thetanew=thetap;Deltanew=Deltap;likhd1=likhd0;nlamnew = nlam;break} 
      Deltanew = betanew[2]
      znew =  G2.fun(thetanew,Deltanew,Rd,prob,h)
      ELnew = eqelm(znew,maxit=25, gradtol=1e-7,svdtol = 1e-9, itertrace=F)
      nlamnew = ELnew$lam
      tempnew = 1 + nlamnew * znew
      dlognew = -llog(tempnew,1/n)
      likhd1 = sum(dlognew) + sum(d)*log(thetanew) - thetanew*sum(y)
      ninner = ninner+1
      if(likhd1 > likhd0) break
    }
    #if(ninner == length(gradient_wts)) break
    
    #cat("ninner=",ninner,"\n")
    if(likhd1 < likhd0) break  ##exclude the situation that 'ninner == length(gradient_wts)'
    else{
      theta = thetanew
      Delta = Deltanew
      beta = c(thetanew,Deltanew);
      nlam = nlamnew;
      if(abs((likhd0-likhd1)/likhd0) <1e-8) break
      #likhd0 = likhd1;
    }
  }
  #cat("ninner=",ninner,"\n")
  #cat("step=",step,"\n")
  return(cbind(likhd0,theta,Delta,nlam,step,ninner)) #not thetanew,nlamnew
}

Two_Grad2<-function(beta,lambda,Rd,prob,h){
  theta = beta[1]
  Delta = beta[2]
  grad = NULL
  hess = matrix(10000,2,2)
  
  X = Rd$X
  A = Rd$A
  delta = Rd$delta
  n = length(X)
  
  y = Rd$Y
  d = Rd$d
  W2 = paiw(X-A,X,delta)
  a = delta / W2
  
  g1 = G2.fun(theta,Delta,Rd,prob,h)
  dg1 = dG2.fun(theta,Delta,Rd,prob,h)   #\partial G / \partial theta
  ddg1 = ddG2.fun(theta,Delta,Rd,prob,h)  #\partial^2 G / \partial theta^2
  
  eta = qexp(prob,rate = theta)
  ss = rep(eta+Delta,n)-X
  d0 = dnorm(ss/h)/h
  d1 = ddnorm(ss/h)/(h^2)
  
  temp1 = d0 * a         #\partial G / \partial Delta
  temp2 = d1 * a          #\partial^2 G / \partial Delta^2
  temp3 = -temp2 * eta/theta          #\partial^2 G / \partial Delta \partial theta
  
  arg = rep(1,n) + lambda*g1
  Ma =  t(t(arg))
  b0 = apply(Ma,2,llogp,1/n)
  c0 = apply(Ma,2,llogpp,1/n)
  
  b1 = b0 * dg1
  grad[1] = -sum(b1)*lambda + sum(d)/theta-sum(y)
  hess[1,1] = -sum(c0*dg1^2)*lambda^2-sum(b0*ddg1)*lambda-sum(d)/(theta^2)
  
  grad[2] = -sum(b0*temp1)*lambda 
  hess[2,2] = -sum(c0*temp1^2)*lambda^2-sum(b0*temp2)*lambda
  
  hess[1,2] = -sum(temp3*b0)*lambda - sum(c0*dg1*temp1)*lambda^2
  hess[2,1] = hess[1,2]
  list(hess=hess,grad=grad)
}


