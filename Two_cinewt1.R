#this program aims to calculate the coverage interval for EL 
#it is for theorem 3
#with estimation equation G1 i.e. weight W1=XSc(X-A)
Two_cinewt1<-function(Delta,thetahat,Deltael,thetael,Rd,prob,h){
  ## Input:
  ##   thetahat -- theta eatimation under H0,i.e. given the ture Delta0
  ##   Deltael,thetael  -- Delta and theta estimation under H1
  y = Rd$Y
  d = Rd$d
  z0 = G1.fun(thetahat,Delta,Rd,prob,h)
  EL0 = eqelm(z0,maxit=25, gradtol=1e-7,svdtol = 1e-9, itertrace=F)
  H0_logelr = EL0$logelr + sum(d)*log(thetahat) - thetahat*sum(y)
  
  z1 = G1.fun(thetael,Deltael,Rd,prob,h)
  EL1 = eqelm(z1,maxit=25, gradtol=1e-7,svdtol = 1e-9, itertrace=F)
  H1_logelr = EL1$logelr + sum(d)*log(thetael) - thetael*sum(y)
  
  weight = Var_weight1(thetael,Deltael,Rd,prob,h)
  adjelr = 2 * (H1_logelr-H0_logelr)*weight
  return(adjelr-chia)
}

