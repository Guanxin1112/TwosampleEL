#
#  The function eqelm() calculates the empirical likelihood
#  ratio for Lagrange multiplier lam
#

#  INPUTS:
#     x              Matrix (or vector) of data
#                     One row (or element) per observation
#     lam          Optional starting value for Lagrange multiplier
#     maxit        Optional maximum number of iterations
#     gradtol      Optional tolerance for convergence test
#     svdtol       Optional tolerance for detecting singularity
#                      while solving equations
#     itertrace   Optional flag to print results during iteration
#


#  OUTPUTS:
#     logelr      Log empirical likelihood
#     lambda      Lagrange multiplier
#     grad        Gradient of log likelihood
#     hess        Hessian of log likelihood
#     wts         Relative observation weights at the solution
#     nits        Number of iterations used

eqelm<-function(x, maxit=25, gradtol=1e-7, svdtol = 1e-9, itertrace=F){
  x <-as.matrix(x)
  n <- nrow(x)
  p <-ncol(x)
  z <- x
  lam <- rep(0,p)

  #
  # Preset the weights for combining Newton and gradient
  # steps at each of 16 inner iterations, starting with
  # the Newton step and progressing towards shorter vectors
  # in the gradient direction.  Most commonly only the Newton
  # step is actually taken, though occasional step reductions
  # do occur.
  #
  
  nwts <- c( 3^-c(0:3), rep(0,12) )
  gwts <- 2^( -c(0:(length(nwts)-1)))
  gwts <- (gwts^2 - nwts^2)^.5
  gwts[12:16] <- gwts[12:16] * 10^(-c(1:5))
  
  #
  #    Iterate, finding the Newton and gradient steps, and
  #    choosing a step that reduces the objective if possible.
  #
  
  nits <- 0
  gsize <- gradtol + 1
  while( nits<maxit && gsize > gradtol){
    arg  <- 1 + z %*% lam
    wts1 <- as.vector(llogp(arg, 1/n) )
    wts2 <- as.vector( -llogpp(arg, 1/n) )^.5
    grad <- as.matrix( -z*wts1 )
    grad <- as.vector( apply( grad, 2, sum ) )
    gsize<- mean( abs(grad))
    hess <- z*wts2
    #                                   -1
    #    The Newton step is -(hess'*hess)  grad,
    #  where the matrix hess is a sqrt of the Hessian.
    #  Use svd on hess to get a stable solution.
    #
    
    svdh <-svd(hess)
    #it is more  simple than hess'*hess;
    if( min(svdh$d) < max(svdh$d)*svdtol + 1e-128)
      svdh$d <- svdh$d + max(svdh$d)*svdtol + 1e-128
    #make it unsingular
    nstep <- svdh$v %*% t(svdh$u/svdh$d)
    nstep <- as.vector( nstep %*% matrix(wts1/wts2,n,1) )
    gstep <- -grad
    if (sum(nstep^2) < sum(gstep^2))
      gstep <- gstep*sum(nstep^2)^.5/sum(gstep^2)^.5
    ologelr <- -sum(llog(arg,1/n))
    ninner <- 0
    for(  i in 1:length(nwts) ){
      nlam <- lam+nwts[i]*nstep+gwts[i]*gstep
      nlogelr <- -sum(llog(1+z %*% nlam,1/n))
      if(nlogelr < ologelr){
        lam <- nlam
        ninner <- i
        break
      }
    }
    nits <- nits+1
    if(  ninner==0  ) nits <- maxit
    if( itertrace )
      print( c(lam, nlogelr, gsize, ninner) )
  }
  list( logelr=nlogelr, lam = lam, grad=grad, hess=t(hess)%*%hess, wts=wts1, nits=nits )
}

###############################################################
#
#    The function llog() is equal to the natural
#  logarithm on the interval from eps >0 to infinity.
#  Between -infinity and eps, llog() is a quadratic.
#  llogp() and llogpp() are the first two derivatives
#  of llog().  All three functions are continuous
#  across the "knot" at eps.
#
#    A variation with a second knot at a large value
#  M did not appear to work as well.
#
#    The cutoff point, eps, is usually 1/n, where n
#  is the number of observations.  Unless n is extraordinarily
#  large, dividing by eps is not expected to cause numerical
#  difficulty.
#
###############################################################
llog<-function(z,eps){
  ans<-z
  lo <-(z<eps)
  ans[ lo] <- log(eps) - 1.5 + 2*z[lo]/eps - 0.5*(z[lo]/eps)^2
  # using taylor expansion :ln(1+x)=x-.5*x^2+o(x^2),this is a technique to submit complexity
  ans[ !lo] <- log( z[!lo] ) 
  return (ans)
}
####################################################################
llogp<-function( z, eps ){
  ans <- z
  lo <-(z<eps)
  ans[ lo ] <- 2.0/eps - z[lo]/eps^2
  ans[ !lo ] <- 1/z[!lo]
  return (ans)
}
###############################################################
llogpp<-function( z, eps ){
  ans <- z
  lo <- (z<eps)
  ans[ lo] <- -1.0/eps^2
  ans[ !lo ] <- -1.0/z[!lo]^2
  return (ans)
}
