blksmd_asian_pca <- function( n = 12 ){
  t.exp <- seq(1/n,1,1/n)
  t.mat <- matrix(t.exp,n,n)
  covar <- pmin(t.mat,t(t.mat))
  var <- t.exp*t.exp
  
  svd<-svd(covar)
  svd$d<-sqrt(svd$d)
  #svd$u <- svd$u * sign(sum(svd$u[,1]))
  svd$v <- NULL
  covar.sqrt <- svd$u %*% diag(svd$d)
  
  grad <- as.vector( t(covar.sqrt) %*% rep(1,n) )
  Q <- blksmd_householder(grad)

  covar.sqrt.q <- covar.sqrt %*% Q   # sigma[,1] == fac1st
  fac1st <- covar.sqrt.q[,1]
  fac1st.norm <- sqrt( as.vector( crossprod( fac1st ) ) )
  svdN <- svd(covar.sqrt.q[,-1])
  svdN$v <- NULL
  svdN$u <- cbind(fac1st/fac1st.norm, svdN$u)
  svdN$d <- c(fac1st.norm, svdN$d)

  crossp1 <- sum(fac1st)/sqrt(n)
  len_scale <- svdN$d / crossp1

  rv<-list()
  rv$cov<-covar
  rv$t.exp <- t.exp
  rv$d <- svdN$d
  rv$u1 <- fac1st
  rv$u2 <- svdN$d[2] * svdN$u[,2]
  rv$u3 <- svdN$d[3] * svdN$u[,3]
  rv$u4 <- svdN$d[4] * svdN$u[,4]
  rv$u1cont <- sqrt(3)*(rv$t.exp-0.5*rv$t.exp*rv$t.exp)
  rv$s1 <- sqrt(2)/(pi*(1-0.5))*sin((1-0.5)*pi*rv$t.exp )
  rv$s2 <- sqrt(2)/(pi*(2-0.5))*sin((2-0.5)*pi*rv$t.exp )
  rv$s3 <- sqrt(2)/(pi*(3-0.5))*sin((3-0.5)*pi*rv$t.exp )
  rv$vartotal <- sum(var)
  rv$len_scale <- len_scale
  
  return(rv)
}

#' Main pricing function for spread and basket options
#' 
#' @param strk strike prices (array)
#' @param spot spot prices (array)
#' @param t.exp time to expiry
#' @param vol volatilities (array)
#' @param corr correlation matrix between assets
#' @param r interest rate
#' @param d dividend rate
#' @param lambda parameter to determin quadrature size
#' @param n.quad if given, ignore lambda and force the given quadrature sizes
#' @param uniform if TRUE, use equally-distanced dense grid instead of quadrature.
#' Computationally heavy. Used only for validation purpose. 
#' @param detail return more information if TRUE
#' @param CV use control-variate using deta if TRUE
#' @param callput call/put for TRUE/FALSE
#' 
#' @return option prices for given strike prices (array)
#' 
#' @references Choi, J. (2018). Sum of all Black-Scholes-Merton models: An efficient 
#' pricing method for spread, basket, and Asian options. Journal of Futures Markets, 
#' 38(6), 627–644. https://doi.org/10.1002/fut.21909
blksmd_basket <- function( strk, spot, t.exp, vol, wts=1/length(vol), corr = 0, 
                           r = 0, d = 0, lambda = 3, n.quad = NA, uniform = F, detail = F, CV = T, callput = T ){
  
  n.var <- length(wts)
  fwd.wts <- exp((r-d)*t.exp)*spot*wts
  fwd.wts.unit <- fwd.wts / sqrt(sum(fwd.wts*fwd.wts))
  
  covar <- vol %o% vol * corr * t.exp
  var <- vol*vol * t.exp
  var.sqrt <- sqrt(var)

  #svd <- svd(covar)
  #svd$d <- sqrt(svd$d)
  #svd$u <- sign(sum(svd$u[,1])) * svd$u
  #svd$v <- NULL
  #covar.sqrt <- svd$u %*% diag(svd$d) 

  covar.sqrt <- t(chol(covar))

  fac1st <- as.vector( covar %*% fwd.wts )
  fac1st <- fac1st / sqrt(sum(fac1st*fwd.wts))
  thres <- 0.01 * var.sqrt
  
  mask.fac.fix <- (sign(fwd.wts)*fac1st < thres)
  fac1st[mask.fac.fix] <- (sign(fwd.wts)*thres)[mask.fac.fix]

  if(any(mask.fac.fix)) {
    grad <- solve(covar.sqrt, fac1st)
  } else {
    grad <- as.vector( t(covar.sqrt) %*% fwd.wts )
  }

  #varN <- var - target * facmat^2
  q.mat <- blksmd_householder(grad)
  covar.sqrt.q <- covar.sqrt %*% q.mat   # sigma[,1] == fac1st, but...

  fac1st <- covar.sqrt.q[,1]
  fac1st.norm <- sqrt(sum(fac1st*fac1st))

  svdN <- svd(covar.sqrt.q[,-1]) # all except the first column
  svdN$v <- NULL
  svdN$d <- c(fac1st.norm, svdN$d)
  svdN$u <- cbind(fac1st/fac1st.norm, svdN$u)
  
  crossp1 <- sum(fwd.wts.unit*fac1st) #crossprod
  len_scale <- svdN$d/crossp1

  if( uniform ) {
    if( CV ) {
      warning('CV turned off in uniform mode')
      CV <- FALSE
    }
    nodes <- seq(-6.5, 6.5, 0.25)
    weights <- dnorm(nodes) * 0.25
    quadSt <- list( nodes = nodes, weights = weights )
    n.quad <- rep( length(nodes), n.var )
  } else if( !is.na( n.quad ) ){
    n.quad <- rep( n.quad, n.var )
  } else {
    n.quad <- pmin( round(len_scale*lambda + 1), 30)  # 1.25
  }
  
  mask <- ( n.quad > 1 ); mask[1] <- TRUE
  n.quad[1] <- 1

  facmat <- svdN$u[,mask] %*% diag(svdN$d[mask], nrow=sum(mask), ncol=sum(mask))
  # In case mask has only 1 TRUE value, diag(...) is still 1x1 matrix as we pass nrow & ncol
  # therefore, facmat is still in matrix form (Nx1). 
  
  #facmat <- cbind(fac1st, svdN$u %*% diag(svdN$d)  should be same but diag(double) fails when n=2
  
  varN <- rowSums(facmat*facmat) - fac1st*fac1st

  idx <- 1:n.var
  idx <- idx[mask]
  
  nodes <- list()

  for (k in idx) {
    if( k == 1L ) {
      ww <- 1 #put a trivial value
      xx <- 0 #put a trivial value
    } else {
      if(!uniform) {
        quadSt <- statmod::gauss.quad.prob(n.quad[k], dist="normal")
      }
      nodes[[k]] <- quadSt$nodes
      
      xx_prev <- xx
      xx <- numeric(0) #empty vector
      for (node in quadSt$nodes) {
        xx <- cbind(xx,rbind(xx_prev,node))
      }
      ww <- as.vector(ww %o% quadSt$weights)
    }
  }
  
  if(length(xx)>1) {
    # length(xx)>1 is the case where no discretization required for the factors j>1
    dimnames(xx) <- list(NULL,NULL)
  }
  
  is.fac.uniform <- ( max(fac1st)-min(fac1st) < .Machine$double.eps * 100 )
  fac1st.uniform.val <- mean(fac1st)
  
  # wrap as.matrix() to avoid sigma[,mask] becomes a vector when mask has only one component
  f_k.mat <- exp( facmat %*% xx - 0.5*varN )
  fwd.mat <- fwd.wts * f_k.mat
  fac1st.exp <- exp( -0.5*fac1st*fac1st )
  
  # statistics on fwd samples
  f_k.err <- f_k.mat %*% ww - 1
  
  zz <- rep(NA,length(ww))
  price.fwd <- rep(NA,length(strk))
  
  for (k in 1:length(strk)){
    price <- rep(0,length(ww))
    roots <- rep(0,length(ww))
    
    for (j in 1:length(ww)){

      if( is.fac.uniform ) {
        root <- log(strk[k]/sum(fwd.mat[,j]))/fac1st.uniform.val + 0.5*fac1st.uniform.val
      } else {
        root <- root_one(fac1st.exp*fwd.mat[,j], fac1st, strk[k])
      }

      price[j] <- blks_N(fwd.mat[,j], fac1st, strk[k], root, callput = callput ) # call: TRUE, put: FALSE
      roots[j] <- root
    }
    
    price.fwd[k] <- sum(price*ww)
    if( CV ) {
      delta <- as.vector( ( fwd.mat * pnorm(outer(fac1st, -roots, "+")) ) %*% ww )  # Delta_k * F_k
      if(!callput) {
        delta <- delta - 1
      }
      price.fwd[k] <- price.fwd[k] - sum(delta * f_k.err)
    }
  }
  
  ret <- list(call=price.fwd * exp(-r*t.exp))
  
  if( detail ){
    ret$dim <- c( n.quad, prod(n.quad,na.rm=TRUE) )
    ret$facmat <- cbind(c(sum(fwd.wts.unit*fac1st), fwd.wts.unit),
                        rbind( svdN$d[mask], facmat ), c(sqrt(sum(var)), sqrt(var)))
    ret$len_scale <- len_scale
    ret$f_k.err <- f_k.err
    ret$zz <- array(zz, n.quad)
  }
  return(ret)
}

#' Main pricing function for Asian options
#' 
#' @param strk strike prices (array)
#' @param spot spot prices (array)
#' @param t.obs observation time (array). The last value is the expiry.
#' @param vol volatility
#' @param wts weights
#' @param r interest rate
#' @param d dividend rate
#' @param CV use control-variate using deta if TRUE
#' @param callput call/put for TRUE/FALSE
#' @param detail return more information if TRUE
#' 
#' @return option prices for given strike prices (array)
#' 
#' @references Choi, J. (2018). Sum of all Black-Scholes-Merton models: An efficient 
#' pricing method for spread, basket, and Asian options. Journal of Futures Markets, 
#' 38(6), 627–644. https://doi.org/10.1002/fut.21909
blksmd_asian <- function( strk, spot, t.obs, vol, wts=1/length(t.obs), r = 0, d = 0, 
                          CV = T, callput = T, detail = F ){
  
  n.obs <- length(t.obs)
  fwd.wts <- exp((r-d)*t.obs)*spot*wts
  fwd.wts.unit <- fwd.wts / sqrt(sum(fwd.wts*fwd.wts))
  
  t.mat <- matrix(t.obs, ncol=n.obs, nrow=n.obs)
  covar <- pmin(t.mat, t(t.mat))*(vol*vol)
  var <- vol^2 * t.obs

  t.diff.sqrt <- sqrt( c( t.obs[1], diff(t.obs) ) )
  covar.sqrt <- matrix(t.diff.sqrt*vol, ncol=n.obs, nrow=n.obs, byrow = T)
  covar.sqrt[ upper.tri(covar.sqrt) ] <- 0

  n.quad <- rep(1, n.obs)
  n.quad[1:8] <- c(99, 3, 3, 3, 3, 1, 1, 1) # 81 evaluations
  ## Anternative candidate is 
  #n.quad[1:6] <- c(99, 3, 3, 3, 3, 3) # 108 evaluations
  mask <- ( n.quad > 1 )
  n.quad[1] <- 1

  grad <- as.vector( t(covar.sqrt) %*% fwd.wts )
  Q <- blksmd_householder(grad)

  covar.sqrt.q <- covar.sqrt %*% Q   # sigma[,1] == fac1st
  fac1st <- covar.sqrt.q[,1]
  fac1st.norm <- sqrt( sum( fac1st*fac1st ) )

  #### because we need only a few first singular values and vectors,
  #### better use irlba package. R's default svd computes all singular values.
  #### svd <- svd(covar.sqrt.q[,-1], nu=sum(mask)-1, nv=0) # in case we use the default svd
  svd <- irlba::irlba(covar.sqrt.q[,-1], nu=sum(mask)-1, nv=0) # all except the first column
  svd$u <- cbind(fac1st/fac1st.norm, svd$u)
  svd$d <- c(fac1st.norm, svd$d)

  crossp1 <- sum(fwd.wts.unit*fac1st)
  len_scale <- svd$d / crossp1
  
  facmat <- svd$u %*% diag(svd$d[mask], nrow=sum(mask), ncol=sum(mask)) 
  #facmat <- cbind(fac1st, svdN$u %*% diag(svdN$d)  should be same but diag(double) fails when n=2
  
  varN <- rowSums(facmat*facmat) - fac1st*fac1st
  
  idx <- 1:n.obs
  idx <- idx[mask]

  for (k in idx) {
    if(k==1L){
      ww <- 1 #put a trivial value
      xx <- 0 #put a trivial value
    } else {
      quadSt <- statmod::gauss.quad.prob(n.quad[k], dist='normal')

      xx_prev <- xx
      xx <- numeric(0)
      for (node in quadSt$nodes) {
        xx <- cbind(xx,rbind(xx_prev,node))
      }
      ww <- as.vector(ww %o% quadSt$weights)
    }
  }
  
  dimnames(xx) <- list(NULL,NULL)
  
  # wrap as.matrix() to avoid sigma[,mask] becomes a vector when mask has only one component
  f_k.mat <- exp( facmat %*% xx - 0.5*varN )
  fwd.mat <- fwd.wts * f_k.mat
  fac1st.exp <- exp( -0.5*fac1st*fac1st )
  
  # statistics on fwd samples
  f_k.err <- f_k.mat %*% ww - 1

  price.fwd <- rep(NA, length(strk))
  
  for (k in 1:length(strk)){
    price <- rep(0,length(ww))
    roots <- rep(0,length(ww))
    
    for (j in 1:length(ww)){
      root <- root_one(fac1st.exp*fwd.mat[,j], fac1st, strk[k])
      roots[j] <- root
      price[j] <- blks_N(fwd.mat[,j], fac1st, strk[k], root, callput = callput)
    }
    
    price.fwd[k] <- sum(price*ww)
    
    if(CV) {
      delta <- as.vector( ( fwd.mat * pnorm(outer(fac1st, -roots, "+")) ) %*% ww )  # Delta_k * F_k
      if(!callput) {
        delta <- delta - 1
      }
      price.fwd[k] <- price.fwd[k] - sum(delta * f_k.err)
    }
  }

  ret <- list(call=price.fwd * exp(-r*max(t.obs)))
  
  if( detail ){
    ret$dim <- c(n.quad[mask], prod(n.quad[mask], na.rm=TRUE))
    ret$len_scale <- len_scale
    ret$d <- svd$d
    ret$u <- svd$u
    ret$vartotal <- sum(var)
    ret$missingvol <- 1-sqrt(sum(varN)/sum(var))
    ret$f_k.err <- f_k.err
  }
  return(ret)
}


#' Margrabe's exchange option formula (Spread Opt with K=0)
#'
#' @param spot two spot prices as an array
#' @param t.exp time to expiry
#' @param vol volatility
#' @param rho correlation
#' @param d dividend rate
blksmd_margrabe <- function( spot, t.exp, vol, rho, d = 0 ){
  # S1 - S2 >= 0
  std <- sqrt(sum(vol*vol)-2*rho*prod(vol))*sqrt(t.exp)
  fwd <- spot * exp(-d*t.exp)
  d_p <- log(fwd[1]/fwd[2])/std + 0.5*std
  d_m <- d_p - std
  call <- fwd[1]*pnorm(d_p)-fwd[2]*pnorm(d_m)
  return(call)
}

#' Spread option pricing method by Bjerksund & Stensland (2014)
#' Implemented for performance comparison. 
#' 
#' @references Bjerksund, P., & Stensland, G. (2014). Closed form spread option valuation. 
#' Quantitative Finance, 14(10), 1785–1794.
blksmd_bjerkspread <- function( strk, spot, t.exp, vol, rho, r = 0, d = 0 ){
  call <- rep(NA, length(strk))
  fwd <- spot * exp((r-d)*t.exp)
  
  for(k in 1:length(strk)){
    a<-fwd[2]+strk[k]
    b<-fwd[2]/a
    
    std <- vol*sqrt(t.exp)

    std11 <- std[1]*std[1]
    std12 <- std[1]*std[2]
    std22 <- std[2]*std[2]
    
    std0 <- sqrt(std11 - 2*b*rho*std12 + b^2*std22)
    
    d1 <- ( log(fwd[1]/a) + 0.5*std11 - b*rho*std12 + 0.5*b^2*std22 ) / std0
    d2 <- ( log(fwd[1]/a) - 0.5*std11 + rho*std12 + (0.5*b^2-b)*std22 ) / std0
    d3 <- ( log(fwd[1]/a) - 0.5*std11 + 0.5*b^2*std22 ) / std0
    call[k] <- (fwd[1]*pnorm(d1)-fwd[2]*pnorm(d2)-strk[k]*pnorm(d3))
  }
  return(exp(-r*t.exp)*call)
}

#' Spread option pricing method by Lo. (2015) Strang's splitting approximation I
#' Implemented for performance comparison. 
#' 
#' @param Kirk if TRUE, return Kirk's spread option formula
#' 
#' @references Lo, C.-F. (2015). Pricing Spread Options by the Operator Splitting 
#' Method. Wilmott, 2015(September), 64–67. https://doi.org/10.1002/wilm.10449
blksmd_splittingpread <- function( strk, spot, t.exp, vol, rho, r = 0, d = c(0,0), Kirk = F ){
  call <- rep(NA, length(strk))
  fwd <- spot * exp((r-d)*t.exp)
  std <- vol*sqrt(t.exp)
  
  for(k in 1:length(strk)){
    ratio <- fwd[2]/(fwd[2]+strk[k])
    std2r <- std[2]*ratio
    std_ <- sqrt(std[1]^2 - 2*rho*std[1]*std2r + std2r^2)

    d1 <- log(fwd[1]/fwd[2]*ratio) / std_ + 0.5*std_
    d2 <- d1 - std_

    call[k] <- fwd[1]*pnorm(d1)-(fwd[2]+strk[k])*pnorm(d2)
    
    if(!Kirk){
      coef <- -0.5*std2r^2*strk[k]*dnorm(d2)
      term1 <- (rho*std[1]-std2r)*std[2]/std_^2
      term2 <- d1*d2 + (1-rho^2)*(std[1]/(rho*std[1]-std2r))^2

      correction <- coef*term1*( d2*(1-rho*std[1]/std2r) - 0.5*std_*term2*term1*(1-ratio) )
      call[k] <- call[k] + correction
    }
  }
  return(exp(-r*t.exp)*call)
}

#' Spread option pricing method by Li Deng Zhou wiht y0 = 0
#' Implemented for performance comparison. 
#' 
#' @references Li, M., Deng, S., & Zhou, J. (2008). Closed-Form Approximations 
#' for Spread Option Prices and Greeks. The Journal of Derivatives, 15(3), 58–80.
blksmd_CalcSpreadOptLDZ <- function( strk, spot, t.exp, vol, rho, r = 0, d = 0 ){
  call <- rep(NA, length(strk))
  fwd <- spot * exp((r-d)*t.exp)
  std <- vol * sqrt(t.exp)
  
  rho.comp <- sqrt(1-rho*rho)*std[1]
  mu1 <- log(fwd[1]) - 0.5*std[1]*std[1]
  mu2.exp <- fwd[2] * exp(-0.5*std[2]*std[2]) # R = exp(mu2) with y0 = 0

  # vector for the rest
  r.plus.k <- mu2.exp + strk # (R+K) in DLZ (vectorized)

  epsilon <- -1/(2*rho.comp) * (std[2]*std[2]*mu2.exp*strk)/(r.plus.k*r.plus.k)
  
  C3 <- ( mu1 - log(r.plus.k) )/rho.comp
  D3 <- ( rho*std[1] - std[2]*mu2.exp/r.plus.k )/rho.comp
  
  C2 <- C3 + std[2]*( D3 + epsilon*std[2] )
  D2 <- D3 + 2*std[2]*epsilon
  
  C1 <- C3 + rho*std[1]*(D3 + epsilon*rho*std[1]) + rho.comp
  D1 <- D3 + 2*rho*std[1]*epsilon

  I_S1 <- blksmd_CalcSpreadOptLDZ_I( C1, D1, epsilon )
  I_S2 <- blksmd_CalcSpreadOptLDZ_I( C2, D2, epsilon )
  I_K  <- blksmd_CalcSpreadOptLDZ_I( C3, D3, epsilon )

  price.fwd <- fwd[1]*I_S1 - fwd[2]*I_S2 - strk*I_K

  return(exp(-r*t.exp)*price.fwd)
}

#' An auxiliry function used by blksmd_CalcSpreadOptLDZ
#' 
blksmd_CalcSpreadOptLDZ_I <- function( u, v, eps ) {
  u2 <- u*u
  u4 <- u2*u2
  v2 <- v*v
  v4 <- v2*v2
  v6 <- v4*v2
  v2.sqrt <- sqrt(1+v2)
  
  arg.norm <- u/v2.sqrt
  J0 <- pnorm(arg.norm)
  J1 <- (1+(1+u2)*v2)/v2.sqrt^5L * dnorm(arg.norm)
  J2 <- (6*(1-u2)*v2 + (21-2*u2-u4)*v4 + 4*(3+u2)*v6 - 3)*u/v2.sqrt^11L * dnorm(arg.norm)
  
  return( J0 + eps*(J1 + 0.5*eps*J2))
}


#' Houseohlder reflection
#' 
#' @param u1 an input vector
blksmd_householder <- function( u1 ) {
  #returns a Householder reflection (orthonormal matrix) which maps e1 into the normalized u1
  v <- u1/sqrt(sum(u1*u1))
  v[1] <- v[1]-1
  
  if( abs(v[1]) < .Machine$double.eps*100 )
    return( diag(length(v)) )
  else {
    R <- diag(length(v)) + v%o%v/v[1]
    return( R )
  }
}

root_basket <- function(f, a, K) {
  # return the interval of x where f*exp(a*x) <= K
  # f >= 0
  
  x.bdd <- c(-7, 7)
  y.vec.bdd <- f*exp(a%o%x.bdd)
  y.bdd <- colSums(y.vec.bdd)

  # we use the fact the f'(x) is increasing function
  if( K<= 0 )
    return( NA )
  else if( y.bdd[1] <= K & y.bdd[2] <= K )
    return(matrix(c(-Inf,Inf),2))
  else if( y.bdd[1] <= K ){
    # yHi > K
    xGuess <- (K - sum(f))/sum(a*f)
    x <- root_guess(f,a,K,xGuess)
    return(matrix(c(-Inf, x),2))
    
  } else if( y.bdd[2] <= K ) {
    # yLo > K
    xGuess <- (K - sum(f))/sum(a*f)
    x <- root_guess(f,a,K,xGuess)
    return(matrix(c(x, Inf),2))
    
  } else {
    #2 root case
    yDerBdd <- colSums(a*y.vec.bdd)
    
    if( yDerBdd[1]>=0 || yDerBdd[2]<=0 )
      return( NA )
    
    xMin <- root_guess(f*a,a,0,0)
    yMin <- sum(f*exp(a*xMin))
    
    if(yMin>=K)
      return( NA ) # min is above K, so no solution
    
    x1 <- root_guess(f,a,K,xLo)
    
    if(x1>xMin)
      stop('unexpected root location')
    
    x2 <- root_guess(f,a,K,xHi)
    
    if(x2<xMin)
      stop('unexpected root location')
    
    return(matrix(c(x1,x2),2))
  }
}

root_one <- function(f, a, K) {
  # return the root of f*exp(a*x) <= K
  # f can be negative, the sign of f and a should be same so f'(x) >= 0 so monotonically increasing
  if( min(f*a) < -.Machine$double.eps*1e4 )
    stop('Invalid argument for root_one')

  x.bdd <- c(-9, 9)
  y.vec.bdd <- f*exp(a%o%x.bdd)
  y.bdd <- colSums(y.vec.bdd)

  # we use the fact the f'(x) >= 0
  if( y.bdd[1] > K )
    return( -Inf )
  else if( y.bdd[2] < K )
    return( Inf )
  else{
    if( min(f)<0 )
      x.guess <- 0 # spread
    else
      x.guess <- min(c((K-sum(f))/sum(f*a), log(K/f)/a))
    
    x <- root_guess(f,a,K,x.guess)
    return(x)
  }
}

root_guess <- function(f, a, K, x.guess) {
  max.iter <- 32
  tol <- 1e-12
  x <- x.guess
  
  use.log <- (min(f)>=0)
  K.log <- log(K)
  
  #print('-----------------------------')
  
  for(k in 1:max.iter) {

    y.vec <- f*exp(a*x)
    
    if( use.log ) {
      y <- log(sum(y.vec)) - K.log
    } else {
      y <- sum(y.vec) - K
    } 
    
    #print( c( x, y ))
    if(abs(y)<tol) {
      break
    }
    
    if( use.log ) {
      dy <- sum(a*y.vec)/sum(y.vec)
    } else {
      dy <- sum(a*y.vec)
    }
    
    #print( c(x,y,dy) )
    
    x <- x -  y / dy;
  }
  
  if( k<max.iter || abs(y) < 10*tol ) {
    return(x)
  } else {
    stop( c('No root;', sprintf(' %f', c(f,a,K)), sprintf('%g', y)) ) 
  }
}

#' Multi-variate black-scholes price for a given the root
#' int( f*exp(a*z-0.5*a^2) - K ) z from root (= -d) to +inf
#'
#' @param f 
#' @param a
#' @param K
#' @param root
#' @param callput
#' 
#' @return option price
blks_N <- function(f, a, K, root, callput = TRUE){
  a.ext <- c(0.0, a)
  f.ext <- c(-K, f)

  # pnorm can handle +/- Inf correctly, so no special case
  p_put <- -sum( f.ext * pnorm(root - a.ext) )
  p_call <- sum( f.ext ) + p_put

  rv <- if( callput ) p_call else p_put
  return(rv)
}

#' Spread option without numerical root finding
#' 
#' @param strk strike price (array)
#' @param spot two spot prices
#' @param t.exp time to expiry
#' @param vol two volatilities
#' @param rho correlation between the two assets
#' @param r interest rate
#' @param d dividend rate
#' @param n.quad the quadrature size
#' @return spread option price (array)
blksmd_spreadquick <- function( strk, spot, t.exp, vol, rho, r = 0, d = 0, n.quad = 9 ){
  fwd <- spot * exp((r-d)*t.exp)
  
  var1 <- vol[1]*vol[1]*t.exp # variance for asset 1
  var2 <- vol[2]*vol[2]*t.exp # variance for asset 2
  facmat <- blksmd_2by2_factormat( A = var1, B = rho*vol[1]*vol[2]*t.exp, C = var2 )
  
  quad <- statmod::gauss.quad.prob(n.quad, dist='normal');
  xx <- quad$nodes
  ww <- quad$weights
  
  varV1 <- exp(0.5*facmat[1,1]*facmat[1,1])
  # solve AA * yy - BB / yy = K where yy = exp(facmat[1,1]*-dd)
  AA <- fwd[1] * exp(facmat[1,2]*xx - 0.5*var1)
  BB <- fwd[2] * exp(facmat[2,2]*xx - 0.5*var2)

  price.fwd <- rep(NA,length(strk))
  
  for( k in 1:length(strk) ) {
    yy <- ( strk[k] + sqrt( strk[k]*strk[k] + 4*AA*BB) )/( 2*AA )
    dd <- -log(yy) / facmat[1,1]
    call <- varV1*AA*pnorm(dd+facmat[1,1]) - varV1*BB*pnorm(dd-facmat[1,1]) - strk[k]*pnorm(dd)
    price.fwd[k] <- sum(call*ww)
    #print( dd )
    #print(call)
  }
  return( exp(-r*t.exp)*price.fwd )
}


#' An auxilirary function used in blksmd_spreadquick()
#' The factor matrix with the first column ( x, -x )
#' 
#' @param covar
#' @param A
#' @param B
#' @param C
#' @param direction
#' @return 2-by-2 factor matrix 
blksmd_2by2_factormat <- function( covar = NA, A = NA, B = NA, C = NA, direction = -1 ) {
  if( !is.na( covar )) {
    A <- covar[1,1]
    B <- covar[1,2]
    C <- covar[2,2]
  }
  
  det <- A*C-B*B
  if( det<=0 ) {
    stop( 'Matrix is not possitive definite, returning NA')
  }
  
  det.sqrt <- sqrt( det )
  
  #chol2 <- matrix( c(sqrt(A), B/sqrt(A), 0, sqrt(C-B*B/A)), ncol=2 )
  #a <- sqrt(A*C-B*B)/sqrt(A*A+A*C+2*A*B)
  #b <- (A+B)/sqrt(A*A+A*C+2*A*B)
  #Q <- matrix( c(a,-b,b,a), ncol=2 )
  
  if( direction < 0 ) {
    # the values in the first column have oppisite signs (for spread option)
    cholV <- matrix( c(det.sqrt, -det.sqrt, 
                       (A+B), (B+C) )/sqrt(A+2*B+C), ncol=2 )
  } else {
    # the values in the first column have same signs (for basket option)
    cholV <- matrix( c(det.sqrt, det.sqrt, 
                       (B-A), (C-B) )/sqrt(A-2*B+C), ncol=2 )
  }
  
  #validity check
  if( sum( ( err <- c(A,B,B,C)-as.vector(cholV %*% t(cholV)) ) * err ) > 1e-12 ) {
    stop('computed factor matrix can not reproduce the variance matrix. check the value')
  }
  
  return( cholV )
}
