# Parameter Set B1
# Krekel et al. (2004), Caldana et al. (2016)

source('blksmd.R')
options(digits=9)

# Table 6: Varying strike prices for the base parameters

n.var <- 4
t.exp <- 5
r <- 0
d <- 0
rho <- 0.5
corr <- diag(n.var)*(1-rho) + rho

vol <- rep(0.4,n.var)
wts <- rep(1/n.var, n.var)
spot <- rep(100, n.var)
strk = seq(50,150,10) #strk = 50

# Convergent Price (CP): n.quad=30
system.time(pCall1 <- blksmd_basket( strk, spot, t.exp, vol, wts, corr, n.quad=30 ))
writeClipboard(as.character(pCall1$call))
print(pCall1$call)

# Fast Price (FP): lambda=9
system.time(pCall2 <- blksmd_basket( strk, spot, t.exp, vol, wts, corr, lambda=9 ))
writeClipboard(as.character(pCall2$call))
print(pCall2$call - pCall1$call)


# Table 7: Varying correlations 

n.var <- 4
t.exp <- 5
vol <- rep(0.4, n.var)
rhos <- c(-10,10,30,50,80,95)/100
wts <- rep(1/n.var ,n.var)
spot <- rep(100, n.var)
pCallVal1 = rep(NA, length(rhos))
pCallVal2 = rep(NA, length(rhos))

for (k in 1:length(rhos)) {
  corr <- diag(n.var)*(1-rhos[k]) + rhos[k]
  
  # CP: n.quad=30
  pCall1 <- blksmd_basket( 100, spot, t.exp, vol, wts, corr, detail=T, n.quad=30 )
  print(pCall1$dim)  
  pCallVal1[k] <- pCall1$call[1]
  
  # FP: lambda=9
  pCall2 <- blksmd_basket( 100, spot, t.exp, vol, wts, corr, detail=T, lambda=9 )
  print(pCall2$dim)
  pCallVal2[k] <- pCall2$call[1]
}

writeClipboard(as.character(pCallVal1))
writeClipboard(as.character(pCallVal2))

print(pCallVal1)
print(pCallVal2-pCallVal1)

# Table 8: assymetric vol with fixed vol_1 = 100%

n.var <- 4
t.exp <- 5
vols <- c(5,10,20,40,60,80,100)/100
rho <- 0.5
corr <- diag(n.var)*(1-rho) + rho

wts <- rep(1/n.var, n.var)
spot <- rep(100, n.var)
pCallVal1 = rep(NA,length(vols))
pCallVal2 = rep(NA,length(vols))

for (k in 1:length(vols)) {
  volvec <- c(rep(vols[k],3),1)  # change first 3 vols only
  #volvec <- rep(vols[k],4)  # all vol simultaneously
  
  # CP: n.quad=30
  pCall1 <- blksmd_basket( 100, spot, t.exp, volvec, wts, corr, detail=T, n.quad=30 )
  pCallVal1[k] <- pCall1$call[1]
  
  # FP: lambda=9
  pCall2 <- blksmd_basket( 100, spot, t.exp, volvec, wts, corr, detail=T, lambda=9 )
  print( pCall2$dim )
  pCallVal2[k] <- pCall2$call[1]
}

writeClipboard(as.character(pCallVal1))
writeClipboard(as.character(pCallVal2))

print(pCallVal1)
print(pCallVal2-pCallVal1)