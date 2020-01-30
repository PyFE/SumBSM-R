# Table 5 (Parameter Set S2)
# Our own SKEW example, not used in literature

source('blksmd.R')

t.exp <- 1
wts <- c(1,-1)

vol <- c(15,30)/100
spot <- c(200,100)
strk <- 100
r <- 0
d <- 0

rhos <- seq(0.9,-0.9,-0.2)

pCallVal2 = rep(NA,length(rhos))
for (k in 1:length(rhos)) {
  corr <- diag(2)*(1-rhos[k]) + rhos[k]
  # CP: lambda=9, FP: lambda=3
  pCall2 <- blksmd_basket( strk, spot, t.exp, vol, wts, corr, r=r, d=d, detail=T, lambda=9)
  print( pCall2$dim )
  pCallVal2[k] <- pCall2$call[1]
}
writeClipboard(as.character(pCallVal2))


#
# Other methods for benchmarking
#
pCallVal1 = rep(NA,length(rhos))
for (k in 1:length(rhos)) {
  # Uncomment only one line below
  
  # BjSt method
  pCallVal1[k] <- blksmd_bjerkspread(strk, spot, t.exp, vol, rhos[k])
  
  # Lo method
  #pCallVal1[k] <- blksmd_splittingpread(strk, spot, t.exp, vol, rhos[k], Kirk=F)
  
  # LDZ method
  #pCallVal1[k] <- blksmd_CalcSpreadOptLDZ(strk, spot, t.exp, vol, rhos[k])
  
  # This is a quick method, NOT used in the paper.
  #pCallVal1[k] <- blksmd_spreadquick(strk, spot, t.exp, vol, rhos[k])
}
writeClipboard(as.character(pCallVal1))
