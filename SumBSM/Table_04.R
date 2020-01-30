# Table 4 (Parameter Set S1)
# Table 1 from Caldana & Fusai, the results identical to Hurd & Zhou with 15 nodes

source('blksmd.R')

t.exp <- 1
rho <- 0.5
corr <- diag(2)*(1-rho) + rho

vol <- c(20,10)/100
wts <- c(1,-1)
spot <- c(100,96)
strk <- seq(0,4,0.4)
r <- 0.1
d <- 0.05

# n.quad = 2, 3, or 4
# CV = T or F
system.time(pCall3 <- blksmd_basket( strk, spot, t.exp, vol, wts, corr, r=r, d=d, n.quad=4, CV=T ))
print(pCall3$call)
writeClipboard(as.character(pCall3$call))
