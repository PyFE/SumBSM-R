# Table 12 (Parameter Set A3)
# Previously tested in Linetsky (2004)

source('blksmd.R')

# Observation time and Simpson's rule weights for T=1 
n.step <- 200
t.obs <- (0:n.step)/n.step
wts <- c(1, 4, rep(c(2,4), n.step/2-1), 1)/(3*n.step)

# Observation time and Simpson's rule weights for T=2 
n.step <- 400
t.obs2 <- 2*(0:n.step)/n.step
wts2 <- c(1, 4, rep(c(2,4), n.step/2-1), 1)/(3*n.step)

# Case 1
spot <- 2
strk <- 2
r <- 0.02
vol <- 0.1

system.time(
  rv1 <- blksmd_asian( strk, spot, t.obs, vol, wts=wts, r=r, detail=F, CV=T )
)

# Case 2
spot <- 2
strk <- 2
r <- 0.18
vol <- 0.3

system.time(
  rv2 <- blksmd_asian( strk, spot, t.obs, vol, wts=wts, r=r, detail=F, CV=T )
)

# Case 3
# Notice t.obs2 and wts2 since T=2
spot <- 2
strk <- 2
r <- 0.0125
vol <- 0.25

system.time(
  rv3 <- blksmd_asian( strk, spot, t.obs2, vol, wts=wts2, r=r, detail=F, CV=T )
)

# Case 4, 5, and 6
spot <- 1
strk <- 2/c(1.9,2.0,2.1)
r <- 0.05
vol <- 0.5

rv456 <- blksmd_asian( strk, spot, t.obs, vol, wts=wts, r=r, detail=F, CV=T )
rv456$call <- rv456$call * c(1.9,2.0,2.1)
  
# Case 7
# Notice t.obs2 and wts2 since T=2
rv7 <- blksmd_asian( 2, 2, t.obs2, vol, wts=wts2, r=r, detail=F, CV=T )

# Exact value from Linetsky (2004)
rv_exact <- c(0.0559860415, 0.2183875466, 0.1722687410, 0.1931737903, 0.2464156905, 0.3062203648, 0.3500952190)

rv_all <- c(rv1$call, rv2$call, rv3$call, rv456$call, rv7$call)
print(rv_all)
writeClipboard(as.character(rv_all))