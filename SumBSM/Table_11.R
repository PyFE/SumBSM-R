# Table 11 (Parameter Set A2)
# Fusai (Levy 2008), Cerny & Kyriakou (Table 4), Cai Li Shi (Table 3&5)

source('blksmd.R')

spot <- 100
strk <- c(90, 100, 110)
r <- 0.0367
vol <- 0.17801
n.obs <- 12  # or 50 or 250
t.obs <- seq(0,1,1/n.obs)

system.time( rv <- blksmd_asian( strk, spot, t.obs, vol, r=r, detail=F, CV=T ) )
print(rv)
writeClipboard(as.character(rv$call))
