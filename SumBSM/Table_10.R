# Table 10 (Parameter Set A1)
# Carverhill&Clewlow(1990) Benhamou FFT(Table 1), Cerny & Kyriakou (Table 3)

source('blksmd.R')

spot <- 100
strk <- c(80, 90, 100, 110, 120)
n.obs <- 50
t.obs <- seq(0,1,1/n.obs)
r <- 0.1
vol <- 0.1  # 0.3 or 0.5

system.time( rv <- blksmd_asian( strk, spot, t.obs, vol, r=r, detail=F, CV=T ) )
print(rv)
writeClipboard(as.character(rv$call))
