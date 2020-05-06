# Table 9 (Parameter Set B2)
# Basket example: Milevsky & Posner (1998), Zhou and Wang (2008)

source('blksmd.R')

# Table 3
corr <- matrix( c(1, 0.35, 0.1, 0.27, 0.04, 0.17, 0.71,
                  0.35, 1, 0.39, 0.27, 0.5, -0.08, 0.15,
                  0.1, 0.39, 1, 0.53, 0.7, -0.23, 0.09,
                  0.27, 0.27, 0.53, 1, 0.46, -0.22, 0.32,
                  0.04, 0.5, 0.7, 0.46, 1, -0.29, 0.13,
                  0.17, -0.08, -0.23, -0.22, -0.29, 1, -0.03,
                  0.71, 0.15, 0.09, 0.32, 0.13, -0.03, 1), 7)

vol <- c(11.55, 20.68, 14.53, 17.99, 15.59, 14.62, 15.68)/100
wts <- c(10, 15, 15, 5, 20, 10, 25)/100
d <- c(1.69, 2.39, 1.36, 1.92, 0.81, 3.62, 1.66)/100
r <- 0.063

spot <- rep(100, 7)
strk <- c(80, 100, 120)

t.exp <- 3 # 0.5, 1, 2, 3

system.time(pCall3 <- blksmd_basket( strk, spot, t.exp, vol, wts, corr, 
                                     r = r, d = d, CV = T, lambda = 3, detail=T ))
writeClipboard(as.character(pCall3$call))
print(pCall3$call)