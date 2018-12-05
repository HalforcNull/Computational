#The goal is to design a RVM 2-classes classifier that can classify accurately all 44 objects


library(mvtnorm)
library(quadprog)
library(Matrix)


load('ink.training.rdata')

## Step 1. Data re-format
dat.training.S1 <- apply(ink.training.rdata[,,,1], 3, as.numeric)
dat.training.S2 <- apply(ink.training.rdata[,,,2], 3, as.numeric)

dat.training.X <- cbind(dat.training.S1, dat.training.S2)
dat.training.Tar <- c(rep(1, 22), rep(2,22))








