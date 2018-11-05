library(mvtnorm)
library(quadprog)
library(Matrix)

set.seed(7675)
x1.tmp <- cbind(rmvnorm(50,c(0,0),matrix(c(0.15,0,0,0.15),2,2)),1)
x2.tmp <- cbind(rmvnorm(50,c(1,1),matrix(c(0.15,0,0,0.15),2,2)),2)
x3.tmp <- cbind(rmvnorm(40,c(-1,-1),matrix(c(0.15,0,0,0.15),2,2)),2)

X <- rbind(x1.tmp,x2.tmp,x3.tmp)

# modify X to have correct class indices
X[X[,3]==1,3] <- -1
X[X[,3]==2,3] <- 1

plot(X[X[,3]==-1,1],X[X[,3]==-1,2],pch=16,col="green",xlim=c(-2,2),ylim=c(-2,2))
points(X[X[,3]==1,1],X[X[,3]==1,2],pch=16,col="lightblue")





