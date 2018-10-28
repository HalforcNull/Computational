library(mvtnorm)

#### HW chapter 7
####
#### Part 1
####
#### Fit a SVM to the following data
####
#### Part 1.1 using k(x,y) = exp(-0.5*t(x-y)%*%(x-y))
####
#### Part 1.2 using k(x,y) = <phi(x),phi(y)> where phi(x) = [dnorm(x,mu1,sig),dnorm(x,mu2,sig)] where mu1 = c(0,0), mu2=c(-1,-1), sig = diag(2)
####
#### Part 1.3 verify that the separation of 1.2 is linear in the space defined by phi()

set.seed(7675)
x1.tmp <- cbind(rmvnorm(50,c(0,0),matrix(c(0.05,0,0,0.05),2,2)),-1)
x2.tmp <- cbind(rmvnorm(50,c(1,1),matrix(c(0.05,0,0,0.05),2,2)),1)
x3.tmp <- cbind(rmvnorm(40,c(-1,-1),matrix(c(0.05,0,0,0.05),2,2)),1)

X <- rbind(x1.tmp,x2.tmp,x3.tmp)

# plot
plot(X[X[,3]==-1,1],X[X[,3]==-1,2],pch=16,col="green",xlim=c(-2,2),ylim=c(-2,2))
points(X[X[,3]==1,1],X[X[,3]==1,2],pch=16,col="lightblue")



