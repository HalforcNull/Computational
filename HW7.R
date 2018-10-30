library(mvtnorm)
library(quadprog)
library(proxy)
library(Matrix)
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

# kenrel 1
ker_1 <- function(x,y){
    return(exp(-0.5 * t(x-y) %*% (x-y)))
}

# Phi 2
Phi_2 <- function(x){
    dnorm_1 <- dmvnorm(x, c(-1,-1), diag(2))
    dnorm_2 <- dmvnorm(x, c(0,0), diag(2))
    return(c(dnorm_1, dnorm_2))
}

# kernel 2
ker_2 <- function(x,y){
    return(Phi_2(x) * Phi_2(y))
}

# Calculate grad matrix
## sub function: calculate the diag
## Only allow calculate one element 
Gram_diag_element <- function(x, distanceFun){
    return(distanceFun(x,x))
}

GramMat <- function(xVec,kernelFun){
    Gram_diag <- diag( apply(xVec, 1, Gram_diag_element, distanceFun = kernelFun) ) # n * n diag
    Gram_dist <- proxy::dist(xVec, method=kernelFun) # n-1 * n-1 mat
    
    return(as.matrix(Gram_dist) + Gram_diag)
}

# Dmat
DMat <- function(Tvec, GramMat){
    return(Tvec %*% t(Tvec) * GramMat)
    # In our case Tvec %*% t(Tvec) provides a 'sign' matrix for GramMat
}

# dvec
dVec <- function(n){
    return(rep(1,n))
}

# Amat
AMat <- function(tn){
    IMat <- diag(length(tn))
    return(cbind(tn, IMat))
}

# bvec
bVec <- function(n){
    return(rep(0,n+1))
}



#### Data:

xVec <- X[,1:2]
tN <- X[,3]

#### Part 1.1 using k(x,y) = exp(-0.5*t(x-y)%*%(x-y))
GramMatrix_1 <- GramMat(xVec, ker_1)
DMatrix_1 <- DMat(tN, GramMatrix_1)
DVector_1 <- dVec(length(tN))
AMatrix_1 <- AMat(tN)
bVector_1 <- bVec(length(tN))

Solution_1 <- solve.QP(nearPD(DMatrix_1)$mat, DVector_1, AMatrix_1, bVector_1, meq=1)

importantPointIndex <- which(Solution_1$solution > 0.001)

