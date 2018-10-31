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
    return(Phi_2(x) %*% Phi_2(y))
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

# calculate B
findB <- function(tN, aN, Gmat){
    bN <- tN - colSums(aN * tN * Gmat)
    return(mean(bN))
}

# prediction
calcY <- function(xNew, tN, aN, xVec, b.num, kenrelFun){
    GramCol <- apply(xVec, 1, kenrelFun, y=xNew)
    return( sum(aN*tN*GramCol) + b.num )
}

#### Data: Training

xVec <- X[,1:2]
tN <- X[,3]

#### Data: Predicting
one.dim <- seq(-2,2,0.01)
all.x <- expand.grid(x = one.dim, y = one.dim)


#### Part 1.1 using k(x,y) = exp(-0.5*t(x-y)%*%(x-y))
GramMatrix_1 <- GramMat(xVec, ker_1)
DMatrix_1 <- DMat(tN, GramMatrix_1)
DVector_1 <- dVec(length(tN))
AMatrix_1 <- AMat(tN)
bVector_1 <- bVec(length(tN))

Solution_1 <- solve.QP(nearPD(DMatrix_1)$mat, DVector_1, AMatrix_1, bVector_1, meq=1)

importantPointIndex <- which(Solution_1$solution > 0.1)
model.1.data <- xVec[importantPointIndex,]

b_1 <- findB(tN[importantPointIndex], Solution_1$solution[importantPointIndex], GramMatrix_1[importantPointIndex,importantPointIndex])
prediction_1 <- apply(  all.x, 1, calcY, tN = tN[importantPointIndex], 
                        aN = Solution_1$solution[importantPointIndex], 
                        xVec = xVec[importantPointIndex,], 
                        b.num = b_1,
                        kenrelFun = ker_1)

predic_1_matrix <- matrix(prediction_1, nrow = length(one.dim))
contour(x=one.dim, y=one.dim, z=predic_1_matrix, level=c(-1,0,1), xlim = c(-2,2), ylim=c(-2,2))
title("Part 1.1 Training Data and Model used data (red circled)")
points(X[X[,3]==-1,1],X[X[,3]==-1,2],pch=16,col="green",xlim=c(-2,2),ylim=c(-2,2))
points(X[X[,3]==1,1],X[X[,3]==1,2],pch=16,col="lightblue")
points(model.1.data[, 1], model.1.data[, 2], pch=1,col='red')


#### Part 1.2 using k(x,y) = <phi(x),phi(y)> where phi(x) = [dnorm(x,mu1,sig),dnorm(x,mu2,sig)] where mu1 = c(0,0), mu2=c(-1,-1), sig = diag(2)
GramMatrix_2 <- GramMat(xVec, ker_2)
DMatrix_2 <- DMat(tN, GramMatrix_2)
DVector_2 <- dVec(length(tN))
AMatrix_2 <- AMat(tN)
bVector_2 <- bVec(length(tN))

Solution_2 <- solve.QP(nearPD(DMatrix_2)$mat, DVector_2, AMatrix_2, bVector_2, meq=1)
importantPointIndex <- which(Solution_2$solution > 0.1)
model.2.data <- xVec[importantPointIndex,]

b_2 <- findB(tN[importantPointIndex], Solution_2$solution[importantPointIndex], GramMatrix_2[importantPointIndex,importantPointIndex])
prediction_2 <- apply(  all.x, 1, calcY, tN = tN[importantPointIndex], 
                        aN = Solution_2$solution[importantPointIndex], 
                        xVec = xVec[importantPointIndex,], 
                        b.num = b_2,
                        kenrelFun = ker_2)

predic_2_matrix <- matrix(prediction_2, nrow = length(one.dim))
contour(x=one.dim, y=one.dim, z=predic_2_matrix, level=c(-1,0,1), xlim = c(-2,2), ylim=c(-2,2))
points(X[X[,3]==-1,1],X[X[,3]==-1,2],pch=16,col="green",xlim=c(-2,2),ylim=c(-2,2))
title("Part 1.2 Training Data and Model used data (red circled)")
points(X[X[,3]==1,1],X[X[,3]==1,2],pch=16,col="lightblue")
points(model.2.data[, 1], model.2.data[, 2], pch=1,col='red')


#### Part 1.3 verify that the separation of 1.2 is linear in the space defined by phi()

PhiSpaceLoc <- apply(xVec,1,Phi_2)
PhiSpaceLoc <- t(PhiSpaceLoc)

model.pts <- apply(model.2.data,1,Phi_2)
model.pts <- t(model.pts)

phi.all.x <- apply(all.x,1,Phi_2)
phi.all.x <- t(phi.all.x)
# all.y <- phi.all.x %*% weight + linear.model.b

greenLineIdx <- which(abs(prediction_2 + 1) < 0.01)
blackLineIdx <- which(abs(prediction_2) < 0.01)
blueLineIdx <- which(abs(prediction_2 - 1) < 0.01)

plot(PhiSpaceLoc[greenIdx,1], PhiSpaceLoc[greenIdx,2], pch=16,col="green",xlim=c(0,0.2),ylim=c(0,0.2))
title("Part 1.2 Training Data and Decision on Phi Space")
points(PhiSpaceLoc[blueIdx,1], PhiSpaceLoc[blueIdx,2], pch=16,col="lightblue")
points(model.pts[, 1], model.pts[, 2], pch=1,col='red')

points(phi.all.x[greenLineIdx,1], phi.all.x[greenLineIdx,2], type='l', col="green")
points(phi.all.x[blackLineIdx,1], phi.all.x[blackLineIdx,2], type='l', col="black")
points(phi.all.x[blueLineIdx,1], phi.all.x[blueLineIdx,2], type='l', col="blue")



