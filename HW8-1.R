library(mvtnorm)
library(quadprog)
library(Matrix)

set.seed(7675)
x1.tmp <- cbind(rmvnorm(50,c(0,0),matrix(c(0.15,0,0,0.15),2,2)),1)
x2.tmp <- cbind(rmvnorm(50,c(1,1),matrix(c(0.15,0,0,0.15),2,2)),2)
x3.tmp <- cbind(rmvnorm(40,c(-1,-1),matrix(c(0.15,0,0,0.15),2,2)),2)

X <- rbind(x1.tmp,x2.tmp,x3.tmp)



# Solution: Most of the HW7 code are reusable. 
# For the reusable part, directly copy from HW7.R
# Changes:
#   1. dvec, since the minimize target function changed
#   2. AMat and bvect, since the constrain is changed

# modify X to have correct class indices
X[X[,3]==1,3] <- -1
X[X[,3]==2,3] <- 1

plot(X[X[,3]==-1,1],X[X[,3]==-1,2],pch=16,col="green",xlim=c(-2,2),ylim=c(-2,2))
points(X[X[,3]==1,1],X[X[,3]==1,2],pch=16,col="lightblue")


# kenrel 1
ker_1 <- function(x,y){
    return(exp(-0.5 * t(x-y) %*% (x-y)))
}

# Calculate grad matrix
## sub function: calculate the diag
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
# This is different from HW7
dVec <- function(n){
    return(rep(0,n))
}

# Amat 
# This is different from HW7
AMat <- function(tn){
    IMat <- diag(length(tn))
    OneVect <- rep(1, length(tn))
    return(cbind(cbind(tn, IMat), cbind(-IMat, OneVect)))
}

# bvec
# This is different from HW7
bVec <- function(n, Nu){
    c1 <- rep(0,n+1)
    c2 <- rep(-1/n, n)
    c3 <- Nu
    return(c(c1,c2,c3))
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

xVec <- X[,1:2]
tN <- X[,3]

#### Data: Predicting
one.dim <- seq(-2,2,0.01)
all.x <- expand.grid(x = one.dim, y = one.dim)


Nu <- 0.1 ## I used several different Nu, and decide pick 0.1

GramMatrix_1 <- GramMat(xVec, ker_1)
DMatrix_1 <- DMat(tN, GramMatrix_1)
DVector_1 <- dVec(length(tN))
AMatrix_1 <- AMat(tN)
bVector_1 <- bVec(length(tN), Nu)
Solution_1 <- solve.QP(nearPD(DMatrix_1)$mat, DVector_1, AMatrix_1, bVector_1, meq=1)


## Model and Model used points
importantPointIndex <- which(abs(Solution_1$solution) > 0.001)
model.1.data <- xVec[importantPointIndex,]
b_1 <- findB(tN[importantPointIndex], Solution_1$solution[importantPointIndex], GramMatrix_1[importantPointIndex,importantPointIndex])
plot(X[X[,3]==-1,1],X[X[,3]==-1,2],pch=16,col="green",xlim=c(-2,2),ylim=c(-2,2))
points(X[X[,3]==1,1],X[X[,3]==1,2],pch=16,col="lightblue")
points(model.1.data[, 1], model.1.data[, 2], pch=1,col='red')


## Draw Prediction and Model used points
## This may consume more than 1 minute
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
