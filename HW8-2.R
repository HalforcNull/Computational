library(mvtnorm)
library(quadprog)
library(proxy)
library(Matrix)

set.seed(76757)
###### define hyperparamters
beta <- 1/0.2^2
x <- seq(0,1,length=1000)

###### generate a few targets and predictors
n.sim <- 20
xN <- runif(n.sim ,0,1)
tN <- sin(2*pi*xN) + rnorm(n.sim,0,sqrt(1/beta))

quartz()
plot(x,sin(2*pi*x),type="l",col="green")
points(xN,tN,col="blue")


# Changes from HW7
#   1. aN. In this case, we have a 2*n length vector, to contain c(aN, aN_hat)
#   2. Dmat. D_mat is changed since the quadratic part of the problem is changed
#   3. dvect. dvect changed as we discussed during class
#   4. Amat and bvect. Since the constrain of the problem is changed
#   5. Kernel function. ker_1 = abs(x1-x2)
#                       ker_2 = min( abs( (x1-0.25)^2 - (x2-0.25)^2 ), abs( (x1-0.75)^2 - (x2-0.75)^2) ))

# kernel 1:
ker_1 <- function(x,y){
    return(abs(x-y))
}

# kernel 2:
ker_2 <- function(x,y){
    return(min( abs( (x-0.25)^2 - (y-0.25)^2 ) , abs( (x-0.75)^2 - (y-0.75)^2 ) ))
}

# kenrel 3
ker_3 <- function(x,y){
    return(exp(-0.5 * 16 * 1.414 * (x-y) * (x-y)))
}


# Gram matrix calculate:
#      since the x is one dim. don't need call apply anymore
# Calculate grad matrix
## sub function: calculate the diag
Gram_diag_element <- function(x, distanceFun){
    return(distanceFun(x,x))
}

GramMat <- function(xVec,kernelFun){
    diag.vec <- Gram_diag_element(xVec, kernelFun)
    if(length(diag.vec) != length(xVec) ){
        Gram_diag <- diag( rep(diag.vec[1], length(xVec)) ) # n * n diag
    }else{
        Gram_diag <- diag( diag.vec ) # n * n diag
    }
    
    Gram_dist <- proxy::dist(xVec, method=kernelFun) # n-1 * n-1 mat
    
    return(as.matrix(Gram_dist) + Gram_diag)
}

# DMat Changed
DMat <- function(GramMat){
    mLen <- dim(GramMat)[1]
    PreMat <- rbind(diag(mLen), -diag(mLen))
    return(PreMat %*% GramMat %*% t(PreMat))
}

# DVec changed
dVec <- function(eps, tN){
    return(c(  tN - rep(eps,length(tN)), -rep(eps,length(tN)) - tN))
}

# AMat 
AMat <- function(mLen){
    IMat <- diag(mLen)
    ZeroMat <- matrix(0, nrow=mLen, ncol=mLen)
    mat1 <- c(rep(1, mLen), rep(-1, mLen))
    mat2 <- rep(-1, 2*mLen)
    mat3 <- cbind(IMat, ZeroMat)
    mat4 <- cbind(ZeroMat, IMat)
    mat5 <- -mat3
    mat6 <- -mat4
    
    result <- rbind(mat1, mat2, mat3, mat4, mat5, mat6)
    rownames(result) <- NULL
    return(t(result))
}

# bvec
bVec <- function(HyperC, mLen, Nu){
    return(c( 0, -HyperC*Nu, rep(0,2*mLen), rep(-HyperC/mLen, 2*mLen)))
}

# calculate B 7.69 and the equaltion we discuss during class
findB <- function(tN, aN, aN_Hat, Gmat, Nu ){
    bN <- tN - Nu - colSums((aN - aN_Hat) * Gmat)
    return(mean(bN))
}

# Predict
calcY <- function(xNew, aN, aN_hat, xVec, b.num, kernelFun){
    return( sum( (aN - aN_hat) * kernelFun(xVec, xNew) ) + b.num )
}

# data 

xVec <- xN
eps <- 0.3
HyperC <- 5
Nu <- 0.3

# kernel 1 = abs(x1-x2)
GramMatrix_1 <- GramMat(xVec, ker_1) 
DMatrix_1 <- DMat(GramMatrix_1)
DVector_1 <- dVec(eps, tN)
AMatrix_1 <- AMat(length(tN))
bVector_1 <- bVec(HyperC, length(tN), Nu)

Solution_1 <- solve.QP(nearPD(DMatrix_1)$mat, DVector_1, AMatrix_1, bVector_1, meq=1)
Solution_1$solution

importantPointIndex <- which(abs(Solution_1$solution) > 0.001)
#importantPointIndex <- ifelse(importantPointIndex > 20, importantPointIndex - 20, importantPointIndex)
importantPointIndex_aN <- importantPointIndex[which(importantPointIndex <= 20)]
importantPointIndex_aN_hat <- importantPointIndex[which(importantPointIndex > 20)] - 20
importantPointIndex <- ifelse(importantPointIndex > 20, importantPointIndex - 20, importantPointIndex)
importantPointIndex <- unique(importantPointIndex)
importantPointIndex

model.1.data <- rbind(xVec[importantPointIndex], tN[importantPointIndex])
model.1.data <- t(model.1.data)

plot(x,sin(2*pi*x),type="l",col="green")
points(xN,tN,col="blue",pch=16)
points(model.1.data[, 1], model.1.data[, 2], pch=1,col='red')


# data 
xVec <- xN
eps <- 0.1
HyperC <- 100
Nu <- 0.1

# kernel 2
GramMatrix_2 <- GramMat(xVec, ker_2) 
DMatrix_2 <- DMat(GramMatrix_2)
DVector_2 <- dVec(eps, tN)
AMatrix_2 <- AMat(length(tN))
bVector_2 <- bVec(HyperC, length(tN), Nu)

Solution_2 <- solve.QP(nearPD(DMatrix_2)$mat, DVector_2, AMatrix_2, bVector_2, meq=1)
Solution_2$solution

importantPointIndex <- which(abs(Solution_2$solution) > 0.001)
#importantPointIndex <- ifelse(importantPointIndex > 20, importantPointIndex - 20, importantPointIndex)
importantPointIndex_aN <- importantPointIndex[which(importantPointIndex <= 20)]
importantPointIndex_aN_hat <- importantPointIndex[which(importantPointIndex > 20)] - 20
importantPointIndex <- ifelse(importantPointIndex > 20, importantPointIndex - 20, importantPointIndex)
importantPointIndex <- unique(importantPointIndex)
importantPointIndex


model.2.data <- rbind(xVec[importantPointIndex], tN[importantPointIndex])
model.2.data <- t(model.2.data)

plot(x,sin(2*pi*x),type="l",col="green")
points(xN,tN,col="blue",pch=16)
points(model.2.data[, 1], model.2.data[, 2], pch=1,col='red')


model2.xVec <- xVec[importantPointIndex]
model2.tN <- tN[importantPointIndex]
model2.aN <- Solution_2$solution[importantPointIndex]
model2.aN.hat <- Solution_2$solution[importantPointIndex+20]
model2.Gmat <- GramMatrix_2[importantPointIndex,importantPointIndex]

b.2.num <- findB(model2.tN, model2.aN, model2.aN.hat, model2.Gmat, Nu)

pred.Y <- lapply( x, calcY, aN = model2.aN, 
                            aN_hat = model2.aN.hat, 
                            xVec = model2.xVec, 
                            b.num = b.2.num, 
                            kernelFun = ker_2 )

points(x, pred.Y, type='l',col='blue')
y.Max <- as.numeric(pred.Y) + eps
points(x, y.Max , type='l',col='blue')
y.Min <- as.numeric(pred.Y) - eps
points(x,y.Min, type='l',col='blue')




# kernel 3
xVec <- xN
eps <- 0.3
HyperC <- 100
Nu <- 0.1


GramMatrix_3 <- GramMat(xVec, ker_3) 
DMatrix_3 <- DMat(GramMatrix_3)
DVector_3 <- dVec(eps, tN)
AMatrix_3 <- AMat(length(tN))
bVector_3 <- bVec(HyperC, length(tN), Nu)

Solution_3 <- solve.QP(nearPD(DMatrix_3)$mat, DVector_3, AMatrix_3, bVector_3, meq=1)
Solution_3$solution

importantPointIndex <- which(abs(Solution_3$solution) > 0.001)
#importantPointIndex <- ifelse(importantPointIndex > 20, importantPointIndex - 20, importantPointIndex)
importantPointIndex_aN <- importantPointIndex[which(importantPointIndex <= 20)]
importantPointIndex_aN_hat <- importantPointIndex[which(importantPointIndex > 20)]
importantPointIndex <- ifelse(importantPointIndex > 20, importantPointIndex - 20, importantPointIndex)
importantPointIndex <- unique(importantPointIndex)
importantPointIndex

model.3.data <- rbind(xVec[importantPointIndex], tN[importantPointIndex])
model.3.data <- t(model.3.data)

plot(x,sin(2*pi*x),type="l",col="green", ylim=(c(-1.5, 1.5)))
points(xN,tN,col="blue",pch=16)
text(xN, tN, labels = seq(1,20,1))
points(model.3.data[, 1], model.3.data[, 2], pch=1,col='red')

model3.xVec <- xVec[importantPointIndex]
model3.tN <- tN[importantPointIndex]
model3.aN <- Solution_3$solution[importantPointIndex]
model3.aN.hat <- Solution_3$solution[importantPointIndex+20]
model3.Gmat <- GramMatrix_3[importantPointIndex,importantPointIndex]

b.3.num <- findB(model3.tN, model3.aN, model3.aN.hat, model3.Gmat, Nu)

pred.Y <- lapply( x, calcY, aN = model3.aN, 
                            aN_hat = model3.aN.hat, 
                            xVec = model3.xVec, 
                            b.num = b.3.num, 
                            kernelFun = ker_3 )

points(x, pred.Y, type='l',col='blue')
y.Max <- as.numeric(pred.Y) + eps
points(x, y.Max , type='l',col='blue')
y.Min <- as.numeric(pred.Y) - eps
points(x,y.Min, type='l',col='blue')



pred.Y2 <- NULL

for(i in 1:length(x)){
    tmp.y <- calcY(x[i], model3.aN, model3.aN.hat,  model3.xVec,b.3.num,ker_3 )
    pred.Y2 <- c(pred.Y2, tmp.y)
}


