
library(mvtnorm)
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


# Gram matrix calculate:
#      since the x is one dim. don't need call apply anymore
# Calculate grad matrix
## sub function: calculate the diag
Gram_diag_element <- function(x, distanceFun){
    return(distanceFun(x,x))
}

GramMat <- function(xVec,kernelFun){
    Gram_diag <- diag( Gram_diag_element(xVec, kernelFun) ) # n * n diag
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
    return(c( rep(eps,length(tN)) - tN, rep(eps,length(tN)) + tN))
}

# AMat 
AMat <- function(mLen){
    IMat <- diag(mLen)
    ZeroMat <- matrix(0, nrow=mLen, ncol=mLen)
    mat1 <- c(rep(1, mLen), rep(-1, mLen))
    mat2 <- cbind(IMat, ZeroMat)
    mat3 <- cbind(ZeroMat, Imat)
    mat4 <- -mat2
    mat5 <- -mat3
    mat6 <- rep(1, 2*mLen)
    result <- rbind(mat1, mat2, mat3, mat4, mat5, mat6)
    rownames(result) <- NULL
    return(result)
}

# bvec
bVec <- function(HyperC, mLen, Nu){
    return(c( rep(0,mLen+1), rep(-HyperC/mLen, mLen), HyperC*Nu))
}


