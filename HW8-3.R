
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

# sigma function
calcSigma <- function(alpha, beta, phi){
    return(solve(diag(alpha) + beta * t(phi) %*% phi))
}

# mu function
calcMu <- function(beta, sigma, phi, tN){
    return(beta * sigma %*% t(Phi) %*% tN)
}


# gamma: used to udpate alpha and beta
calcGammaVec <- function(sigma, currentAlpha){
    diagSigma <- rowSums( sigma * diag(dim(sigma)[1]) )
    gamma <- 1 - currentAlpha * diagSigma
    return(gamma)
}

# update alpha 
updateAlpha <- function(mu, gamma){
    return(gamma / mu^2)
}

# update beta
updateBetaInv <- function(t, phi, mu, N, gamma ){
    dis.sqr <- dist(rbind(t, phi %*% mu))^2
    return(dis.sqr / (N - sum(gamma)))
}

