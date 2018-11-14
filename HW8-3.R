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
    return(solve(nearPD(diag(alpha) + beta * t(phi) %*% phi)$mat))
}

# mu function
calcMu <- function(beta, sigma, Phi, tN){
    return(beta * sigma %*% t(Phi) %*% as.matrix(tN, nrow = length(tN)))
}


# gamma: used to udpate alpha and beta
calcGammaVec <- function(sigma, currentAlpha){
    diagSigma <- rowSums( sigma * diag(dim(sigma)[1]) )
    gamma <- 1 - currentAlpha * diagSigma
    return(gamma)
}

# update alpha 
calcAlpha <- function(mu, gamma){
    return(gamma / mu^2)
}

# update beta
calcBeta <- function(tN, Phi, mu, N, gamma ){
    dis.sqr <- dist(rbind(tN, t(Phi %*% mu)))^2
    return( (N - sum(gamma))/ dis.sqr )
}

calcY <- function(xNew, mu, xVec, kenrelFun){
    return( sum(mu * kernelFun(xVec, xNew)))
}


# main function :
xVec <- xN
GramMatrix <- GramMat(xVec, ker_3) 
my.N <- length(xN)
alpha = rep(0.1, dim(GramMatrix)[1])
beta = 0.1
sigma <- calcSigma(alpha, beta, GramMatrix)
mu <- calcMu(beta, sigma, GramMatrix, tN)

threshold.alpha <- 10000000


for(i in 1:1000){    
    tmp.gamma <- calcGammaVec(sigma, alpha)
    tmp.alpha <- calcAlpha(mu, tmp.gamma)
    #tmp.alpha <- ifelse(alpha > threshold.alpha , alpha, tmp.alpha)
    
 
    tmp.beta <- calcBeta(tN, GramMatrix, mu, my.N, tmp.gamma)
    tmp.beta <- as.numeric(tmp.beta)
    tmp.sigma <- calcSigma(tmp.alpha, tmp.beta, GramMatrix)
    tmp.sigma <- matrix(tmp.sigma, nrow = nrow(tmp.sigma), ncol=ncol(tmp.sigma))
    tmp.mu <- calcMu(tmp.beta, tmp.sigma, GramMatrix, tN)

    alpha <- tmp.alpha
    beta <- tmp.beta
    sigma <- tmp.sigma
    mu <- tmp.mu
    print(alpha)
}


pred.Y <- lapply( x, calcY, aN = model3.aN, 
                            aN_hat = model3.aN.hat, 
                            xVec = model3.xVec, 
                            b.num = b.3.num, 
                            kernelFun = ker_3 )




