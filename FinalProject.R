#The goal is to design a RVM 2-classes classifier that can classify accurately all 44 objects


library(mvtnorm)
library(quadprog)
library(Matrix)


load('ink.training.rdata')

## Step 1. Data re-format
dat.training.S1 <- apply(ink.training.dat[,,,1], 3, as.numeric)
dat.training.S2 <- apply(ink.training.dat[,,,2], 3, as.numeric)

dat.training.X <- cbind(dat.training.S1, dat.training.S2)
dat.training.Tar <- c(rep(1, 22), rep(2,22))


dim(dat.training.X)
length(dat.training.Tar)


## Step 2. Support Functions

### 2.1 Kernel
### ker_1 is the one I used in previous HW
ker_1 <- function(x,y){
    return(exp(-0.5 * 16 * 1.414 * t(x-y) %*% (x-y)))
}

### ker_2 is the one I used in final project of Multi var
ker_2 <- function(x,y){
  return(cor(x,y))
}

### 2.2 Gram Matrix
Gram_diag_element <- function(x, distanceFun){
    return(distanceFun(x,x))
}

GramMat <- function(xVec,kernelFun){
    Gram_diag <- diag( apply(xVec, 1, Gram_diag_element, distanceFun = kernelFun) ) # n * n diag
    Gram_dist <- proxy::dist(xVec, method=kernelFun) # n-1 * n-1 mat
    
    return(as.matrix(Gram_dist) + Gram_diag)
}

### 2.3 sigma function
calcSigma <- function(alpha, beta, phi){
    return(solve(nearPD(diag(alpha) + beta * t(phi) %*% phi)$mat))
}

### 2.4 mu function
calcMu <- function(beta, sigma, Phi, tN){
    return(beta * sigma %*% t(Phi) %*% as.matrix(tN, nrow = length(tN)))
}

### 2.5 gamma: used to udpate alpha and beta
calcGammaVec <- function(sigma, currentAlpha){
    diagSigma <- rowSums( sigma * diag(dim(sigma)[1]) )
    gamma <- 1 - currentAlpha * diagSigma
    return(gamma)
}

### 2.6 update alpha 
calcAlpha <- function(mu, gamma){
    return(gamma / mu^2)
}

### 2.7 update beta
calcBeta <- function(tN, Phi, mu, N, gamma ){
    dis.sqr <- dist(rbind(tN, t(Phi %*% mu)))^2
    return( (N - sum(gamma))/ dis.sqr )
}

### 2.8 predict function
calcY <- function(xNew, mu, xVec, kenrelFun){
    return( sum(mu * kernelFun(xVec, xNew)))
}

## Step 3. Training Function and Predict Function
TrainingModel <- function(xN, tN){
    GramMatrix <- GramMat(xN, ker_1)
    my.N <- dim(xN)[1]+1
    phi <- cbind(rep(1, my.N-1), GramMatrix)
    alpha = rep(0.1, my.N)
    beta = 1
    sigma <- calcSigma(alpha, beta, phi)
    mu <- calcMu(beta, sigma, phi, tN) ## mu is not a value anymore

    threshold.alpha <- 1e7

    for(i in 1:1000){
        tmp.gamma <- calcGammaVec(sigma, alpha)
        tmp.alpha <- calcAlpha(mu, tmp.gamma)
        if(all(tmp.alpha > threshold.alpha)){
            break
        }
        tmp.alpha <- ifelse(alpha < threshold.alpha , tmp.alpha, alpha)
        
    
        tmp.beta <- calcBeta(tN, phi, mu, my.N, tmp.gamma)
        tmp.beta <- as.numeric(tmp.beta)
        tmp.sigma <- calcSigma(tmp.alpha, tmp.beta, phi)
        tmp.sigma <- matrix(tmp.sigma, nrow = nrow(tmp.sigma), ncol=ncol(tmp.sigma))
        tmp.mu <- calcMu(tmp.beta, tmp.sigma, phi, tN)

        alpha <- tmp.alpha
        beta <- tmp.beta
        sigma <- tmp.sigma
        mu <- tmp.mu

        print(max(alpha))
    }


    return(RVMModel)
}


DoPrediction <- function(model, xNew){


    return(results)
}

## Step 4. Function call

xN <- t(dat.training.X)
tN <- dat.training.Tar
