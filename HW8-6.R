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


# main function :
xVec <- X[,1:2]
tN <- X[,3]

GramMatrix <- GramMat(xVec, ker_1) 
my.N <- length(xN)
alpha = rep(0.00001, dim(GramMatrix)[1])
beta = 0.00001
sigma <- calcSigma(alpha, beta, GramMatrix)
mu <- calcMu(beta, sigma, GramMatrix, tN)

for(i in 1:1000){
    tmp.gamma <- calcGammaVec(sigma, alpha)
    tmp.alpha <- calcAlpha(mu, tmp.gamma)
    tmp.beta <- calcBeta(tN, GramMatrix, mu, my.N, tmp.gamma)
    tmp.beta <- as.numeric(tmp.beta)
    tmp.sigma <- calcSigma(tmp.alpha, tmp.beta, GramMatrix)
    tmp.sigma <- matrix(tmp.sigma, nrow = nrow(tmp.sigma), ncol=ncol(tmp.sigma))
    tmp.mu <- calcMu(tmp.beta, tmp.sigma, GramMatrix, tN)
    if( sum( abs(tmp.mu - mu) ) < 0.0001){
        writeLines('found result')
        break
    }
    
    alpha <- as.vector(tmp.alpha)
    beta <- tmp.beta
    sigma <- tmp.sigma
    mu <- tmp.mu
    print(alpha)
}


