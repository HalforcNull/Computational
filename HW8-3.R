
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

