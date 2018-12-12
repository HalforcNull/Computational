#The goal is to design a RVM 2-classes classifier that can classify accurately all 44 objects


library(mvtnorm)
library(quadprog)
library(Matrix)


load('ink.training.rdata')

## Step 1. Data re-format
dat.training.S1 <- apply(ink.training.dat[,,,1], 3, as.numeric)
dat.training.S2 <- apply(ink.training.dat[,,,2], 3, as.numeric)

dat.training.X <- cbind(dat.training.S1, dat.training.S2)
dat.training.Tar <- c(rep(0, 22), rep(1,22)) #  To start with, we consider two-class problems
                                             #  with a binary target variable t ∈ {0, 1}.


dim(dat.training.X)
length(dat.training.Tar)

## Step 2. Support Functions

### 2.1 Kernel
### ker_1 is the one I used in previous HW
ker_1 <- function(x,y){
    return(exp(-0.5 * t(x-y) %*% (x-y)))
}

### ker_2 is the one I used in final project of Multi var
ker_2 <- function(x,y){
  return(cor(x,y))
}

ker_3 <- function(x,y){
    return(1/(cor(x,y))
}

ker_4 <- function(x,y){
    return(exp(10-cor(x,y)) - exp(9))
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

### 2.3 Sigmoid Function
Sigmoid.fun <- function(w, Phi_X){
  a = t(w) %*% Phi_X
  return(as.vector(1/(1+exp(-a))))
}

### 2.4 optimization target (used to be 4.142, now is 7.109)
#### 2.4.1 main ln p function
ln.w.t <- function(w.vec, A.Mat, t.vec, PhiMat){
    y.vec <- Sigmoid.fun(w.vec, PhiMat)
    part.1 = sum( t.vec * log2(y.vec) )
    part.2 = sum( (1-t.vec) * log2(1-y.vec))
    part.3 = -0.5 * t(w.vec) %*% A.Mat %*% w.vec

    return(part.1+part.2+part.3)
}

opmtTarget <- function(w.vec, A.Mat, t.vec, PhiMat){
  return(-ln.w.t(w.vec, A.Mat, t.vec, PhiMat))
}


#### 2.4.2 A.mat function
calAMat <- function(alpha){
    return(diag(alpha))
}


### 2.5 Gradient function
gradientFunc <- function(w.vec, t.vec, A.Mat, PhiMat){
    y.vec <- Sigmoid.fun(w.vec, PhiMat)
    return(t(PhiMat) %*% (t.vec-y.vec) - A.Mat%*%w.vec)
}

### 2.7 gamma: used to udpate alpha and beta
calcGammaVec <- function(sigma, currentAlpha){
    diagSigma <- rowSums( sigma * diag(dim(sigma)[1]) )
    gamma <- 1 - currentAlpha * diagSigma
    return(gamma)
}

### 2.8 update alpha 
calcAlpha <- function(w.new, gamma){
    return(gamma / w.new^2)
}

### 2.9 predict function
calcY <- function(xNew, xModel, wModel, HassianModel, kernelFun){
    K.new <- apply(xModel, 1, kernelFun, y = xNew)
    mu.a <- t(wModel) %*% K.new
    var.a <- t(K.new) %*% -HassianModel %*% K.new
    kappa <- 1/sqrt(1 + pi*var.a/8)
    return(1/(1+exp(-kappa %*% t(mu.a))))
}

## Step 3. Training Function and Predict Function
TrainingModel <- function(xN, tN){
    # calc gram matrix and size of sample
    GramMatrix <- GramMat(xN, ker_4) / 100
    my.N <- dim(xN)[1]

    # init 
    my.alpha <- rep(1, my.N)
    threshold.alpha <- 1e7
    converge.threshold <- 1e-6
    
    last.alpha <- rep(0, my.N)  ## used for checking converge
    my.weight <- rep(0, my.N)

    for(i in 1:5000){
        # calc/init optim need parms
        tmp.A.Mat <- calAMat(my.alpha)
        
        optimResult <- optim(par=my.weight, fn=opmtTarget, hessian=TRUE, gr=gradientFunc, A.Mat = tmp.A.Mat, t.vec = tN, PhiMat = GramMatrix)

        weight.new <- optimResult$par
        sigma.new <- solve( -optimResult$hessian )
        gamma.new <- calcGammaVec(sigma.new, my.alpha)
        alpha.new <- calcAlpha(weight.new, gamma.new)

        if(sum(abs(alpha.new-last.alpha)/last.alpha) < converge.threshold){
            break
        }else{
            last.alpha <- alpha.new
        }
        
        my.alpha <- ifelse(abs(alpha.new) > threshold.alpha, my.alpha, alpha.new)
        #my.alpha <- alpha.new
        #my.alpha
        #my.weight <- weight.new # don't need update weight
    }

    selectedNodes <- which(weight.new > 0.001 )

    RVMModel=list()
    RVMModel$selectX <- xN[selectedNodes,]
    RVMModel$selectW <- my.weight[selectedNodes]
    RVMModel$Hessian <- optimResult$hessian[selectedNodes,selectedNodes]
    RVMModel$selectT <- tN[selectedNodes]
    return(RVMModel)
}


DoPrediction <- function(model, xNew){

    Y.vec <- calcY(xNew, model$selectX, model$selectW, model$Hessian, ker_4)

    return(ifelse( Y.vec  < 0.5, 1,2))
}

## Step 4. Function call


xN <- t(dat.training.X)
tN <- dat.training.Tar
## 

my.model <- TrainingModel(xN, tN)
dim(my.model$selectX)

i = 1
result <- NULL
for( i in 1:44){
    result <- c( result, as.numeric(DoPrediction(my.model, xN[i,])))
}

result
##xNew <- xN[1,]
##xModel <- my.model$selectX
##wModel <- my.model$selectW
##HassianModel <- my.model$Hessian
##kernelFun <- ker_4



