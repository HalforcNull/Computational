### Midterm 721 - FA2018
###
### Take home part
###
### 1) Use the training data below to build a Gaussian Process classifier Section 6.4.5 and 6.4.6 in the text book
###    You do not need to optimise the parameters of the kernel function using equations 6.89 and following
###    The kernel function that you need to use is k(a,b) = exp(-0.5*t(a-b)%*%(a-b)), where a,b are multidimensional vectors
###    The main deliverable is a function that accepts two data frames as arguments.
###    The first data frame contains training data in 3 columns (x,y,class) and the
###    second data frame contains any number of new samples in 2 columns (x,y)
###    The function has to return the probability OF CLASS 1 for all the points in the second data frame
###
###    You can use the optim() function in R. The rest of the algorithm needs to be coded from scratch apart from trivial R functions 
###
### 2) Provide the COMPLETE derivation for exercise 6.22 in the text book
###
### DEADLINE: Friday October 19 @ 0900
### There will be no extension unless medically related
###
### WORK BY YOURSELF!!!!
###
### Code will be audited
###
### Students with duplicated code or analytical solution for the derivation will receive a F

### make sure you set the seed correctly
set.seed(7677)

### generate data
n <- 50
w <- rgamma(n,shape=6,rate=1)
x <- w + rgamma(n,10,1)
y <- rgamma(n,shape=3,rate=0.3)
plot(density(y))
z <- 20-y + rgamma (n,8,2)

### look at the data
quartz()
par(mfrow=c(2,2))
plot(w,x,col="red",pch=16,xlim=c(0,15),ylim=c(0,30))
points(y,z,col="blue",pch=16)
plot(density(x),col="red",lwd=2)
lines(density(z),col="blue", lwd=2)
plot(density(w),col="red",lwd=2)
lines(density(y),col="blue", lwd=2)

### put the data together
data.df <- data.frame(x=c(w,y), y=c(x,z), class = rep(c(1,2),each=n))

### use data.df

### Below is my solution
### Start @ 10/18/2018 1:00 pm
### Code won't be pushed to https://github.com/HalforcNull/Computational untill 10/19/2018 10:00 am

library(Matrix)

# Kernel function
## As told above
fun.kernel <- function(x,b){
    return(as.numeric(exp(-0.5*as.matrix(x-b)%*%as.matrix(t(x-b)))))
}

# Gram Matrix
## Instead of calculate gram matrix directly, I need a support function
## This support function takes training data and 2 index, then find the kernel of indexed data
fun.gram.support <- function(training, x){
    a <- training[x[1],]
    b <- training[x[2],]
    return(ifelse(x[1]==x[2], 1,fun.kernel(a,b)))
}
## Recall our K matrix is divided by alpha
fun.gram <- function(training, alpha){
    indexs <- seq(1, nrow(training), 1)
    combIndex <- expand.grid(indexs, indexs)
    my.result <- apply(combIndex, 1, fun.gram.support, training=training)
    my.result <- matrix(my.result, nrow=nrow(training))
    return(my.result/alpha)
}

# Calculate CN
fun.CN <- function(Gram, beta){
    gram.dim <- nrow(Gram)
    return(Gram + diag(1/beta, nrow=gram.dim))
}

# Psi(aN) and related functions
## functions need for Psi(aN)
fun.sigmoid <- function(aN){
    return(1/(1+exp(-aN)))
}
## by 6.80 
## remove the parts does not contain aN ( we are optimize Psi wrt. aN, so any part without aN can be treated as constant )
fun.Psi.aN <- function(aN, tN, CN){
    my.N <- length(tN)
    return(as.numeric(-1/2*t(aN) %*% solve(nearPD(CN)$mat) %*% aN + t(tN) %*% aN))
}

fun.Psi.gridiant <- function(aN, tN, CN){
    sigma.aN <- fun.sigmoid(aN)
    return(tN - sigma.aN - solve(nearPD(CN)$mat) %*% aN)
}

## WN function
fun.WN <- function(aN.mode){
    my.N <- length(aN.mode)
    sigma.aN <- fun.sigmoid(aN.mode)
    WN <- sigma.aN * (1 - sigma.aN)
    return(diag(WN))
}

## Hassian
fun.Hassian <- function(WN, CN ){
    return(WN + solve(nearPD(CN)$mat))
}

# mu_a_new and var_a_new 
## Function to find k vector (distance of new observation wrt. training)
## Note: this function also need devided by alpha
fun.k.vec <- function(dat.new, training, alpha){
    return( apply(training, 1, fun.kernel, b=dat.new) / alpha)
}
## Function to find c 
fun.c <- function(dat.new, alpha, beta){
    return(as.numeric(fun.kernel(dat.new, dat.new)/alpha + 1/beta))
}

## by 6.87
fun.mean.new <- function(k.vec, tN, aN.mode){
    sigma.aN <- fun.sigmoid(aN.mode)
    return(as.numeric(k.vec%*%as.matrix(tN-sigma.aN)))
}

## byt 6.88
fun.var.new <- function(c.value, k.vec, WN, CN ){
    return(c.value - as.numeric(t(k.vec) %*% solve(solve(nearPD(WN)$mat) + CN) %*% k.vec))
}

# kappa solution
fun.kappa.solution <- function(mu.a, var.a){
    kappa <- 1/sqrt(1 + pi*var.a/8)
    return(1/(1+exp(-kappa %*% t(mu.a))))
}

# predict function
fun.getModel <- function(training, tN, alpha, beta){
    ## find CN
    m.Gram <- fun.gram(training, alpha)
    m.CN <- fun.CN(m.Gram, beta)
    ## find aN.mode
    aN.mode <- optim(rep(0,100), fn=fun.Psi.aN, gr=fun.Psi.gridiant, tN=tN, CN=m.CN)
    ## WN
    m.WN <- fun.WN(aN.mode$par)
    return(list(training=training, alpha=alpha, beta=beta, tN = tN, aN.mode = aN.mode$par,WN = m.WN, CN=m.CN))
}

fun.predict <- function(dat.new, my.model){
    k.vec <- fun.k.vec(dat.new, my.model$training, my.model$alpha)
    c.value <- fun.c(dat.new, my.model$alpha, my.model$beta)
    mu.a <- fun.mean.new(k.vec, my.model$tN, my.model$aN.mode)
    var.a <- fun.var.new(c.value, k.vec, my.model$WN, my.model$CN)
    return(fun.kappa.solution(mu.a, var.a))
}

### above is the function we used for classification






########################################################################################################
###       Delivered Function: 
###           fun.predict.prob.in.class1
###       Input: 
###           trainingData: training data,     N * 3 matrix, first 2 column is data, 3rd column is class
###           newData: testing data,      n * 2 matrix, first 2 column is data
###       Output:
###           A vector contains the prob of being in class 1
########################################################################################################

fun.predict.prob.in.class1 <- function(trainingData, newData){
    if(ncol(trainingData) != 3){
      print('We need have a N * 3 matrix for training')
      return(NULL)
    }
  
    if(ncol(newData) != 2){
      print('We need have a n * 2 matrix for testing')
      return(NULL)
    }
    
    # Step 1. training Data
    training.dat <- trainingData[, -3]
    training.t <- trainingData[, 3] - 1
    
    # Step 2. Hyper parms
    alpha = 1
    beta = 1
    
    # Step 3. get model
    model.trained <- fun.getModel(training.dat, training.t, alpha, beta)
    
    # Step 4. do prediction on training data
    # pre.label <- apply(training.dat, 1, fun.predict, my.model=model.trained)
    # for some reason, this apply function is not work well
    pre.prob <- NULL
    for( i in 1:nrow(newData) ){
      pre.prob <- c(pre.prob, fun.predict(newData[i,], model.trained))
    }
    
    return(1-pre.prob)
}

###                 END OF FUNCTION       




### try classify training data
pred.prob <- fun.predict.prob.in.class1(data.df, data.df[,-3])
ifelse(pred.prob > 0.5, 1, 2)
data.df[,3]

## Finish prediction code @ 2:49. Start debug now
## Finish Debug @ 4:25
## Finish delivered function @ 4:49

