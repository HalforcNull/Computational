### data for Least square multiclass discriminant function and Fisher multiclass discriminant function
library(mvtnorm)

# define number of objects in each class
N <- 20

# define common covariance matrix 
Sigma <- matrix(c(10,5,5,4),2,2)

# class 1
c1.dat <- rmvnorm(N,mean=c(0,0),Sigma )

# class 2
c2.dat <- rmvnorm(N,mean=c(10,10),Sigma )

# class 3
c3.dat <- rmvnorm(N,mean=c(-5,5),Sigma )

# create a data.frame of the data 
dat <- data.frame(X = c(c1.dat[,1],c2.dat[,1],c3.dat[,1]), Y = c(c1.dat[,2],c2.dat[,2],c3.dat[,2]), Class = as.factor(rep(c(1,2,3),each=N)))

# plot the data 
plot(Y~X,data=dat, col=Class)

### use dat to solve the exercise 



head(dat)

# Step 1 Solve by least square:
# Step 1.1 find weight matrix
LSD.W <- function(dat){
  X <- as.matrix(
    cbind( rep(1, nrow(dat)), dat[,1:2] )
  )
  
  # if we want more classes, we need modify following code
  TarIndex <- as.numeric(dat[, 3])
  Tar <- matrix(0, nrow = nrow(dat), ncol = 3)
  for( i in 1:nrow(dat)){
    Tar[i,TarIndex[i]] = 1
  }
  # if we want more classes, we need modify code above
  
  return(solve(t(X)%*%X)%*%t(X)%*%Tar)
}

lsdw <- LSD.W(dat)
lsdw


# confirm weight we got
# 
# X <- as.matrix(
#   cbind( rep(1, nrow(dat)), dat[,1:2] )
# )
# 
# Y <-  X%*%lsdw 
# 
# for(i in 1:nrow(Y)){
#   print(which.max(Y[i,]))
# }
## looks good

# Step 1.2 plot decision 
a <- seq(-15,20,0.01)

f_c1.c2 <- function(lsdw,x){
  w <- lsdw[,1]-lsdw[,2]
  x <- cbind(1, x)
  return( x %*% w[1:2] / ( -1 * w[3]))
}

f_c2.c3 <- function(lsdw,x){
  w <- lsdw[,2]-lsdw[,3]
  x <- cbind(1, x)
  return( x %*% w[1:2] / ( -1 * w[3]))
}

f_c3.c1 <- function(lsdw,x){
  w <- lsdw[,3]-lsdw[,1]
  x <- cbind(1, x)
  return( x %*% w[1:2] / ( -1 * w[3]) )
}

# find the point 3 lines join together
ydiffs <- f_c3.c1(lsdw, a)-f_c2.c3(lsdw, a)
zeroPoint<- which.min(abs(ydiffs))

h <- a[1:zeroPoint]
t <- a[zeroPoint:length(a)]

points(t, f_c1.c2(lsdw, t), type='l', col = 'green')
points(t, f_c2.c3(lsdw, t), type='l', col = 'green')
points(h, f_c3.c1(lsdw, h), type='l', col='green')

# Step 1.3 classification function
  # now give the classification function, by using the weight we calculated
findClassLSD <- function(x, lsdw){
  x <- cbind(1, as.matrix(t(x)))
  y <- x %*% lsdw
  return(apply(y, 1, which.max))
}



# Step 2: Fisher 
# plot the data 
plot(Y~X,data=dat, col=Class)
x <- dat[,1:2]
class <- as.integer(dat[,3])

# weight matrix size: dim(x) * 1
# Step 2.1 sub funcions
# projection: y <- t(w) %*% x 
f.projection <- function(x,w){
  return(as.vector(w %*% t(x) / sqrt(sum(w^2))))
}

f.muK <- function(y,class){
  muK <- NULL
  for(c in unique(class)){
    muK <- c(muK, mean(y[class==c]))
  }
  return(muK)
}

f.mu <- function(y){
  return(mean(y))
}

f.Sw <- function(muK, y, class){
  result = 0
  for(i in 1:length(class)){
    result = result + (y[i]-muK[class[i]])^2
  }
  return(result)
}

f.Sb <- function(muK, mu, class){
  result = 0
  for(c in unique(class)){
    result = result + sum(class==c) * (muK[c] - mu)^2
  }
  return(result)
}

f.Jw <- function(Sw, Sb){
  return(sum((Sw^-1 * Sb)))
}

#Step 2.2: Define the function we want optimize 
# since optim will find min value, I multiple Jw by '-1' before return
optFun <- function(w){
  y <- f.projection(x, w)
  muK <- f.muK(y, class)
  mu <- f.mu(y)
  Sw <- f.Sw(muK, y, class)
  Sb <- f.Sb(muK, mu, class)
  return(-f.Jw(Sw, Sb))
}


# Step 2.3: optimze by optim function
optResult <- optim(par=t( c(1,-1) ), optFun)

# Step 2.4: draw the new space (line)
points(a, optResult$par[2]*a/optResult$par[1], type='l')

# Step 2.5: project data onto this space
w <- optResult$par
y <- f.projection(x, w)
muK <- f.muK(y, class)
y1 <- w[1]/sqrt(w[1]^2 + w[2]^2) * y
y2 <- w[2]/sqrt(w[1]^2 + w[2]^2) * y
points(y1,y2 , col = class, pch=4)

# Step 2.6 classification function:
findClassFisher <- function(x, w, muK){
  y <- f.projection(x,w)
  return(which.min(abs(muK-y)))
}


# Step 3. using both classification function:

findClassLSD(matrix(c(5,0,10,10,0,10), nrow=2,ncol=3), lsdw)


findClassFisher(t(c(5,0)), w, muK)
findClassFisher(t(c(10,10)), w, muK)
findClassFisher(t(c(0,10)), w, muK)



