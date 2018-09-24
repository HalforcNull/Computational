# STAT721
# Data for homework Chapter 4
library(mvtnorm)
library(rgl)
library(Matrix)

set.seed(7675)
x1.tmp <- cbind(rmvnorm(70,c(0,0),matrix(c(0.2,0,0,0.1),2,2)),1)
x2.tmp <- cbind(rmvnorm(40,c(1,1),matrix(c(0.2,0,0,0.1),2,2)),2)
x3.tmp <- cbind(rmvnorm(40,c(-1,-1),matrix(c(0.2,0,0,0.1),2,2)),2)

X <- rbind(x1.tmp,x2.tmp,x3.tmp)

plot(X[X[,3]==1,1],X[X[,3]==1,2],pch=16,col="green",xlim=c(-2,2),ylim=c(-2,2))
points(X[X[,3]==2,1],X[X[,3]==2,2],pch=16,col="blue")



# Frequenist Logistic Reg.
my.df <- as.data.frame(X)
colnames(my.df) <- c('x.cod','y.cod','grp')
my.df$grp <- my.df$grp - 1
model.fit <- glm(grp~x.cod+y.cod, data=my.df, family = binomial)

my.newData <- as.data.frame(X[,-3])
colnames(my.newData) <- c('x.cod','y.cod')
freq.prob <- predict(model.fit,
                     newdata = my.newData,
                     type='response')


freq.predict <- ifelse(freq.prob > 0.5, 1, 2)
Y <- cbind(my.newData, freq.predict)

plot(X[X[,3]==1,1],X[X[,3]==1,2],pch=16,col="green",xlim=c(-2,2),ylim=c(-2,2), main = 'Use glm() function directly')
points(X[X[,3]==2,1],X[X[,3]==2,2],pch=16,col="blue")
points(Y[Y[,3]==1,1],Y[Y[,3]==1,2],pch=0,col="green")
points(Y[Y[,3]==2,1],Y[Y[,3]==2,2],pch=0,col="blue")

table(Y[,3], X[,3])

# Apply some 'kernel' function 
my.df <- as.data.frame(X)
colnames(my.df) <- c('x.cod','y.cod','grp')
my.df$new.cod <- my.df$x.cod ^ 2 + my.df$y.cod^2
my.df$grp <- my.df$grp - 1
model.fit <- glm(grp~new.cod, data=my.df, family = binomial)

my.newData <- as.data.frame(X[,1]^2 + X[,2]^2)
colnames(my.newData) <- c('new.cod')
freq.prob <- predict(model.fit,
                     newdata = my.newData,
                     type='response')


freq.predict <- ifelse(freq.prob < 0.5, 1, 2)
Y <- cbind(X[,-3], freq.predict)

plot(X[X[,3]==1,1],X[X[,3]==1,2],pch=16,col="green",xlim=c(-2,2),ylim=c(-2,2), main = 'Use glm() function directly')
points(X[X[,3]==2,1],X[X[,3]==2,2],pch=16,col="blue")
points(Y[Y[,3]==1,1],Y[Y[,3]==1,2],pch=0,col="green")
points(Y[Y[,3]==2,1],Y[Y[,3]==2,2],pch=0,col="blue")

table(Y[,3], X[,3])

# Baysian solution

# basis function
basis.fun <- function(x){
  phi_1 <- dmvnorm(x, mean=rep(0,2),sigma = diag(2))
  phi_2 <- dmvnorm(x, mean=rep(-1,2), sigma = diag(2))
  
  return(rbind(phi_1, phi_2))
}

#4.59 4.87
sigma.fun <- function(w, Phi_X){
  a = t(w) %*% Phi_X
  return(as.vector(1/(1+exp(-a))))
}

#4.98
matrix.R.fun <- function(y){
  return(diag(y*(1-y)))
}

#4.100
matrix.Z.fun <- function(w, Phi_X, matrix.R, y, t){
  return(t(w %*% Phi_X)- solve(matrix.R) %*% (y-t))
}

# update w 4.99
w.update <- function(Phi_X, matrix.R, matrix.Z){
  return(solve(Phi_X%*%matrix.R%*%t(Phi_X)) %*% Phi_X %*% matrix.R %*% matrix.Z)
  #             2*150   150*150    150*2         2*150         150*150     150*1  
}

# run iteration
Phi_X <- basis.fun(X[,-3])
target.T <- X[,3] -1
N <- 2

w.old <- c(1,1)

for(i in 1:10000){
  print(paste0('Run ', i, '' ) )
  print('Old weight:')
  print(w.old)
  old.y <- sigma.fun(w.old, Phi_X)
  old.R <- matrix.R.fun(old.y)
  old.Z <- matrix.Z.fun(w.old, Phi_X, old.R, old.y, target.T)
  w.new <- as.vector( w.update(Phi_X, old.R, old.Z) )
  print('New weight:')
  print(w.new)
  if( abs(w.old[1] - w.new[1]) < 0.00000001 &
      abs(w.old[2] - w.new[2]) < 0.00000001 )
  {
    print('find acceptable weight')
    break
  }
  w.old <- w.new 
}

new.y <- sigma.fun(w.new, Phi_X)
pred.class <- ifelse(new.y < 2, 1, 2)
plot(X[X[,3]==1,1],X[X[,3]==1,2],pch=16,col="green",xlim=c(-2,2),ylim=c(-2,2), main = 'Use glm() function directly')
points(X[X[,3]==2,1],X[X[,3]==2,2],pch=16,col="blue")
points(X[pred.class==1,1],X[pred.class==1,2],pch=0,col="green")
points(X[pred.class==2,1],X[pred.class==2,2],pch=0,col="blue")


###  Try other basis function
{
basis.fun2 <- function(x){
  phi_1 <- dmvnorm(x, mean=rep(0,2), sigma = diag(c(0.2, 0.1)))
  phi_2 <- dmvnorm(x, mean=rep(-1,2), sigma = diag(c(0.2, 0.1)))
  return(rbind(phi_1, phi_2))
}

Phi_X <- basis.fun2(X[,-3])
target.T <- X[,3]

w.old <- c(1,1)

for(i in 1:10000){
  print(paste0('Run ', i, '' ) )
  print('Old weight:')
  print(w.old)
  old.y <- sigma.fun(w.old, Phi_X)
  old.R <- matrix.R.fun(old.y)
  old.Z <- matrix.Z.fun(w.old, Phi_X, old.R, old.y, target.T)
  w.new <- as.vector( w.update(Phi_X, old.R, old.Z) )
  print('New weight:')
  print(w.new)
  if( abs(w.old[1] - w.new[1]) < 0.00000001 &
      abs(w.old[2] - w.new[2]) < 0.00000001 )
  {
    print('find acceptable weight')
    break
  }
  w.old <- w.new 
}

new.y <- sigma.fun(w.new, Phi_X)
pred.class <- ifelse(new.y < 2, 1, 2)
plot(X[X[,3]==1,1],X[X[,3]==1,2],pch=16,col="green",xlim=c(-2,2),ylim=c(-2,2), main = 'Use glm() function directly')
points(X[X[,3]==2,1],X[X[,3]==2,2],pch=16,col="blue")
points(X[pred.class==1,1],X[pred.class==1,2],pch=0,col="green")
points(X[pred.class==2,1],X[pred.class==2,2],pch=0,col="blue")

}
