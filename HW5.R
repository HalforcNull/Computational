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


# basis function
basis.fun <- function(x){
  phi_1 <- dmvnorm(x, mean=rep(0,2),sigma = diag(2)*0.5)
  phi_2 <- dmvnorm(x, mean=rep(-1,2), sigma = diag(2)*0.5)
  phi_3 <- dmvnorm(x, mean=rep(1,2), sigma = diag(2)*0.5)
  
  return(rbind(phi_1, phi_2, phi_3))
}

#4.59 4.87
sigma.fun <- function(w, Phi_X){
  a = t(w) %*% Phi_X
  return(as.vector(1/(1+exp(-a))))
}


# Create grid
library(plotly)
one.dim <- seq(-2,2,0.01)

x <- combn(one.dim, m = 2)
rev.x <-  rbind(x[2,],x[1,])
all.x <- cbind(x,rev.x)
all.x <- expand.grid(x = one.dim, y = one.dim)
all.phi.x <- basis.fun(all.x)

# Frequenist Logistic Reg.
# Find weight value and use it to do predict

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
target.T <- X[,3] - 1
N <- 3

w.old <- c(1,1,1)

for(i in 1:10000){
  # writeLines(paste0('Run ', i, '' ) )
  # writeLines('Old weight:')
  # writeLines(w.old)
  old.y <- sigma.fun(w.old, Phi_X)
  old.R <- matrix.R.fun(old.y)
  old.Z <- matrix.Z.fun(w.old, Phi_X, old.R, old.y, target.T)
  w.new <- as.vector( w.update(Phi_X, old.R, old.Z) )
  # writeLines('New weight:')
  # writeLines(w.new)
  if( abs(w.old[1] - w.new[1]) < 0.00000001 &
      abs(w.old[2] - w.new[2]) < 0.00000001 &
      abs(w.old[3] - w.new[3]) < 0.00000001 )
  {
    writeLines('Frequenist Logistic Reg. :')
    writeLines(paste0('find acceptable weight using ', i, ' Runs'))
    writeLines('New weight:')
    print(w.new)
    writeLines("")
    writeLines("")
    break
  }
  w.old <- w.new 
}

new.y <- sigma.fun(w.new, Phi_X)
pred.class <- ifelse(new.y < 0.5, 1, 2)
all.y <- sigma.fun(w.new, all.phi.x)
zmatrix <-  matrix( as.vector(all.y), nrow=length(one.dim), ncol=length(one.dim))
contour(x=one.dim, y=one.dim, z=zmatrix, levels=c(0.1,0.3,0.5,0.7,0.9), xlim = c(-2,2), ylim=c(-2,2), main = 'Frequenist Logistic Reg.')

points(X[X[,3]==1,1],X[X[,3]==1,2],pch=16,col="green")
points(X[X[,3]==2,1],X[X[,3]==2,2],pch=16,col="blue")
points(X[pred.class==1,1],X[pred.class==1,2],pch=0,col="green")
points(X[pred.class==2,1],X[pred.class==2,2],pch=0,col="blue")

writeLines('Frequenist Logistic Reg. Result:')
print(table(pred.class, target.T))
writeLines("")
writeLines("")



# Baysian Logistic Reg.
# 1. Find the distribution of weight
# 2. Do convolution (if possible, otherwise, sample from weight and calculate mean)


# Step 1. Laplace Approximation to approximate the distribution of weight 
# # W_new = w_old + t( t(Phi_X) * (t-y) * Sn_old) <- by first derivation = 0
# find.w.map <- function(w.old, Phi_X, Training.T, y.old, Sn.old){
#   #return( w.old + t( t(Phi_X) %*% (Training.T-y.old) %*% Sn_old) )
#   #         3*1        3*150       150 * 1                3*3
#   return(w.old + t( t(Training.T-y.old) %*% t(Phi_X) %*%  Sn.old))
#   #       3*1     t(  1 * 150                150*3      3*3   )
# }

# 4.142
ln.w.t <- function(w.map, m0, S0, t.Vector, Phi_X){
  y.Vector <- sigma.fun(w.map, Phi_X)
  part.1 = -0.5 * t(w.map - m0) %*% solve(S0) %*% (w.map-m0)
  part.2 = sum( t.Vector * log2(y.Vector) )
  part.3 = sum( (1-t.Vector) * log2(1-y.Vector))
  return(part.1+part.2+part.3)
}

# Sn_new = Sn_old^(-1) + t(Phi_X) * R * Phi_X
find.Sn <- function(Sn.old, Phi_X, Updated.R){
#  return(solve(Sn.old) + t(Phi_X) %*% Updated.R %*% Phi_X)
#         3*3             150*3        150*150       3*150
  return(solve(Sn.old) + Phi_X %*% Updated.R %*% t(Phi_X))
  #         3*3             150*3        150*150       3*150
}

# Run Iteration 
m0 <- c(1,1,1)
S0 <- diag(3)
Phi_X <- basis.fun(X[,-3])
target.T <- X[,3] - 1
N <- 3
init.w.map <- c(2,3,5)

opmtTarget <- function(w.map){
  return(-ln.w.t(w.map, m0, S0, target.T, Phi_X))
}

gradientFunc <- function(w.map){
  y.Vector <- sigma.fun(w.map, Phi_X)
  return((m0-w.map)%*%S0+(target.T-y.Vector)%*%t(Phi_X))
}

optimResult <- optim(init.w.map, opmtTarget, hessian=TRUE, gr=gradientFunc)
w.new <- optimResult$par

updated.y <- sigma.fun(w.new, Phi_X)
updated.R <- matrix.R.fun(updated.y)
Sn.new <- find.Sn(S0, Phi_X, updated.R)


####################### Compare hessian matrix given by optim function 
####################### And Sn calculated by 4.143
writeLines('Baysian Logistic Logistic Reg. :')
writeLines('')
writeLines('Weight Map :')
print(w.new)
writeLines('')
writeLines('Hessian v.s. Sn: ')
print(-optimResult$hessian)
print(Sn.new)
writeLines('')
writeLines('')

# Step 2. Find Covlution

## Solution 1: using Sampling
N = 100

baysianSamplingSolution <- function(Phi_X_new, w.map, Sn, N){
  sampled.w <- rmvnorm(N, w.map, Sn)
  sampled.y <- apply( sampled.w, 1, sigma.fun, Phi_X = Phi_X_new)
  return(rowMeans(sampled.y))
}


new.y <- baysianSamplingSolution(Phi_X, w.new, Sn.new, N)
pred.class <- ifelse(new.y < 0.5, 1, 2)


all.y <- baysianSamplingSolution(all.phi.x, w.new, Sn.new, N)
zmatrix <-  matrix( as.vector(all.y), nrow=length(one.dim), ncol=length(one.dim))
contour(x=one.dim, y=one.dim, z=zmatrix, levels=c(0.1,0.3,0.5,0.7,0.9),xlim=c(-2,2),ylim=c(-2,2), main = 'Baysian Logistic Reg. (Sampling weight)')

points(X[X[,3]==1,1],X[X[,3]==1,2],pch=16,col="green")
points(X[X[,3]==2,1],X[X[,3]==2,2],pch=16,col="blue")
points(X[pred.class==1,1],X[pred.class==1,2],pch=0,col="green")
points(X[pred.class==2,1],X[pred.class==2,2],pch=0,col="blue")

writeLines('Baysian Logistic Reg.(Sampling) Result:')
print(table(pred.class, target.T))
writeLines("")
writeLines("")

## Solution 2: using Approximation 
baysianApproximationSolution <- function(Phi_X_new, w.map, Sn){
  # 4.149
  mu.a <- t(w.new) %*% Phi_X_new
  # 4.150
  var.a <- t(Phi_X_new) %*% Sn %*% Phi_X_new
  # 4.154
  kappa <- 1/sqrt(1 + pi*var.a/8)
  return(1/(1+exp(-kappa %*% t(mu.a))))
}

new.y <- apply(Phi_X, 2, baysianApproximationSolution, w.map=w.new, Sn=Sn.new)
pred.class <- ifelse(new.y < 0.5, 1, 2)
all.y <- apply(all.phi.x, 2, baysianApproximationSolution, w.map=w.new, Sn=Sn.new)
zmatrix <-  matrix( as.vector(all.y), nrow=length(one.dim), ncol=length(one.dim))
dim(zmatrix)
contour(x=one.dim, y=one.dim, z=zmatrix, levels=c(0.1,0.3,0.5,0.7,0.9), xlim = c(-2,2), ylim=c(-2,2), main = 'Baysian Logistic Reg. (Kappa)')
points(X[X[,3]==1,1],X[X[,3]==1,2],pch=16,col="green")
points(X[X[,3]==2,1],X[X[,3]==2,2],pch=16,col="blue")
points(X[pred.class==1,1],X[pred.class==1,2],pch=0,col="green")
points(X[pred.class==2,1],X[pred.class==2,2],pch=0,col="blue")

writeLines('Baysian Logistic Reg.(Kappa) Result:')
print(table(pred.class, target.T))
writeLines("")
writeLines("")






# 
# pred.area <- ifelse( all.y < 0.5, 1, 2 )
# plot(all.x[pred.area==1,1], all.x[pred.area==1,2],col='purple', pch=16, xlim=c(-2,2),ylim=c(-2,2), main = 'Baysian Logistic Reg. (Kappa)')
# points(X[X[,3]==1,1],X[X[,3]==1,2],pch=16,col="green")
# points(X[X[,3]==2,1],X[X[,3]==2,2],pch=16,col="blue")
# points(X[pred.class==1,1],X[pred.class==1,2],pch=0,col="green")
# points(X[pred.class==2,1],X[pred.class==2,2],pch=0,col="blue")
# 
# 
# my.a <- t(w.new) %*% Phi_X
# plot(my.a[X[,3]==1], y = rep(0, sum(X[,3]==1)), pch=16,col="green", xlim=c(-15,15),ylim=c(-0.5,0.5), main = 'Data on decision space (one dim).')
# points(my.a[X[,3]==2], y = rep(0, sum(X[,3]==2)), pch=16,col="blue")
# 
# my.plot.ly.data <- cbind(t(Phi_X), X[,3] )
# colnames(my.plot.ly.data) <- c('x1','x2','x3','class')
# my.plot.ly.data <- as.data.frame(my.plot.ly.data)
# plot_ly(my.plot.ly.data, x = ~x1, y = ~x2, z = ~x3,
#         marker = list(color = ~class, colorscale = c('red','green', 'blue'), showscale = FALSE))


