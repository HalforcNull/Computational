remove(list=ls())
### data for Least square multiclass discriminant function and Fisher multiclass discriminant function
library(mvtnorm)
library(scales)

# define number of objects in each class
N <- 40

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
### least square multiclass discriminant function
# prepare the data by adding a column of 1 and changing the class indicator
dat <- data.frame(intercept = 1, dat)
T.mat <- t(sapply(dat$Class, function(x){T.vect <- rep(0,3); T.vect[x]<-1; return(T.vect)}))
X.mat <- as.matrix(dat[,1:3])
W.tilda.LSD <- solve(t(X.mat)%*%X.mat)%*%t(X.mat)%*%T.mat

# get the prediction
dat$pred.LSD <- apply(X.mat%*%W.tilda.LSD,1,which.max)

#plot
quartz()
plot(Y~X,data=dat, col=Class, pch=20)
points(Y~X,data=dat, col=pred.LSD, pch=2,lwd=1)

### get the boundaries
xs <- seq(min(dat$X),max(dat$X),length.out=500)
ys <- seq(min(dat$Y),max(dat$Y),length.out=500)

# boundary between C1 and C2
tmp <- W.tilda.LSD[,2]-W.tilda.LSD[,1]
tmp2 <- W.tilda.LSD[,3]-W.tilda.LSD[,1]
C12.y.bound <- -1*tmp[2]/tmp[3]*xs - tmp[1]/tmp[3]
C12.y.constrain <- -1*tmp2[2]/tmp2[3]*xs - tmp2[1]/tmp2[3]
ix <- which(C12.y.bound<C12.y.constrain)
lines(xs[ix],C12.y.bound[ix])

# boundary between C1 and C3
tmp <- W.tilda.LSD[,3]-W.tilda.LSD[,1]
tmp2 <- W.tilda.LSD[,2]-W.tilda.LSD[,1]
C13.y.bound <- -1*tmp[2]/tmp[3]*xs - tmp[1]/tmp[3]
C13.y.constrain <- -1*tmp2[2]/tmp2[3]*xs - tmp2[1]/tmp2[3]
ix <- which(C13.y.bound<C13.y.constrain)
lines(xs[ix],C13.y.bound[ix])

# boundary between C2 and C3
tmp <- W.tilda.LSD[,3]-W.tilda.LSD[,2]
tmp2 <- W.tilda.LSD[,1]-W.tilda.LSD[,2]
C23.y.bound <- -1*tmp[2]/tmp[3]*xs - tmp[1]/tmp[3]
C23.y.constrain <- -1*tmp2[2]/tmp2[3]*xs - tmp2[1]/tmp2[3]
ix <- which(C23.y.bound>C23.y.constrain)
lines(xs[ix],C23.y.bound[ix])


### Fisher discriminant function
library(Matrix)
# prepare the data by adding a column of 1 and changing the class indicator
# create a data.frame of the data
dat <- data.frame(X = c(c1.dat[,1],c2.dat[,1],c3.dat[,1]), Y = c(c1.dat[,2],c2.dat[,2],c3.dat[,2]), Class = as.factor(rep(c(1,2,3),each=N)))
dat <- data.frame(intercept = 1, dat)

# get the mean vectors for the classes
mean.vects <- by(dat[,2:3],dat$Class, colMeans)
# get the grand mean
grand.mean.vect <- colMeans(dat[,2:3])

# get the within covariance
S.within <- Reduce('+',lapply(1:3, function(x, dat, mean.vects){ (t(dat[dat$Class==x,2:3])-mean.vects[[x]])%*%t(t(dat[dat$Class==x,2:3])-mean.vects[[x]])}, dat=dat, mean.vects = mean.vects))

# get the between covariance
S.between <- Reduce('+',lapply(1:3, function(x, mean.vects, grand.mean.vect){ (mean.vects[[x]]-grand.mean.vect)%*%t(mean.vects[[x]]-grand.mean.vect)}, mean.vects = mean.vects, grand.mean.vect=grand.mean.vect)) * N

S.total <- (t(dat[,2:3])-grand.mean.vect)%*%t(t(dat[,2:3])-grand.mean.vect)
S.within + S.between

# calculate criteria
eigen.vectors <- eigen(solve(nearPD(S.within)$mat)%*%(S.between))$vect

# get the prediction
dat$pred.FISHER <- as.matrix(dat[,2:3])%*%eigen.vectors
plot(dat$pred.FISHER[,1],dat$pred.FISHER[,2],col=dat$Class,pch=20)

# project the means
proj.mean.vects <- lapply(mean.vects, function(x,eigen.vectors){x%*%eigen.vectors},eigen.vectors=eigen.vectors)
proj.mean.vects <- matrix(unlist(proj.mean.vects),3,2,byrow=TRUE)
points(proj.mean.vects[,1],proj.mean.vects[,2],col=1:3,pch=11)

pred.labs <- apply(dat$pred.FISHER, 1, function(x, proj.mean.vects){which.min(diag(t(x-proj.mean.vects)%*%(x-proj.mean.vects)))}, proj.mean.vects=t(proj.mean.vects))
points(dat$pred.FISHER[,1],dat$pred.FISHER[,2],col=pred.labs,lwd=2)

table(pred.labs, dat$Class)





# initialise W.tilda
#W.tilda.init <- t(W.tilda.LSD[2:3,])
W.tilda.init <- matrix(runif(6,0,1),3,2)
# create the objective function
objective.fun <- function(par, dat,mean.vects, grand.mean.vect, N){
   W.tilda <- matrix(par,3,2,byrow=FALSE)
   # get the within covariance
   S.within <- Reduce('+',lapply(1:3, function(x, dat, mean.vects){ (t(dat[dat$Class==x,2:3])-mean.vects[[x]])%*%t(t(dat[dat$Class==x,2:3])-mean.vects[[x]])}, dat=dat, mean.vects = mean.vects))
   # get the between covariance
   S.between <- Reduce('+',lapply(1:3, function(x, mean.vects, grand.mean.vect){ (mean.vects[[x]]-grand.mean.vect)%*%t(mean.vects[[x]]-grand.mean.vect)}, mean.vects = mean.vects, grand.mean.vect=grand.mean.vect)) * N
   #J <- sum(diag(as.matrix(solve(nearPD(W.tilda%*%S.within%*%t(W.tilda))$mat)%*%(W.tilda%*%S.between%*%t(W.tilda)))))
   J <- sum(diag(as.matrix(solve(W.tilda%*%S.within%*%t(W.tilda)+diag(3)*10)%*%(W.tilda%*%S.between%*%t(W.tilda)))))
   print(c(S.within,S.between,J))
   return(J*-1)
}

param <- optim(par = as.vector(W.tilda.init), fn =objective.fun,  dat=dat, mean.vects=mean.vects, grand.mean.vect= grand.mean.vect,N=N, method='BFGS',control = list(maxit=1000))
W.optim <- t(matrix(param$par,3,2,byrow=FALSE))

# get the prediction
dat$pred.FISHER.optim <- as.matrix(dat[,2:3])%*%W.optim
quartz()
plot(dat$pred.FISHER.optim[,1],dat$pred.FISHER.optim[,2],col=dat$Class, pch=20)

# project the means
proj.mean.vects.optim <- lapply(mean.vects, function(x,W.optim){x%*%W.optim},W.optim=W.optim)
proj.mean.vects.optim <- matrix(unlist(proj.mean.vects.optim),3,3,byrow=TRUE)[,1:2]
points(proj.mean.vects.optim[,1],proj.mean.vects.optim[,2],col=1:3,pch=11)

pred.labs <- apply(dat$pred.FISHER.optim[,1:2], 1, function(x, tmp){which.min(diag(t(x-tmp)%*%(x-tmp)))}, tmp=t(proj.mean.vects.optim))
points(dat$pred.FISHER.optim[,1],dat$pred.FISHER.optim[,2],col=pred.labs,lwd=2)

table(pred.labs, dat$Class)
