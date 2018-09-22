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


freq.predict <- ifelse(freq.prob > 0.5, 1, 2)
Y <- cbind(X[,-3], freq.predict)

plot(X[X[,3]==1,1],X[X[,3]==1,2],pch=16,col="green",xlim=c(-2,2),ylim=c(-2,2), main = 'Use glm() function directly')
points(X[X[,3]==2,1],X[X[,3]==2,2],pch=16,col="blue")
points(Y[Y[,3]==1,1],Y[Y[,3]==1,2],pch=0,col="green")
points(Y[Y[,3]==2,1],Y[Y[,3]==2,2],pch=0,col="blue")



