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






