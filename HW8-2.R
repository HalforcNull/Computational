
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


# Changes from HW7
#   1. aN. In this case, we have a 2*n length vector, to contain c(aN, aN_hat)
#   2. Dmat. D_mat is changed since the quadratic part of the problem is changed
#   3. dvect. dvect changed as we discussed during class
#   4. Amat and bvect. Since the constrain of the problem is changed
#   5. Kernel function. ker_1 = abs(x1-x2)
#                       ker_2 = min( abs( (x1-0.25)^2 - (x2-0.25)^2 ), abs( (x1-0.75)^2 - (x2-0.75)^2) ))




