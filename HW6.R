n = 10       # sample size


# Step 1. Sampling from sin(x)
x <- runif(n) # sample x using uniform distribution. By default, min = 0, max = 1.
t <- sin(2*pi*x) + rnorm(n, sd = 0.1) # calculate t

t.Mean <- mean(t)
t <- t-t.Mean


# Step 2. Kernel Function : k(x, x_n) = { pdf of N(x_n, sd^2) @ x } / { sum  of all pdf}
g.fun <- function(x, xn, hyperSd){
  return(dnorm(x,xn,hyperSd))
}

kernelFun <- function(gfun.list){
  return(gfun.list/sum(gfun.list))
}

# Step 3. E(t) = sum( k(x,xn) * tn ), variance = sum( k(x,xn) * tn^2 ) - [E(t)]^2
calcMean <- function(k_the_pdf, training_t){
  return(sum(k_the_pdf * training_t))
}

calcVar <- function(k_the_pdf, training_t, hyperSd){
  new_t_mean <- calcMean(k_the_pdf, training_t)
  return( sum( k_the_pdf * ( training_t^2 + hyperSd ^ 2) ) - new_t_mean^2 )
}


# Step 4. calculate
hyperSd = 0.1
a <- seq(0, 1, 0.01)

gfun.matrix <- apply( as.matrix(a),1, g.fun, xn = x, hyperSd = hyperSd)

k_the_pdf <- apply( gfun.matrix,2, kernelFun )
predictMean <- apply( k_the_pdf, 2, calcMean, training_t = t) + t.Mean
predictVar <- apply( k_the_pdf, 2, calcVar, training_t = t, hyperSd = hyperSd)
predictSd <- sqrt(predictVar)
PredMax <- predictMean + 1.95 * predictSd
PredMin <- predictMean - 1.95 * predictSd

# Step 5. Plot 
plot(a, sin(2 * pi * a), type="l", col = 'red', xlim=range(c(0,1)), ylim=range(c(-2,2)))
points(x, t)

points(a, predictMean, type="l", col = 'green')

polygon(c(a, rev(a)), c(PredMax, rev(PredMin)), col=rgb(0,1,1,0.4))

