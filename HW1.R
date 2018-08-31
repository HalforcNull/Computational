n = 10        # sample size


# Step 1. Sampling from sin(x)
x <- runif(n) # sample x using uniform distribution. By default, min = 0, max = 1.
t <- sin(2*pi*x) + rnorm(n, sd = 0.1) # calculate t



findW <- function(m){
  # Step 2. Build X matrix
  X <- rep.int(1, n)
  for(i in 1:(m-1)){
    X <- cbind(X, x^i)
  }
  X <- as.matrix(X)
  
  
  # Step 3. By calculation, we need solve t(w) %*% t(X) = t(T)
  w <- qr.solve(X, t)
  return(w)
}


# m = 1
m = 1
a <- seq(0, 1, 0.01)
plot(a, sin(2 * pi * a), type="l", col = 'red', xlim=range(c(0,1)), ylim=range(c(-1.5,1.5)))
points(x, t)
y_1 = findW(m)[1]
abline(y_1)


# function m > 1
drawModule <- function(m){
  a <- seq(0, 1, 0.01)
  plot(a, sin(2 * pi * a), type="l", col = 'red', xlim=range(c(0,1)), ylim=range(c(-1.5,1.5)))
  points(x, t)
  weight <- findW(m)
  y <- 0
  for(i in 1:m) {
    y = y + weight[i] * a ^i
  }
  par(new=TRUE)
  plot(a, y, type="l", col = 'green',xlim=range(c(0,1)), ylim=range(c(-1.5,1.5)))
}

# m = 2 3 5 7
drawModule(2)
drawModule(3)
drawModule(5)
drawModule(7)

