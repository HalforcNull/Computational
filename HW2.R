n = 1        # sample size
N = 9         # Module size ( how many normal distribution is contained)

# Step 0. Sampling from sin(x)
x <- runif(n) # sample x using uniform distribution. By default, min = 0, max = 1.
t <- sin(2*pi*x) + rnorm(n, sd = 0.1) # calculate t, plug gaussian noise



# Step 1. function: apply basis functions on X
applyBasicFuntion <- function(x, mu, sd, N){
  if(length(mu)!=N || length(sd)!=N){
    return(NULL)
  }
  
  Phi_X = NULL
  for(i in 1:N){
    Phi_X <- cbind(Phi_X, dnorm(x, mean=mu[i],sd=sd[i]))
  }
  
  return(as.matrix(Phi_X))
}

# Step 2. funciton: solve problem as we used in HW1 : t(w) %*% t(Phi_X) = t(T)
solveW <- function(Phi_X, t){
  return(qr.solve(Phi_X, t))
}

# Step 3. function: Basis functions
Predictions <- function(x, W, mu, sd, N, WMax = NULL, WMin = NULL){
  t$predic <- W %*% applyBasicFuntion(x,mu,sd,N)
  if(!is.null(WMax) && !is.null(WMin) 
     && length(WMax) == N && length(WMin) == N ){
    tmp_1 <- WMax %*% applyBasicFuntion(x,mu,sd,N)
    tmp_2 <- WMin %*% applyBasicFuntion(x,mu,sd,N)
    t$predicMax <- pmax(tmp_1,tmp_2)
    t$predicMin <- pmin(tmp_1,tmp_2)
  }
  else
  {
    t$predicMax <- NULL
    t$predicMin <- NULL
  }
  return(t)
}

# Step 4. using above functions draw graphs
mu <- seq(0,1,0.125)
sd <- rep.int(x=0.5, times=9)

Phi_X <- applyBasicFuntion(x, mu, sd, 9)
W <- solveW(Phi_X, t)

# 
# I tried alpha = 2, it is toooooo small, Set alpha = 2 which means W ~ N(mu=W, var=0.5) ????????????????????
# Why book using alpha = 2???????????????
#
# In my case I use alpha = 16

# Then find the max and min of each W (CI of 80%)
# Then the prediction range = function(x, Wmax, mu, sd, N) and function(x, Wmin, mu, sd, N)
alpha = 64
offset <- qnorm(p = 0.9,sd=sqrt(alpha^-1))
W_Max <- W + offset
W_Min <- W - offset


a <- seq(0,1,0.01)

notMergedT <- Predictions(a, W, mu, sd, N, WMax=W_Max, WMin=W_Min )

TODO: NOT DOne Yet.


mergedMax <- rowSums(Pre_Max_by_basis)
plot(a,(Pre_Min),type='l', ylim = c(-4,4))
par(new=TRUE)
plot(a,rowSums(Pre_Max),type='l', ylim = c(-4,4))
polygon(rowSums(Pre_Min), rowSums(Pre_Max), col='blue')





