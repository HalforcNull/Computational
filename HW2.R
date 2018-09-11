
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

# Step 2. 
solveW <- function(Phi_X, t){
  return(qr.solve(Phi_X, t))
}
# Calculate Sn and Mn based on eq 3.53 3.54
# Remember we set M0 = 0, and S0^-1 = alpha * I
CalculateWeightSn <- function(alpha, beta, Phi_X, N){
  return( solve(alpha * diag(N) + beta * t(Phi_X) %*% Phi_X) )
}

CalculateWeightMn <- function(beta, Sn, Phi_X, t){
  return( beta * Sn %*% t(Phi_X) %*% t)
}

# Step 3. Calculate mean and variance of t given x, based on 3.58 and 3.59
CalculateTMn <- function(WMn, Phi_X){
  return(t(WMn) %*% t(Phi_X))
}

CalculateTVar <- function(WSn, Phi_X, beta){
  myVar <- NULL
  for(i in 1:nrow(Phi_X)){
    myVar <- cbind(myVar,beta^-1 + t(Phi_X[i,]) %*% WSn %*% Phi_X[i,])
  }
    
  return(myVar)
}



# Step 4. using above functions draw graphs
drawFigure <- function(x, t, mu, sd, N, alpha, beta)
{
  Phi_X <- applyBasicFuntion(x, mu, sd, N)

  W <- solveW(Phi_X, t)
  
  W$Sn <- CalculateWeightSn(alpha, beta, Phi_X, N)
  W$Mn <- CalculateWeightMn(beta, W$Sn, Phi_X, t)
  
  a <- seq(0,1,0.01)
  Phi_A <- applyBasicFuntion(a, mu, sd, N)
  
  Pred_Mn <- CalculateTMn(W$Mn, Phi_A)[1,]
  Pred_Var <- CalculateTVar(W$Sn, Phi_A, beta)[1,]
  
  
  PredMax <- Pred_Mn + sqrt(Pred_Var)
  PredMin <- Pred_Mn - sqrt(Pred_Var)
  
  plot(a, Pred_Mn, xlim=c(0,1), ylim = c(-2,2), type='l', main = paste0('Result of training size = ', as.character(length(x))))
  points(a, sin(2*pi*a), xlim=c(0,1), ylim=c(-2,2), type='l', col='red')
  points(a, PredMax, xlim=c(0,1), ylim = c(-2,2), type='l')
  points(x,t)
  points(a, PredMin, xlim=c(0,1), ylim=c(-2,2), type='l')
  polygon(c(a, rev(a)), c(PredMax, rev(PredMin)), col=rgb(0,1,1,0.4))
}

# Step 5. Run test given different size of training data

# Sampling from sin(x)
n = 300        # sample size (max)
N = 9         # Module size ( how many normal distribution is contained)

x <- runif(n) # sample x using uniform distribution. By default, min = 0, max = 1.
t <- sin(2*pi*x) + rnorm(n, sd = 0.1) # calculate t, plug gaussian noise


# setup basis function mean and sd
mu <- seq(0,1,0.125)
sd <- rep.int(x=0.1, times=N)

# setup alpha beta
alpha <- 2
beta <- 25

# n = 1
datasetX <- x[1:1]
datasetT <- t[1:1]
drawFigure(datasetX, datasetT, mu, sd, N, alpha, beta)


# n = 2
datasetX <- x[1:2]
datasetT <- t[1:2]
drawFigure(datasetX, datasetT, mu, sd, N, alpha, beta)

# n = 4
datasetX <- x[1:4]
datasetT <- t[1:4]
drawFigure(datasetX, datasetT, mu, sd, N, alpha, beta)

# n = 20
datasetX <- x[1:20]
datasetT <- t[1:20]
drawFigure(datasetX, datasetT, mu, sd, N, alpha, beta)


