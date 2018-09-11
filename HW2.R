n = 40        # sample size
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

# Step 2. Calculate Sn and Mn based on eq 3.53 3.54
# Remember we set M0 = 0, and S0^-1 = alpha * I
CalculateWeightSn <- function(alpha, beta, Phi_X, N){
  return( (alpha * diag(N) + beta * t(Phi_X) %*% Phi_X)^-1 )
}

CalculateWeightMn <- function(beta, Sn, Phi_X, t){
  return( beta * Sn %*% t(Phi_X) %*% t)
}

# Step 3. Calculate mean and variance of t given x, based on 3.58 and 3.59
CalculateTMn <- function(WMn, Phi_X){
  return(t(WMn) %*% t(Phi_X))
}

CalculateTVar <- function(WSn, Phi_X, beta){
  return(beta^-1 + t(Phi_X) %*% WSn %*% Phi_X)
}



# Step 4. using above functions draw graphs
mu <- seq(0,1,0.125)
sd <- rep.int(x=0.1, times=N)

Phi_X <- applyBasicFuntion(x, mu, sd, N)

alpha <- 2
beta <- 100

W$Sn <- CalculateWeightSn(alpha, beta, Phi_X, N)
W$Mn <- CalculateWeightMn(beta, W$Sn, Phi_X, t)

a <- seq(0,1,0.01)
Phi_A <- applyBasicFuntion(a, mu, sd, N)

Predict$Mn <- CalculateTMn(W$Mn, Phi_A)
Predict$Var <- CalculateTVar(W$Sn, Phi_A, beta)

plot(a, t$Mn, ylim = c(-2,2), type='l')

library('mvtnorm')





# # Step 3. function: Basis functions
# Predictions <- function(x, W, mu, sd, N, WMax = NULL, WMin = NULL){
#   if(!is.null(WMax) && !is.null(WMin) 
#      && length(WMax) == N && length(WMin) == N ){
#     CalcMaxAndMin = TRUE
#   }else{
#     CalcMaxAndMin = FALSE
#   }
#   
#   W <- diag(W)
#   WMax <- diag(WMax)
#   WMin <- diag(WMin)
#   t$predic <- applyBasicFuntion(x,mu,sd,N) %*% W
#   if(CalcMaxAndMin){
#     tmp_1 <- applyBasicFuntion(x,mu,sd,N) %*% WMax
#     tmp_2 <- applyBasicFuntion(x,mu,sd,N) %*% WMin
#     t$predicMax <- pmax(tmp_1,tmp_2)
#     t$predicMin <- pmin(tmp_1,tmp_2)
#   }
#   else
#   {
#     t$predicMax <- NULL
#     t$predicMin <- NULL
#   }
#   return(t)
# }

# Step 4. using above functions draw graphs
mu <- seq(0,1,0.125)
sd <- rep.int(x=0.1, times=9)

Phi_X <- applyBasicFuntion(x, mu, sd, 9)
W <- solveW(Phi_X, t)

# 
# I tried alpha = 2, it is toooooo small, Set alpha = 2 which means W ~ N(mu=W, var=0.5) ????????????????????
# Why book using alpha = 2???????????????
#
# In my case I use alpha = 16

# Then find the max and min of each W (CI of 80%)
# Then the prediction range = function(x, Wmax, mu, sd, N) and function(x, Wmin, mu, sd, N)
alpha = 4
offset <- qnorm(p = 0.9,sd=sqrt(alpha^-1))
W_Max <- W * ( 1 + offset )
W_Min <- W * ( 1 - offset )


a <- seq(0,1,0.01)

notMergedT <- Predictions(a, W, mu, sd, N, WMax=W_Max, WMin=W_Min )

MergedMaxPredict <- rowSums(notMergedT$predicMax)
MergedMinPredict <- rowSums(notMergedT$predicMin)


points(a, MergedMinPredict, ylim = c(-2,2), type = 'l')

polygon(c(a, rev(a)), c(MergedMaxPredict, rev(MergedMinPredict)), col='skyblue')
points(x,t,ylim=c(-4,4))


TODO: NOT DOne Yet.


 mergedMax <- rowSums(Pre_Max_by_basis)
plot(a,(Pre_Min),type='l', ylim = c(-4,4))
par(new=TRUE)
plot(a,rowSums(Pre_Max),type='l', ylim = c(-4,4))
polygon(rowSums(Pre_Min), rowSums(Pre_Max), col='blue')





