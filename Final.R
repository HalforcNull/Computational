
load("C:/GitRepo/Computational/ink.training.rdata")

x1 <- NULL
x2 <- NULL
for(i in 1:22){
    tmp.x1 <- as.numeric(ink.training.dat[,,i,1]) 
    tmp.x2 <- as.numeric(ink.training.dat[,,i,2])
    x1 <- rbind(x1, tmp.x1)
    x2 <- rbind(x2, tmp.x2)
}

dat.x <- rbind(x1,x2)
dat.y <- c(rep(1,22), rep(0,22))


# kenrel
# kenrel = exp[-norm2(x,y)^2]
ker <- function(x,y){
    return(exp(-0.5 * 16 * 1.414 * (x-y) * (x-y)))
}








