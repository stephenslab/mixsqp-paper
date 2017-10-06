library(rjulia)
library(mixopt)
julia_init()
data(normmix.data)
L   <- normmix.data$L
k   <- ncol(L)
x0  <- rep(1,k)/k
out <- mixsqp(L,x0,verbose = FALSE)
