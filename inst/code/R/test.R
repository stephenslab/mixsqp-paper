library(mixopt)
data(normmix.data)
L   <- normmix.data$L
k   <- ncol(L)
x0  <- rep(1,k)/k
out <- mixsqp(L,x0)
