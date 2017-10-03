data(normmix.data)
L   <- normmix.data$L
k   <- ncol(L)
x0  <- rep(k)/k
out <- mixsqp(L,x0)
