# TO DO: Explain what this function does, and how to use it.
#
# NOTES:
#
#   - Need to add CVXR dependency to DESCRIPTION.
#
library(CVXR)
load("../../../data/normmix.data.RData")
L <- normmix.data$L
L <- L[,18:20]
k <- ncol(L)
x <- Variable(k)
objective   <- Maximize(sum(log(L %*% x)))
constraints <- list(x > 0,sum(x) == 1)
prob        <- Problem(objective,constraints)
out         <- solve(prob)
out$getValue(x)
