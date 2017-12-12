# TO DO: Explain what this function does, and how to use it.
#
# NOTES:
#
#   - Need to add CVXR dependency to DESCRIPTION.
#
library(CVXR)
load("../../../data/normmix.data.RData")
L <- normmix.data$L
k <- ncol(L)
x <- Variable(k)
objective   <- Maximize(sum(log(L %*% x)))
constraints <- list(x > 0,sum(x) == 1)
prob        <- Problem(objective,constraints)
out         <- solve(prob,solver = "SCS")
out$getValue(x)

n = dim(L)[1]; k = dim(L)[2];
t_rmosek = system.time(res <- REBayes::KWDual(L, rep(1,k), rep(1,n)/n))[3]

set.seed(2017)
n = 1000
m = 1.2
n = floor(n/2)
z = c(rnorm(n),4*rt(n,df=6))
z = z[order(abs(z))]
data = ashr::set_data(z,1)
grid = ashr:::autoselect.mixsd(data, mult=m, mode=0,
    mixcompdist = "normal", grange = c(-Inf,Inf))
grid = c(0,grid)
k = length(grid)
g  = ashr::normalmix(rep(1/k,k),rep(0,k),grid)
llik <- t(ashr:::log_comp_dens_conv.normalmix(g,data))
L = llik - apply(llik, 1, max)
L = exp(L)

print(sum(log(L %*% y)),digits=12)
