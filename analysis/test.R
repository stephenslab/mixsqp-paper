# Small script for testing Rcpp implementation of mixsqp method.
#
# Commands to update Rcpp code and install package:
#
#   devtools::document()
#   devtools::install_local("mixopt")
#
library(Rcpp)
sourceCpp("src/mixsqp.cpp")
source("R/mixopt.R")
set.seed(1)

# Load the data.
load("data/normmix.data.RData")
L <- normmix.data$L
k <- ncol(L)
n <- nrow(L)

# Set initial solution.
x <- runif(k)
x <- x/sum(x)

# Run the Rcpp implementation.
out <- mixsqp(L,x = x,pqrtol = 0,algorithm.version = "Rcpp")

