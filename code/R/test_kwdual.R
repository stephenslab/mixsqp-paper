# Small script to show how the MOSEK interior-point solver runtime can
# be assessed more precisely using the output from Rmosek::mosek.
# 
# To use this script, the KWDual function in the REBayes package needs
# to be modified to output the Rmosek problem specification. This can
# be done by changing the last line of the KWDual function definition
# in KWDual.R to
#
#   list(P = P,f = f, g = g, status = status)
#
# then installing REBayes locally from source, e.g.,
#
#   devtools::install_local("REBayes")
#
library(Rmosek)
library(REBayes)

# Load the conditional likelihood matrix.
L <- as.matrix(read.table("../../data/sample10000x40.txt.gz"))
colnames(L) <- NULL
n <- nrow(L)
k <- ncol(L)

# Solve using MOSEK.
out.kwdual <- KWDual(L,rep(1,k),rep(1,n)/n)
r <- system.time(out.mosek <- mosek(out.kwdual$P,
                                    opts = list(getinfo = TRUE,verbose = 0)))
cat("Total time including R overhead:",r["elapsed"],"\n")
cat("Time spent in MOSEK solver:     ",out.mosek$dinfo$OPTIMIZER_TIME,"\n")
