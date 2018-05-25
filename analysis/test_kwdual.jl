# Small script to check overhead of calling MOSEK solver from Rmosek
# which is called from Julia. Compare the first "elasped" time against
# the OPTIMIZER_TIME output from calling Rmosek::mosek.
using Distributions
using LowRankApprox
using RCall
include("../code/julia/datasim.jl");
include("../code/julia/likelihood.jl");

srand(1);

x = normtmixdatasim(round(Int,25e4));

sd = autoselectmixsd(x,nv = 100);
L  = normlikmatrix(x,sd = sd);
size(L)

function REBayes(L)
  @rput L;
  R"library(REBayes)
  library(Rmosek)
  n   <- nrow(L)
  k   <- ncol(L)
  res <- KWDual(L,rep(1,k),rep(1,n)/n)
  print(system.time(out.mosek <- mosek(res$P,
                      opts = list(getinfo = TRUE,verbose = 0))))
  print(out.mosek$dinfo$OPTIMIZER_TIME)
  x = res$f/sum(res$f)"
  @rget x
  return x
end

out = REBayes(L);
