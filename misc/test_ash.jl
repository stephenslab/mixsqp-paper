# Short script to test implementation of the adaptive shrinkage
# ("ash") function.
using Distributions
using LowRankApprox
using RCall
include("../code/datasim.jl");
include("../code/likelihood.jl");
include("../code/mixSQP.jl");
include("../code/REBayes.jl");
include("../code/ash.jl");

# Script parameters.
gridmult = 1.1;
n        = round(Int,1e5);

# Initialize the pseudorandom number generator.
srand(1);

# Generate data set.
x = normtmixdatasim(n);
s = ones(n);

# Run adaptive shrinkage with mix-SQP algorithm.
out_mixsqp = ash(x,s,gridmult = gridmult,method = "mixSQP",lowrank = "qr");

# Run adaptive shrinkage with interior point method (KWDual/MOSEK).
out_rebayes = ash(x,s,gridmult = gridmult,method = "REBayes");
