# Small illustration of convergence issue in mixSQP algorithm when QR
# approximation to likelihood matrix yields negative entries in
# likelihood matrix.

# Setup.
using Distributions
using LowRankApprox
include("../code/julia/datasim.jl");
include("../code/julia/likelihood.jl");
include("../code/julia/mixSQP.jl");
include("../code/julia/REBayes.jl")

# Initialize the sequence of pseudorandom numbers.
srand(1);

# Create a data set with 50,000 data points.
x = normtmixdatasim(round(Int,5e4));

# Compute the likelihood matrix.
sd = autoselectmixsd(x,gridmult = 1.425);
L  = normlikmatrix(x,sd = sd);
size(L)

# Compute the solution using the SQP algorithm.
out = mixSQP(L,lowrank = "nothing");
x   = out["x"];

# Compute an approximate solution using the SQP algorithm.
srand(1);
out = mixSQP(L,pqrtol = 1e-10);
x2  = out["x"];

# Compute an approximate solution using the SQP algorithm. If we let
# the eps argument be the default setting (1e-8), the solver will not
# converge.
srand(1); 
out = mixSQP(L,pqrtol = 1e-6,eps = 2e-6);
x3  = out["x"];
