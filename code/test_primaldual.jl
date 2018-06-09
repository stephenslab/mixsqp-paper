# TO DO: Explain here what this script is for.
using Distributions
using Mosek
using JuMP
include("datasim.jl");@time x_nonnegsqp = nonnegSQP(L);
include("likelihood.jl");
include("primaldual.jl");

# Initialize the pseudorandom number generator.
srand(1);

# Generate a data set.
n = round(Int,5e4);
x = normtmixdatasim(n);

# Compute the likelihood matrix.
sd = autoselectmixsd(x,nv = 20);
L  = normlikmatrix(x,sd = sd);

# Fit the mixture model using the different formulations
# (simplex-constrained, non-negatively-constrained and dual) and
# different solvers (IP and SQP).
@time x_simplexip  = simplexIP(L);
@time x_nonnegip   = nonnegIP(L);
@time x_dual       = dualIP(L);
@time x_simplexsqp = sqp_simp(L);

@time x_nonnegsqp = nonnegSQP(L);

full([x_simplexip x_nonnegip x_dual x_simplexsqp])
