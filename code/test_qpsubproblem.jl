# Short script to test the implementations of the SQP algorithm in
# which the QP subproblem is solved using an active-set method or
# an interior-point method (MOSEK).
using Distributions
using Mosek
using JuMP
include("datasim.jl");
include("likelihood.jl");
include("QPsubprob.jl");
include("mixSQP.jl");

# Initialize the pseudorandom number generator.
srand(1);

# Generate a data set.
@printf "Creating data set.\n"
n = round(Int,5e4);
x = normtmixdatasim(n);

# Compute the likelihood matrix.
@printf "Computing likelihood matrix.\n"
sd = autoselectmixsd(x,nv = 20);
L  = normlikmatrix(x,sd = sd);

# Fit the mixture model.
@printf "Fitting mixture model.\n"
out1 = QPsubprob(L,method = "activeset");
out2 = QPsubprob(L,method = "mosek2");
out3 = QPsubprob(L,method = "mosek");

