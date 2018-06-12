# Short script to test the implementations of the SQP algorithm in
# which the QP subproblem is solved using an active-set method or
# an interior-point method (MOSEK).
using Distributions
using Mosek
using JuMP
include("datasim.jl");
include("likelihood.jl");
include("mixSQP.jl");

# Initialize the pseudorandom number generator.
srand(1);

# Generate a data set.
@printf "Creating data set.\n"
n = round(Int,5e4);
z = normtmixdatasim(n);

# Compute the likelihood matrix.
@printf "Computing likelihood matrix.\n"
sd = autoselectmixsd(z,nv = 20);
L  = normlikmatrix(z,sd = sd);

# Fit the mixture model.
@printf "Fitting mixture model using active-set method.\n"
out1 = mixSQP(L,qpsubprob = "activeset",lowrank = "none",sptol = 1e-6,
              verbose = false);
@printf "Fitting mixture model using MOSEK.\n"
out2 = mixSQP(L,qpsubprob = "mosek",lowrank = "none",sptol = 1e-6,
              verbose = false);

# Check quality of the solutions.
f1    = mixobjective(L,out1["x"]);
f2    = mixobjective(L,out2["x"]);
fbest = minimum([f1 f2]);
@printf "Difference between objective at given solution, and best solution:\n"
@printf "active-set: %0.2e\n" f1 - fbest
@printf "MOSEK:      %0.2e\n" f2 - fbest

