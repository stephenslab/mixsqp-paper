# Short script to test that the 2 x 3 = 6 different solvers (compared
# in Fig. 1 of the manuscript) work on a medium-size simulated data
# set.
using Distributions
using Mosek
using JuMP
include("datasim.jl");
include("likelihood.jl");
include("primaldual.jl");
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

# Fit the mixture model using the different formulations
# (simplex-constrained, non-negatively-constrained and dual) and
# different solvers (IP and SQP).
@printf "Fitting model...\n"
@printf "simplex, IP:  "
@time x_simplexip  = simplexIP(L);
@printf "non-neg, IP:  "
@time x_nonnegip   = nonnegIP(L);
@printf "dual,    IP:  "
@time x_dualip     = dualIP(L);
@printf "simplex, SQP: "
@time x_simplexsqp = simplexSQP(L);
# @time x_nonnegsqp = nonnegSQP(L);
@time x_dualsqp = sqp_dual(L);

# Check quality of the solutions.
f_simplexip  = mixobjective(L,x_simplexip);
f_nonnegip   = mixobjective(L,x_simplexip);
f_dualip     = mixobjective(L,x_dualip);
f_simplexsqp = mixobjective(L,x_simplexsqp);
f_best       = minimum([f_simplexip f_nonnegip f_dualip f_simplexsqp]);
@printf "Difference between objective at given solution, and best solution:\n"
@printf "simplex, IP:  %0.2e\n" f_simplexip - f_best
@printf "non-neg, IP:  %0.2e\n" f_nonnegip - f_best
@printf "dual,    IP:  %0.2e\n" f_dualip - f_best
@printf "simplex, SQP: %0.2e\n" f_simplexsqp - f_best

