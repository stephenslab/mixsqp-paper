# Short script to test the projected gradient method against other
# algorithms.
using Distributions
using LowRankApprox
include("datasim.jl");
include("likelihood.jl");
include("mixEM.jl");
include("mixGD.jl");
include("mixSQP.jl");

# Initialize the sequence of pseudorandom numbers.
srand(1);

# Generate a data set with n = 50,000.
z = normtmixdatasim(round(Int,5e4));

# Compute the 50,000 x 20 likelihood matrix.
sd = autoselectmixsd(z,nv = 20);
L  = normlikmatrix(z,sd = sd);

# Run the mix-SQP algorithm.
outsqp = mixSQP(L,lowrank = "none",eps = 1e-6);

# Run the EM algorithm.
xem, fem = mixEM(L,maxiter = 1000,tol = 1e-4);

# Run the projected gradient descent method.
xpgd, fpgd = mixGD(L,maxiter = 1000);

# Compare the quality of the solutions.
@printf "Objective at SQP solution: %0.12e\n" mixobjective(L,outsqp["x"])
@printf "Objective at PGD solution: %0.12e\n" mixobjective(L,xpgd)
@printf "Objective at EM solution:  %0.12e\n" mixobjective(L,xem)
