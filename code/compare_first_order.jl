# Generate the data sets used in the comparisons of mix-SQP vs. the
# first-order methods (EM and projected gradient).
n = 20000;
m = 200;

using Distributions
using LowRankApprox
include("datasim.jl");
include("likelihood.jl");
include("mixEM.jl");
include("mixGD.jl");
include("mixSQP.jl");

# Generate the matrix.
srand(2019);
z  = normtmixdatasim(n);
sd = autoselectmixsd(z,nv = m);
L  = normlikmatrix(z,sd = sd);

# Run the methods.
xem, fem, tem = mixEM(L,maxiter = 10000,tol = 1e-6);
outsqp1 = mixSQP(L,lowrank = "none", eps = 1e-8, maxqpiter = 200,
                 maxiter = 200, verbose = false, nullprior = 0);
outsqp2 = mixSQP(L,lowrank = "qr", eps = 1e-8, maxqpiter = 200,
                 maxiter = 200, verbose = false, nullprior = 0);

# Compare the quality of the solutions.
@printf "Objective at SQP1 solution: %0.12e\n" mixobjective(L,outsqp1["x"])
@printf "Objective at SQP2 solution: %0.12e\n" mixobjective(L,outsqp2["x"])
@printf "Objective at EM   solution: %0.12e\n" mixobjective(L,xem)

# Save the matrix to a CSV file.
writecsv("simdata-n=20000-m=200.csv",L);
