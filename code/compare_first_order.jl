# This is the script used to generate the results comparing mix-SQP
# vs. the EM algorithm. The results for the projected gradient method
# are generated separately in run_pg.m after first running this
# script. To run on the RCC cluster, I set up my environment with the
# following commands:
#
# sinteractive --partition=broadwl --mem=8G
# module load julia/0.6.2
#
# Note that we ran this script with m = 20 and m = 200.
#
n = 20000;
m = 2000;
matrixfile      = "simdata-n=20000-m=2000.csv";
outfile_mixsqp1 = "mixsqp-exact-n=20000-m=2000.csv";
outfile_mixsqp2 = "mixsqp-approx-n=20000-m=2000.csv";
outfile_em      = "em-n=20000-m=2000.csv";

using Distributions
using LowRankApprox
include("datasim.jl");
include("likelihood.jl");
include("mixEM.jl");
include("mixGD.jl");
include("mixSQP.jl");

# Generate the matrix.
srand(1);
@printf "Generating %d x %d data matrix.\n" n m
z  = normtmixdatasim(n);
sd = autoselectmixsd(z,nv = m);
L  = normlikmatrix(z,sd = sd);

# Save the matrix to a CSV file.
@printf "Writing data to CSV file.\n"
writecsv(matrixfile,L);

# Run the methods for a small number of iterations to precompile the
# code.
@printf "Precompiling EM and mixSQP code.\n"
xem, fem, tem = mixEM(L,maxiter = 10);
outsqp1 = mixSQP(L,lowrank = "none",maxiter = 4,maxqpiter = 4,verbose = false);
outsqp2 = mixSQP(L,lowrank = "qr",maxiter = 4,maxqpiter = 4,verbose = false);

# Run mix-SQP with no approximation to the input matrix.
@printf "Fitting model using mix-SQP with exact L.\n"
@time outsqp1 = mixSQP(L,lowrank = "none",eps = 1e-8,sptol = 1e-6,
                       maxqpiter = m,maxiter = 200,verbose = false);
outsqp1["obj"] = outsqp1["obj"]/n;
writecsv(outfile_mixsqp1,[outsqp1["obj"] outsqp1["timing"]]);

# Run the mix-SQP with a low rank (truncated QR) approximation to the
# input matrix.
@printf "Fitting model using mix-SQP with approximate L.\n"
@time outsqp2 = mixSQP(L,lowrank = "qr",eps = 1e-8,sptol = 1e-6,maxqpiter = m,
                 maxiter = 200,verbose = false);
outsqp2["obj"] = outsqp2["obj"]/n;
writecsv(outfile_mixsqp2,[outsqp2["obj"] outsqp2["timing"]]);

# Run the EM algorithm.
@printf "Fitting model using EM.\n"
@time xem, fem, dem, tem = mixEM(L,maxiter = 1000,tol = 1e-6);
fem = fem/n;
writecsv(outfile_em,[fem tem]);

# Compare the quality of the solutions. Note that we need to divide
# the objective by n to match the objective used in the manuscript.
@printf "Objective at SQP1 sol.: %0.12f\n" mixobjective(L,outsqp1["x"],1e-8)/n
@printf "Objective at SQP2 sol.: %0.12f\n" mixobjective(L,outsqp2["x"],1e-8)/n
@printf "Objective at EM sol.:   %0.12f\n" mixobjective(L,xem,1e-8)/n

