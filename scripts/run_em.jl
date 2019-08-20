# This is the script used to generate the results comparing mix-SQP
# vs. the EM algorithm. The results for the projected gradient method
# are generated separately with run_pg.m.
matfile         = "simdata-n=20000-m=20.csv"
outfile_em      = "em-n=20000-m=20.csv";
outfile_mixsqp1 = "mixsqp-exact-n=20000-m=20.csv";
outfile_mixsqp2 = "mixsqp-approx-n=20000-m=20.csv";

using Printf
using DelimitedFiles
using LinearAlgebra
using SparseArrays
using LowRankApprox
using Statistics
include("../code/mixEM.jl");
include("../code/mixSQP.jl");

# Load data..
@printf("Loading data.\n")
L = readdlm(matfile,',');
n = size(L,1);
m = size(L,2);
@printf "Loaded %d x %d matrix.\n" n m

# Run the methods for a small number of iterations to precompile the
# code.
@printf "Precompiling EM and mixSQP code.\n"
xem, fem, tem = mixEM(L,maxiter = 10);
outsqp1 = mixSQP(L,lowrank = "none",maxiter = 4,maxqpiter = 4,verbose = false);
outsqp2 = mixSQP(L,lowrank = "qr",maxiter = 4,maxqpiter = 4,verbose = false);

# Run the EM algorithm.
@printf "Fitting model using EM.\n"
@time xem, fem, dem, tem = mixEM(L,maxiter = 1000,tol = 1e-6);
fem = fem/n;
writedlm(outfile_em,[fem tem],',');

# Run mix-SQP with no approximation to the input matrix.
@printf "Fitting model using mix-SQP with exact L.\n"
@time outsqp1 = mixSQP(L,lowrank = "none",eps = 1e-8,sptol = 1e-6,
                       maxqpiter = m,maxiter = 200,verbose = false);
outsqp1["obj"] = outsqp1["obj"]/n;
writedlm(outfile_mixsqp1,[outsqp1["obj"] outsqp1["timing"]],',');

# Run the mix-SQP with a low rank (truncated QR) approximation to the
# input matrix.
@printf "Fitting model using mix-SQP with approximate L.\n"
@time outsqp2 = mixSQP(L,lowrank = "qr",eps = 1e-8,sptol = 1e-6,
                       maxqpiter = m,maxiter = 200,verbose = false);
outsqp2["obj"] = outsqp2["obj"]/n;
writedlm(outfile_mixsqp2,[outsqp2["obj"] outsqp2["timing"]],',');

# Compare the quality of the solutions. Note that we need to divide
# the objective by n to match the objective used in the manuscript.
@printf "Objective at EM sol.:   %0.12f\n" mixobjective(L,xem,1e-8)/n
@printf "Objective at SQP1 sol.: %0.12f\n" mixobjective(L,outsqp1["x"],1e-8)/n
@printf "Objective at SQP2 sol.: %0.12f\n" mixobjective(L,outsqp2["x"],1e-8)/n
