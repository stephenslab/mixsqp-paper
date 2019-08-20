# Generate the data sets used in the comparisons of mix-SQP and the
# first-order methods (EM and projected gradient descent).

# SCRIPT PARAMETERS
# -----------------
n        = 20000;
m        = 20;
filename = "simdata-n=20000-m=20.csv";

# SET UP ENVIRONMENT
# ------------------
using Random
using Distributions
using DelimitedFiles
include("../code/datasim.jl");
include("../code/likelihood.jl");

# Initialize the sequence of pseudorandom numbers.
Random.seed!(1);

# SIMULATE DATA
# -------------
# Generate the n x m conditional likelihood matrix.
z  = normtmixdatasim(n);
sd = autoselectmixsd(z,nv = m);
L  = normlikmatrix(z,sd = sd);

# Save the matrix to a CSV file.
writedlm(filename,L,',');
